#' Seizure prediction models
#'
#' Create and run potential seizure prediction models
#' For all neontates OR HIE only

library(librarian)
librarian::shelf(ggplot2, dplyr, epiDisplay, scales, lubridate,
                 e1071, caret, readr, rpart,
                 rpart.plot, rsample, ipred,
                 randomForest, ranger, Rcpp, h2o, ROCR, bit64)

## Source dataset from initial matrix creation script
message("If using an already created dataset, load the desired file as 'matrix_new' and skip to splitting train/test data.")

source("neonate_seizure_predict_matrix")
source("../utils.R")

print("Source data loaded, creating train-test matrices...")

### Create train/test matrices for analyses with various models ----
# This includes removing missing data
# And building a table to store model results

matrix_new <- day_1_matrix %>%
  mutate(across(
  .cols = everything(),
  .fns = ~factor(.),
  .names = NULL)) # Converting all numeric columns into binary categorical variables

matrix_new <- matrix_new %>%
  dplyr::select(-c("hie_th", "events"))

matrix_new <- na.omit(matrix_new,
                      cols = c("sex", "race", "graphoelements", "voltage", "continuity",
                               "transients", "variability", "reactivity", "impression",
                               "seizures", "subday_sz"))

matrix_new <- matrix_new %>%
  dplyr::select(-c("END_EXAM_DTTM", "month", "year")) # Remove for now, only needed for iterative model analyses

## Split data
set.seed(123)
matrix_idx <- createDataPartition(matrix_new$subday_sz, p = 0.80, list = FALSE)
matrix_trn <- matrix_new[matrix_idx, ]
matrix_tst <- matrix_new[-matrix_idx, ]

print("Build non-H2o Models...")

## Build table to store any categorical model results for later review
results_tab <- setNames(data.frame(matrix(ncol = 7, nrow = 1)),
                        c("accuracy", "precision", "recall", "F1", "AUC", "AUCPR", "kappa"))

warning("Ensure write.csvs are NOT commented out where wanted below in order to save results")

### Logistic regression model w/ Caret ----
print("Building logistic regression model...")

set.seed(123)
matrix_glm_mod <- train(
  form = factor(subday_sz) ~ ., # Convert to factor for binomial classification
  data = matrix_trn,
  trControl = trainControl(method = "cv", number = 10),
  method = "glm",
  family = "binomial"
)
matrix_glm_mod
matrix_glm_mod$results
matrix_glm_mod$finalModel
summary(matrix_glm_mod)

print("...logistic regression model built...")

# Predict test data
head(predict(matrix_glm_mod, newdata = matrix_tst))

# Get probabilities
head(predict(matrix_glm_mod, newdata = matrix_tst, type = "prob"))

## Test model for accuracy using test dataset
calc_acc(actual = matrix_tst$subday_sz,
         predicted = predict(matrix_glm_mod,
                             newdata = matrix_tst))

pred_raw <- predict(matrix_glm_mod, newdata = matrix_tst, type = "raw")
cfm <- table(matrix_tst$subday_sz, (pred_raw == "1") * 1, dnn = c("Actual", "Predicted"))

## Save model results
precision <- precision(predict(matrix_glm_mod, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(matrix_glm_mod, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(matrix_glm_mod, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(matrix_glm_mod, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
cfm <- confusionMatrix(predict(matrix_glm_mod, matrix_tst), matrix_tst$subday_sz)
accuracy <- cfm$overall["Accuracy"]
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
kappa <- cfm$overall[[2]]
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])

row.names(results_tab) <- c("log_regress_caret")
results_tab$precision <- precision
results_tab$recall <- recall
results_tab$F1 <- f1
results_tab$AUC <- auc
results_tab$AUCPR <- aucpr
results_tab$accuracy <- accuracy
results_tab$kappa <- kappa

print("...logstic regression model completed.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

### Regression tree model using rpart ----

print("Building rpart regression tree...")

set.seed(123)
data_input <- day_1_matrix %>%
  mutate(across(
  .cols = everything(),
  .fns = ~factor(.),
  .names = NULL)) # Converting all numeric columns into binary categorical variables

accuracy_test <- NULL
rpt_models <- list()
for (i in 1:5) { # Number of iterations of train/test split to perform from data

  print(i)

  shuffle_index <- sample(seq_len(NROW(data_input)))

  data_input <- data_input[shuffle_index, ]

  # Clean data
  clean_data_input <- data_input

  # Convert NAs to factor level
  clean_data_input[is.na(clean_data_input)] <- factor(0)

  glimpse(clean_data_input)

  clean_data_input <- clean_data_input %>%
    dplyr::select(seizures, voltage, graphoelements, continuity,
                  variability, transients, impression, reactivity,
                  race, sex, subday_sz)

  # Training/testing data function
  create_train_test <- function(data, size = 0.8, train = TRUE) {
    n_row <- nrow(data)
    total_row <- size * n_row
    train_sample <- 1:total_row
    if (train == TRUE) {
      return(data[train_sample, ])
    } else {
      return(data[-train_sample, ])
    }
  }

  data_train <- create_train_test(clean_data_input, 0.8, train = TRUE)
  data_test <- create_train_test(clean_data_input, 0.8, train = FALSE)
  dim(data_train)

  # Check randomization
  prop.table(table(data_train$subday_sz))
  prop.table(table(data_test$subday_sz))

  # Building the Model
  fit <- rpart(subday_sz ~ .,
               data = data_train,
               method = "class",
               minsplit = 10
               )
  rpart.plot(fit, box.palette = "auto")


  # Testing the Model
  predict_unseen <- predict(fit, data_test, type = "class")

  table_mat <- table(data_test$subday_sz, predict_unseen)
  table_mat

  accuracy_test[i] <- sum(diag(table_mat)) / sum(table_mat)
  print(paste("Accuracy for test: ", accuracy_test[i]))


  # Alternative model creation: grid of models with different split and depths to find best tuned model
  hyper_grid <- expand.grid(
    minsplit = seq(5, 20, 1), # Min data points to attempt a split
    maxdepth = seq(1, 10, 1) # Max nodes between root and terminal
  )
  test <- c(1, 2, 3)
  for (i in seq_len(NROW(hyper_grid))) {

    # Get minsplit, maxdepth values at row i
    minsplit <- hyper_grid$minsplit[i]
    maxdepth <- hyper_grid$maxdepth[i]

    # Train a model and store in the list
    rpt_models[[i]] <- rpart(
      formula = subday_sz ~ .,
      data = data_train,
      method = "class",
      control = list(minsplit = minsplit, maxdepth = maxdepth)
    )
  }

}

hist(accuracy_test)
median(accuracy_test)

test3 <- as.data.frame(accuracy_test)

ggplot(test3,
       aes(x = accuracy_test,
           group = 1)) +
  geom_density()

print("...rpart regression model completed.")


### Alternative Logistic Regression model comparison ----

# Logistic Regression - first day predictors, seizures after day one
set.seed(123)
glm_fit <- glm(factor(subday_sz) ~ seizures + voltage + graphoelements + continuity + variability + transients + impression + reactivity + race + sex, data = clean_data_input,
               family = binomial)
summary(glm_fit)

# Using epiDisplay to get OR
glm1 <- glm(subday_sz ~ seizures + voltage + graphoelements + continuity + variability + transients + impression + reactivity + race + sex,
            data = clean_data_input,
            family = binomial)

logistic.display(glm1)

### Bagged regression tree using caret ----

print("Building caret regression tree model...")

set.seed(123)

## First, run a bagged regression tree (slightly different from random forest)
matrix_tree_mod <- train(
  form = subday_sz ~ ., # Convert to factor for binomial classification
  data = matrix_trn,
  trControl = trainControl(method = "cv", number = 10),
  method = "treebag"
)
matrix_tree_mod

print("...regression tree model built...")

matrix_tree_mod$results
matrix_tree_mod$finalModel
summary(matrix_tree_mod)

## Predict test data
head(predict(matrix_tree_mod, newdata = matrix_tst))

## Get probabilities
head(predict(matrix_tree_mod, newdata = matrix_tst, type = "prob"))

## Test for model accuracy using test dataset
calc_acc(actual = matrix_tst$subday_sz,
         predicted = predict(matrix_tree_mod,
                             newdata = matrix_tst))

## Save model results
precision <- precision(predict(matrix_tree_mod, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(matrix_tree_mod, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(matrix_tree_mod, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(matrix_tree_mod, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
accuracy <- confusionMatrix(predict(matrix_tree_mod, matrix_tst), matrix_tst$subday_sz)$overall["Accuracy"]
kappa <- confusionMatrix(predict(matrix_tree_mod, matrix_tst), matrix_tst$subday_sz)$overall["Kappa"]
cfm <- confusionMatrix(predict(matrix_tree_mod, matrix_tst), matrix_tst$subday_sz)
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])

results_tab["regresstree_caret", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[2] <- precision
results_tab$recall[2] <- recall
results_tab$F1[2] <- f1
results_tab$AUC[2] <- auc
results_tab$AUCPR[2] <- aucpr
results_tab$accuracy[2] <- accuracy
results_tab$kappa[2] <- kappa

print("...Caret regression tree model completed.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

### Traditional random forest models (i.e. random forest package) ----

print("Building random forest model...")

set.seed(123)

model_1 <- randomForest(
  formula = subday_sz ~ .,
  data = matrix_trn
)

model_1
plot(model_1)
which.min(model_1$mse)

print("...random forest model built...")

## Save model results
precision <- precision(predict(model_1, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(model_1, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(model_1, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(model_1, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
accuracy <- confusionMatrix(predict(model_1, matrix_tst), matrix_tst$subday_sz)$overall["Accuracy"]
kappa <- confusionMatrix(predict(model_1, matrix_tst), matrix_tst$subday_sz)$overall["Kappa"]
cfm <- confusionMatrix(predict(model_1, matrix_tst), matrix_tst$subday_sz)
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])


results_tab["random_model1", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[3] <- precision
results_tab$recall[3] <- recall
results_tab$F1[3] <- f1
results_tab$AUC[3] <- auc
results_tab$AUCPR[3] <- aucpr
results_tab$accuracy[3] <- accuracy
results_tab$kappa[3] <- kappa

print("Random forest initial model complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

## Run again but with validation datasets (basically w/ built in predictor tests)
# This one won't be saved
set.seed(123)

## Create training data
valid_split <- initial_split(matrix_trn, .8)
train_v2 <- analysis(valid_split) # 80% of the data from the trn matrix

## Create validation data and run model
matrix_valid <- assessment(valid_split) # The remaining 20% from the trn matrix
x_test <- matrix_valid[setdiff(names(matrix_valid), "subday_sz")]
y_test <- matrix_valid$subday_sz

rf_oob_comp <- randomForest(
  formula = subday_sz ~ .,
  data    = train_v2,
  xtest   = x_test,
  ytest   = y_test
)

rf_oob_comp$err.rate
# Extract OOB & test set errors for comparison (note: will not generation mse with data in certain class types)
oob <- rf_oob_comp$err.rate
validation <- rf_oob_comp$err.rate
rf_oob_comp$confusion
tail(rf_oob_comp$err.rate)

tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree) %>%
  gather(Metric, RMSE, -ntrees) %>%
  ggplot(aes(ntrees, RMSE, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

## Random forest model again but with tuning
print("Building random forest with mtry tuning...")

set.seed(123)

features <- setdiff(names(matrix_trn), "subday_sz")

# Just adjusting mtry (i.e. number of features sample at each split)
set.seed(123)
m2 <- tuneRF(
  x = matrix_trn[features],
  y = matrix_trn$subday_sz,
  ntreeTry = 500,
  mtryStart = 10,
  stepFactor = 1.5,
  improve = 0.01,
  trace = FALSE
)
as.data.frame(m2) %>% arrange(OOBError)

set.seed(123)
model_2 <- randomForest(
  formula = subday_sz ~ .,
  data = matrix_trn,
  mtry = 7
)

model_2$predicted

print("...random forest model with mtry tuning built...")

## Save model results
precision <- precision(predict(model_2, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(model_2, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(model_2, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(model_2, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
accuracy <- confusionMatrix(predict(model_2, matrix_tst), matrix_tst$subday_sz)$overall["Accuracy"]
kappa <- confusionMatrix(predict(model_2, matrix_tst), matrix_tst$subday_sz)$overall["Kappa"]
cfm <- confusionMatrix(predict(model_2, matrix_tst), matrix_tst$subday_sz)
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])


results_tab["random_model_mtry", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[4] <- precision
results_tab$recall[4] <- recall
results_tab$F1[4] <- f1
results_tab$AUC[4] <- auc
results_tab$AUCPR[4] <- aucpr
results_tab$accuracy[4] <- accuracy
results_tab$kappa[4] <- kappa

print("...mtry random forest model complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

## Another model with broader random forest tuning (not saving this iterative test)
print("Building random forest model with additional hyperparameters...")

hyper_grid <- expand.grid(
  mtry = seq(1, 10, by = 1),
  node_size = seq(3, 9, by = 2),
  ntree = seq(250, 500, by = 50),
  oob_rmse = 0
)
NROW(hyper_grid)

## Test for smallest RMSE or error rate
set.seed(123)
for (i in seq_len(NROW(hyper_grid))){

  model <- randomForest(
    formula = subday_sz ~ .,
    data = matrix_trn,
    mtry = hyper_grid$mtry[i],
    nodesize = hyper_grid$node_size[i],
    #sampsize = hyper_grid$sampe_size[i], # Remove for categorical if imbalanced and go with default
    ntree = hyper_grid$ntree[i]
  )
  hyper_grid$oob_rmse[i] <- model$err.rate # sqrt of MSE for continuous, err.rate for categorical
}
hyper_grid %>% arrange(oob_rmse) %>% head(10)

## Random Forest, now with "optimal" settings
set.seed(123)
model_3 <- randomForest(
  formula = subday_sz ~ .,
  data = matrix_trn,
  mtry = 4,
  nodesize = 9,
  #sampsize = .55, # Again, go with default if categorical and imbalanced dataset
  ntree = 500)

print("...random forest model with additional hyperparameters built...")

model_3$predicted
head(model_3$err.rate)
#sqrt(model_3$mse[which.min(model_3$mse)]) # Use for numeric rather than categorical data
model_3$confusion

## Save model results
precision <- precision(predict(model_3, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(model_3, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(model_3, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(model_3, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
accuracy <- confusionMatrix(predict(model_3, matrix_tst), matrix_tst$subday_sz)$overall["Accuracy"]
kappa <- confusionMatrix(predict(model_3, matrix_tst), matrix_tst$subday_sz)$overall["Kappa"]
cfm <- confusionMatrix(predict(model_3, matrix_tst), matrix_tst$subday_sz)
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])


results_tab["random_model_opt", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[5] <- precision
results_tab$recall[5] <- recall
results_tab$F1[5] <- f1
results_tab$AUC[5] <- auc
results_tab$AUCPR[5] <- aucpr
results_tab$accuracy[5] <- accuracy
results_tab$kappa[5] <- kappa

print("...Random forest model with hyperparameter tuning complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

### RANGER Random Forest ----

print("Building ranger random forest model...")

# Ranger provides more efficient scaling/testing of features
set.seed(123)

system.time(
  matrix_ranger <- ranger(
    formula = subday_sz ~ .,
    data = matrix_trn,
    num.trees = 500,
    mtry = 3
  )
)

## Create and test potential tuning hyperparameters for forest
hyper_grid <- expand.grid(
  mtry = seq(1, 10, by = 1),
  node_size = seq(3, 9, by = 2),
  sampe_size = c(.55, .632, .70, .80),
  ntree = seq(250, 500, by = 50),
  oob_rmse = 0
)
nrow(hyper_grid)

for (i in seq_len(NROW(hyper_grid))) {

  model <- ranger(
    formula = subday_sz ~ .,
    data = matrix_trn,
    num.trees = hyper_grid$ntree[i],
    mtry = hyper_grid$mtry[i],
    min.node.size = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sampe_size[i], # Unlike random forest function parameters, ranger seems to run this even imbalanced
    seed = 123
  )

  hyper_grid$oob_rmse[i] <- sqrt(model$prediction.error)

}
hyper_grid %>% arrange(oob_rmse) %>% head(10)

## Build actual ranger model based on "optimal" parameters
range_1 <- ranger(
  formula = subday_sz ~ .,
  data = matrix_trn,
  num.trees = 250,
  mtry = 8,
  min.node.size = 3,
  sample.fraction = .632,
  seed = 123
)
range_1

print("...ranger random forest model built...")

## Save model results
precision <- precision(predict(range_1, matrix_tst)$predictions, matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(range_1, matrix_tst)$predictions, matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(range_1, matrix_tst)$predictions, matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(range_1, matrix_tst)$predictions), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
accuracy <- confusionMatrix(predict(range_1, matrix_tst)$predictions, matrix_tst$subday_sz)$overall["Accuracy"]
kappa <- confusionMatrix(predict(range_1, matrix_tst)$predictions, matrix_tst$subday_sz)$overall["Kappa"]
cfm <- confusionMatrix(predict(range_1, matrix_tst)$predictions, matrix_tst$subday_sz)
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])


results_tab["range_1", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[6] <- precision
results_tab$recall[6] <- recall
results_tab$F1[6] <- f1
results_tab$AUC[6] <- auc
results_tab$AUCPR[6] <- aucpr
results_tab$accuracy[6] <- accuracy
results_tab$kappa[6] <- kappa

print("Ranger random forest model complete...")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

## Alternative attempt at one hot encoding of categorical variables for Ranger
set.seed(123)

temp_matrix <- day_1_matrix %>%
  dplyr::select(-c("END_EXAM_DTTM", "month", "year", "hie_th", "events"))

temp_matrix$subday_sz <- as.numeric(as.character(temp_matrix$subday_sz)) # Convert to numeric to avoid dummyvars splitting into dummy prediction variables

one_hot <- dummyVars(~., temp_matrix, fullRank = FALSE)

matrix_train_hot <- predict(one_hot,
                            newdata = temp_matrix) %>% as.data.frame()

matrix_train_hot <- na.omit(matrix_train_hot)

matrix_idx_hot <- createDataPartition(matrix_train_hot$subday_sz, p = 0.80, list = FALSE)
matrix_trn_hot <- matrix_train_hot[matrix_idx, ]
matrix_tst_hot <- matrix_train_hot[-matrix_idx, ]

## Make ranger compatible names and build hyperparameter grid
names(matrix_trn_hot) <- make.names(names(matrix_trn_hot), allow_ = FALSE)

hyper_grid2 <- expand.grid(
  mtry = seq(1, 10, by = 1),
  node_size = seq(3, 9, by = 2),
  sampe_size = c(.55, .632, .70, .80),
  ntree = seq(250, 500, by = 50),
  oob_rmse = 0
)
NROW(hyper_grid2)

for (i in seq_len(NROW(hyper_grid2))) {

  model <- ranger(
    formula = subday_sz ~ .,
    data = matrix_trn_hot,
    num.trees = 500,
    mtry = hyper_grid2$mtry[i],
    min.node.size = hyper_grid2$node_size[i],
    sample.fraction = hyper_grid2$sampe_size[i],
    seed = 123
  )

  hyper_grid2$oob_rmse[i] <- sqrt(model$prediction.error)

}
hyper_grid2 %>% arrange(oob_rmse) %>% head(10)

## Use optimal settings in ranger to help build model
oob_rmse <- vector(mode = "numeric", length = 100)

for (i in seq_along(oob_rmse)){
  opt_range_1 <- ranger(
    formula = subday_sz ~ .,
    data = matrix_trn_hot,
    num.trees = 250,
    mtry = 3,
    min.node.size = 3,
    sample.fraction = .55,
    importance = "impurity"
  )
  oob_rmse[i] <- sqrt(opt_range_1$prediction.error)
}
hist(oob_rmse, breaks = 20)

opt_range_1$variable.importance %>%
  tidy() %>%
  arrange(desc(x)) %>%
  top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables")

print("Non-H2o Models built.")

### H2o Random Forest Models - Set up and initial model ----

print("Building initial h2o random forest model...")

## In case version is outdated, use installation process below:
# if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
# if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
# pkgs <- c("RCurl","jsonlite")
# for (pkg in pkgs) {
#   if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
# }
# install.packages("h2o"), type = "source", repos = (c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))
# library(h2o)


## Initialize h2o model building and hyperparameter grid for model tuning
set.seed(123)

h2o.init(max_mem_size = "10g")

# In case h2o is getting buggy, try shutting down or performing cleanup loop.
# Be sure to use removeAll() if you plan to re-run/recreate dataframes

h2o.shutdown()
h2o.removeAll()

y <- "subday_sz"
x <- setdiff(names(matrix_trn), y)

# Convert dataframes to h2o object
train_h2o <- as.h2o(matrix_trn)
test_h2o <- as.h2o(matrix_tst)

# Hyperparameters for any h2o grid searches
# Note this first grid is fairly exhaustive and may not be necessary for these analyses

hyper_grid_h2o_full <- list(
  ntrees      = seq(200, 500, by = 100),
  mtries      = seq(1, 10, by = 1),
  max_depth   = seq(1, 40, by = 5),
  min_rows    = seq(1, 5, by = 2),
  nbins       = seq(10, 30, by = 5),
  sample_rate = c(.55, .632, .80)
)

# Balance models can use this one instead of the more exhaustive hyperparameters above
hyper_grid_h2o <- list(
  ntrees      = seq(200, 500, by = 100),
  mtries      = seq(1, 10, by = 1),
  sample_rate = c(.55, .632, 0.70, .80)
)

## Build grid search - note this will take longer than ranger
# Grid performance tests can be sped up with a "randomdiscrete" search for optimal performance rather than "Cartesian"
# which will not run all options in the hypergrid, but will get you close.  Not recommended if also using CV
grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid",
  x = x,
  y = y,
  training_frame = train_h2o,
  hyper_params = hyper_grid_h2o,
  search_criteria = list(strategy = "Cartesian"),
  seed = 123
)

grid_perf <- h2o.getGrid(
  grid_id = "rf_grid",
  sort_by = "aucpr",
  decreasing = TRUE # Switch to false if sorting by MSE
)
print(grid_perf)

# Example of test view of "inflection point" for a hyperparamter to see when increases stop "helping" model
ggplot(as.data.frame(sapply(grid_perf@summary_table, as.numeric))) +
  geom_point(aes(mtries, aucpr)) +
  geom_line(aes(mtries, aucpr, group = 1)) +
  labs(x = "Features", y = "aucpr", title = "Grid Search for Single Tree Models")

# Save "best" model, given criteria above (note this uses the full Cartesian search grid from higher up, NOT the random discrete one)
h2o_model_id <- grid_perf@model_ids[[1]]
h2o_model <- h2o.getModel(h2o_model_id)
model_path <- h2o.saveModel(object = h2o_model,
                            path = "your_path_here",
                            filename = "rf_h2o_1",
                            force = TRUE)

# Run model against test set
h2o_model_perf <- h2o.performance(model = h2o_model,
                                  newdata = test_h2o)

# RMSE of best model
h2o.mse(h2o_model_perf) %>% sqrt()

## Save model results
# Note: for h2o, based all metric thresholds based off of F1 score maximization, unlike caret which tries to average everything
h2o.F1(h2o_model_perf) %>% arrange(desc(f1)) # To find max f1 threshold
h2o.confusionMatrix(h2o_model_perf) # Use confusion matrix to find max f1 threshold more accurately/quickly
precision <- h2o.precision(h2o_model_perf, thresholds = model_threshold(h2o_model_perf))[[1]]
recall <- h2o.recall(h2o_model_perf, thresholds = model_threshold(h2o_model_perf))[[1]]
f1 <- h2o.F1(h2o_model_perf, thresholds = model_threshold(h2o_model_perf))[[1]]
auc <- h2o.auc(h2o_model_perf)
aucpr <- h2o.aucpr(h2o_model_perf)
accuracy <- h2o.accuracy(h2o_model_perf, thresholds = model_threshold(h2o_model_perf))[[1]]
kappa <- cohen_kappa(h2o_model_perf)


results_tab["h2o_1", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[7] <- precision
results_tab$recall[7] <- recall
results_tab$F1[7] <- f1
results_tab$AUC[7] <- auc
results_tab$AUCPR[7] <- aucpr
results_tab$accuracy[7] <- accuracy
results_tab$kappa[7] <- kappa

print("...Initial h2o model complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

### H2o balanced and stratified models ----

print("Building h2o models with balance and stratification...")

## Build grid search - this time balanced, stratified and with cross fold validation of training dataset
# Note: test using random search to speed up and not run literally all models (but don't use in conjunction with CV)

# Random grid search criteria
rand_search_criteria2 <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "logloss",
  stopping_tolerance = 0.001, # 0.001 is the default
  stopping_rounds = 5,
  max_runtime_secs = 600,
  seed = 123
)

grid_b <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid_balance",
  balance_classes = TRUE,
  fold_assignment = "Stratified",
  nfolds = 10,
  keep_cross_validation_fold_assignment = FALSE,
  keep_cross_validation_predictions = FALSE,
  x = x,
  y = y,
  training_frame = train_h2o,
  hyper_params = hyper_grid_h2o,
  search_criteria = list(strategy = "Cartesian"),
  seed = 123
)

 grid_b_perf <- h2o.getGrid(
  grid_id = "rf_grid_balance",
  sort_by = "aucpr",
  decreasing = TRUE # Switch to false if looking at MSE or logloss
)
print(grid_b_perf)

# Create model from best balanced hyperparameters and view parameters
h2o_balanced <- h2o.getModel(grid_b_perf@model_ids[[1]])

print(h2o_balanced@model[["model_summary"]])

# Save the model
model_path <- h2o.saveModel(object = h2o_balanced,
                            path = "your_path_here",
                            filename = "rf_h2o_balanced",
                            force = TRUE)


# Test balanced model
h2o_bal_best <- h2o.performance(model = h2o_balanced,
                                newdata = test_h2o)


## Save model results
h2o.confusionMatrix(h2o_bal_best) # Use confusion matrix to find max f1 threshold
precision <- h2o.precision(h2o_bal_best, thresholds = model_threshold(h2o_bal_best))[[1]]
recall <- h2o.recall(h2o_bal_best, thresholds = model_threshold(h2o_bal_best))[[1]]
f1 <- h2o.F1(h2o_bal_best, thresholds = model_threshold(h2o_bal_best))[[1]]
auc <- h2o.auc(h2o_bal_best)
aucpr <- h2o.aucpr(h2o_bal_best)
accuracy <- h2o.accuracy(h2o_bal_best, thresholds = model_threshold(h2o_bal_best))[[1]]
kappa <- cohen_kappa(h2o_bal_best)

results_tab["h2o_balanced", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[8] <- precision
results_tab$recall[8] <- recall
results_tab$F1[8] <- f1
results_tab$AUC[8] <- auc
results_tab$AUCPR[8] <- aucpr
results_tab$accuracy[8] <- accuracy
results_tab$kappa[8] <- kappa

print("...H2o model with balance and stratified complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

### H2o single models with various balancing metrics ----

print("Building alternative balanced h2o with over under sampling...")

# Re-iniating H2o to make additional memory space
h2o.removeAll()

set.seed(123)

# Tune model using h2o

y <- "subday_sz"
x <- setdiff(names(matrix_trn), y)

# convert to h2o object
train_h2o <- as.h2o(matrix_trn)
test_h2o <- as.h2o(matrix_tst)

grid_cusb <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid_cus_balance",
  balance_classes = TRUE,
  fold_assignment = "Stratified",
  class_sampling_factors = c(0.5, .9), # Undersample majority and oversample minority
  nfolds = 10,
  keep_cross_validation_fold_assignment = FALSE,
  keep_cross_validation_predictions = FALSE,
  x = x,
  y = y,
  training_frame = train_h2o,
  hyper_params = hyper_grid_h2o,
  search_criteria = list(strategy = "Cartesian"),
  seed = 123
)

grid_cusb_perf <- h2o.getGrid(
  grid_id = "rf_grid_cus_balance",
  sort_by = "aucpr",
  decreasing = TRUE # Switch to false if looking at MSE or logloss
)
print(grid_cusb_perf)

# Create model from best balanced hyperparameters and view parameters
h2o_cus_balanced <- h2o.getModel(grid_cusb_perf@model_ids[[1]])

print(h2o_cus_balanced@model[["model_summary"]])

# Save the model
model_path <- h2o.saveModel(object = h2o_cus_balanced,
                            path = "your_path_here",
                            filename = "rf_h2o_custom_balanced",
                            force = TRUE)


# Test balanced model
h2o_cust_bal <- h2o.performance(model = h2o_cus_balanced,
                                newdata = test_h2o)


## Save model results
h2o.confusionMatrix(h2o_cust_bal)
precision <- h2o.precision(h2o_cust_bal, thresholds = model_threshold(h2o_cust_bal))[[1]]
recall <- h2o.recall(h2o_cust_bal, thresholds = model_threshold(h2o_cust_bal))[[1]]
f1 <- h2o.F1(h2o_cust_bal, thresholds = model_threshold(h2o_cust_bal))[[1]]
auc <- h2o.auc(h2o_cust_bal)
aucpr <- h2o.aucpr(h2o_cust_bal)
accuracy <- h2o.accuracy(h2o_cust_bal, thresholds = model_threshold(h2o_cust_bal))[[1]]
kappa <- cohen_kappa(h2o_cust_bal)

results_tab["h2o_custom_bal", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[9] <- precision
results_tab$recall[9] <- recall
results_tab$F1[9] <- f1
results_tab$AUC[9] <- auc
results_tab$AUCPR[9] <- aucpr
results_tab$accuracy[9] <- accuracy
results_tab$kappa[9] <- kappa

print("...H2o custom balanced model complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

### H2o weighted models ----

## Build grid search - this time using weighted columns as alternative to balancing, stratified and with cross fold validation of training dataset
print("Start loop of all weighted H2o models...")

h2o.removeAll()

# Set up for for loop and dataframe for various confusion matrices
results_row <- 10
weights_row <- 1
cm_list <- list()

# "Perfectly" weighted, then other tests of weight variation
weight_tbl <- data.frame(w0 = c(0.6068152, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                         w1 = c(2.840491, 1.5, 2.0, 3.0, 4.0, 5.0, 10, 15.0)
                         )

# Iterate through each weight combination and run in the h2o model
for (i in seq_len(NROW(weight_tbl))) {

  print(paste0("starting weighted model: ", weights_row))

  set.seed(123)

  #h2o.init(max_mem_size = "5g")

  y <- "subday_sz"
  x <- setdiff(names(matrix_trn), y)

  # Convert to h2o object
  train_h2o <- as.h2o(matrix_trn)
  test_h2o <- as.h2o(matrix_tst)

  # Create hyperparameters
  hyper_grid_h2o <- list(
    ntrees      = seq(200, 500, by = 100),
    mtries      = seq(1, 10, by = 1),
    sample_rate = c(.55, .632, 0.70, .80)
  )

  # Add a weights column to give each class roughly equal weight or intentionally "imbalanced" weights
  # (basically tells the algorithm how often to repeat a certain row)
  w0 <- weight_tbl$w0[weights_row] # For exact balance...sum[w1+w0] / (2*[w1])
  w1 <-  weight_tbl$w1[weights_row] # For exact balance...sum[w1+w0] / (2*[w0])

  weight_trn <- as.data.frame(train_h2o) %>% mutate(weights =
                                       case_when(subday_sz == "0" ~ w0,
                                                 subday_sz == "1" ~ w1))
  weight_trn <- as.h2o(weight_trn)

  weight_tst <- as.data.frame(test_h2o) %>% mutate(weights =
                                                      case_when(subday_sz == "0" ~ w0,
                                                                subday_sz == "1" ~ w1))
  weight_tst <- as.h2o(weight_tst)


  grid_w <- h2o.grid(
    algorithm = "randomForest",
    grid_id = paste0("rf_grid_weighted_", w0, "_", w1),
    balance_classes = FALSE,
    weights_column = "weights",
    fold_assignment = "Stratified",
    nfolds = 10,
    keep_cross_validation_fold_assignment = FALSE,
    keep_cross_validation_predictions = FALSE,
    x = x,
    y = y,
    training_frame = weight_trn,
    hyper_params = hyper_grid_h2o,
    search_criteria = list(strategy = "Cartesian"),
    seed = 123
  )

  grid_w_perf <- h2o.getGrid(
    grid_id = paste0("rf_grid_weighted_", w0, "_", w1),
    sort_by = "aucpr",
    decreasing = TRUE # Switch to false if looking at MSE or logloss
  )
  print(grid_w_perf)

  # Create model from best balanced hyperparameters and view parameters
  h2o_weighted <- h2o.getModel(grid_w_perf@model_ids[[1]])

  print(h2o_weighted@model[["model_summary"]])

  # Save the model
  model_path <- h2o.saveModel(object = h2o_weighted,
                              path = "your_path_here",
                              force = TRUE)

  # Test balanced model
  h2o_weight_best <- h2o.performance(model = h2o_weighted,
                                  newdata = weight_tst)


  ## Save model results
  warning("Remember to recreate/rerun random forests for each weighted set before adding results")
  cm_list <- append(cm_list, h2o.confusionMatrix(h2o_weight_best))

  precision <- h2o.precision(h2o_weight_best, thresholds = model_threshold(h2o_weight_best))[[1]]
  recall <- h2o.recall(h2o_weight_best, thresholds = model_threshold(h2o_weight_best))[[1]]
  f1 <- h2o.F1(h2o_weight_best, thresholds = model_threshold(h2o_weight_best))[[1]]
  auc <- h2o.auc(h2o_weight_best)
  aucpr <- h2o.aucpr(h2o_weight_best)
  accuracy <- h2o.accuracy(h2o_weight_best, thresholds = model_threshold(h2o_weight_best))[[1]]
  kappa <- cohen_kappa(h2o_weight_best)

  results_tab[paste0("h2o_weighted_", w0, "_", w1), ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
  results_tab$precision[results_row] <- precision
  results_tab$recall[results_row] <- recall
  results_tab$F1[results_row] <- f1
  results_tab$AUC[results_row] <- auc
  results_tab$AUCPR[results_row] <- aucpr
  results_tab$accuracy[results_row] <- accuracy
  results_tab$kappa[results_row] <- kappa

  print(paste0("finishing weighted model: ", weights_row))

  print("saving results...")

  results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
  write.csv(results_tab_save, "results/test_tab.csv")

  print("...results saved.")

  h2o.removeAll()

  results_row <- results_row + 1
  weights_row <- weights_row + 1

}
print("...H2o weighted model iterations complete.")

write.csv(cm_list, "conf_mats_weightedh2o.csv")

## If desired, Check models manually for accuracy
model_path <- "your_model_path_here"
temp <- h2o.loadModel(model_path)
temp@parameters
temp2 <- h2o.performance(temp, newdata = test_h2o)
h2o.confusionMatrix(temp2)
cohen_kappa(temp2)
h2o.F1(temp2, thresholds = model_threshold(temp2))[[1]]
h2o.precision(temp2, thresholds = model_threshold(temp2))[[1]]
h2o.recall(temp2, thresholds = model_threshold(temp2))[[1]]

### Weighted logistic regression ----

print("Building weighted logistic regression models...")

# Added for confidence interval generation, which skips initial results row generation
if (!exists("results_row")) {
  results_row <- 0
}

## Logistic regression model w/ Caret

print("Building balanced logistic regression models...")

weights <- ifelse(matrix_trn$subday_sz == "0",
                  0.6068152,
                  2.840491)

set.seed(123)

matrix_glm_mod_w <- train(
  form = factor(subday_sz) ~ ., # Convert to factor for binomial classification
  data = matrix_trn,
  trControl = trainControl(method = "cv", number = 10),
  method = "glm",
  family = "binomial",
  weights = weights
)
matrix_glm_mod_w
matrix_glm_mod_w$results
matrix_glm_mod_w$finalModel
summary(matrix_glm_mod_w)

print("...balanced logistic regression model built...")

# Predict test data
head(predict(matrix_glm_mod_w, newdata = matrix_tst))

# Get probabilities
head(predict(matrix_glm_mod_w, newdata = matrix_tst, type = "prob"))

# Test model for accuracy using test dataset
calc_acc(actual = matrix_tst$subday_sz,
         predicted = predict(matrix_glm_mod_w, newdata = matrix_tst))

pred_raw <- predict(matrix_glm_mod_w, newdata = matrix_tst, type = "raw")
cfm <- table(matrix_tst$subday_sz, (pred_raw == "1") * 1, dnn = c("Actual", "Predicted"))

## Save model results
precision <- precision(predict(matrix_glm_mod_w, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(matrix_glm_mod_w, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(matrix_glm_mod_w, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(matrix_glm_mod_w, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
cfm <- confusionMatrix(predict(matrix_glm_mod_w, matrix_tst), matrix_tst$subday_sz)
accuracy <- cfm$overall["Accuracy"]
kappa <- cfm$overall["Kappa"]
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])

results_tab["log_regress_caret_wb", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[results_row] <- precision
results_tab$recall[results_row] <- recall
results_tab$F1[results_row] <- f1
results_tab$AUC[results_row] <- auc
results_tab$AUCPR[results_row] <- aucpr
results_tab$accuracy[results_row] <- accuracy
results_tab$kappa[results_row] <- kappa

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

results_row <- results_row + 1

#
print("Building w3 logistic regression model...")

weights <- ifelse(matrix_trn$subday_sz == "0",
                  0.5,
                  3)

set.seed(123)

matrix_glm_mod_w3 <- train(
  form = factor(subday_sz) ~ ., # Convert to factor for binomial classification
  data = matrix_trn,
  trControl = trainControl(method = "cv", number = 10),
  method = "glm",
  family = "binomial",
  weights = weights
)

print("...w3 logistic regression model built...")

## Save model results
precision <- precision(predict(matrix_glm_mod_w3, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(matrix_glm_mod_w3, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(matrix_glm_mod_w3, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(matrix_glm_mod_w3, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
cfm <- confusionMatrix(predict(matrix_glm_mod_w3, matrix_tst), matrix_tst$subday_sz)
accuracy <- cfm$overall["Accuracy"]
kappa <- cfm$overall["Kappa"]
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])


results_tab["log_regress_caret_w3", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[results_row] <- precision
results_tab$recall[results_row] <- recall
results_tab$F1[results_row] <- f1
results_tab$AUC[results_row] <- auc
results_tab$AUCPR[results_row] <- aucpr
results_tab$accuracy[results_row] <- accuracy
results_tab$kappa[results_row] <- kappa

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

results_row <- results_row + 1

#
print("Building w5 logistic regression model...")

weights <- ifelse(matrix_trn$subday_sz == "0",
                  0.5,
                  5)

set.seed(123)

matrix_glm_mod_w5 <- train(
  form = factor(subday_sz) ~ ., # Convert to factor for binomial classification
  data = matrix_trn,
  trControl = trainControl(method = "cv", number = 10),
  method = "glm",
  family = "binomial",
  weights = weights
)

print("...w5 logistic regression model built...")

## Save model results
precision <- precision(predict(matrix_glm_mod_w5, matrix_tst), matrix_tst$subday_sz, relevant = "1")
recall <- recall(predict(matrix_glm_mod_w5, matrix_tst), matrix_tst$subday_sz, relevant = "1")
f1 <- F_meas(predict(matrix_glm_mod_w5, matrix_tst), matrix_tst$subday_sz, relevant = "1", beta = 1, na.rm = TRUE)
pred <- prediction(as.numeric(predict(matrix_glm_mod_w5, matrix_tst)), as.numeric(matrix_tst$subday_sz)) # Standardize data for evaluation
auc <- as.numeric(performance(pred, "auc")@y.values)
aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
cfm <- confusionMatrix(predict(matrix_glm_mod_w5, matrix_tst), matrix_tst$subday_sz)
accuracy <- cfm$overall["Accuracy"]
kappa <- cfm$overall["Kappa"]
fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])


results_tab["log_regress_caret_w5", ] <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA")
results_tab$precision[results_row] <- precision
results_tab$recall[results_row] <- recall
results_tab$F1[results_row] <- f1
results_tab$AUC[results_row] <- auc
results_tab$AUCPR[results_row] <- aucpr
results_tab$accuracy[results_row] <- accuracy
results_tab$kappa[results_row] <- kappa

print("...Weighted logistic regression models complete.")

print("saving results...")

results_tab_save <- results_tab %>% rownames_to_column(var = "model_type")
write.csv(results_tab_save, "results/test_tab.csv")

print("...results saved.")

print("seizure_predict_models.R file run complete")
