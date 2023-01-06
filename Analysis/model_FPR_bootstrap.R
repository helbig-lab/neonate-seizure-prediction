#' False Positive Rate and Confidence Interval generator
#'
#' Loads models and generates their performance metrics
#' Including FPR for desired H2o models and CIs for all models

library(librarian)
librarian::shelf(ggplot2, tidyverse, scales, lubridate, h2o, glue, ROCR,
                 bit64, caret, rpart, rpart.plot, ranger, randomForest,
                 svglite, boot, purrr, rsample)

### Load raw data to recreate test/train sets ----
h2o.init()

message("when switching between all neonates and HIE, go to neonate source file and change 'switch' to on/off respectively")
message("If using an already created dataset, load the desired file as 'matrix_new' and skip to splitting train/test data.")
source("neonate_seizure_predict_matrix.R")
source("../utils.R")

## Generate dataset for all models
matrix_new <- day_1_matrix %>%
  mutate(across(
  .cols = everything(),
  .fns = ~factor(.),
  .names = NULL)) # Converting all numeric columns into binary categorical variables

matrix_new <- matrix_new %>% dplyr::select(-c("hie_th", "events"))
matrix_new <- na.omit(matrix_new, cols = c("sex", "race", "graphoelements", "voltage", "continuity",
                                         "transients", "variability", "reactivity", "impression",
                                         "seizures", "subdaySz"))
matrix_new <- matrix_new %>%
  dplyr::select(-c("END_EXAM_DTTM", "month", "year")) # Only needed for iterative full model analyses

# Split data
set.seed(123)
matrix_idx <- createDataPartition(matrix_new$subdaySz, p = 0.80, list = FALSE)
matrix_trn <- matrix_new[matrix_idx, ]
matrix_tst <- matrix_new[-matrix_idx, ]

### Recreate matrix to get H2o test data only ----
h2o_model_path <- "your_path_here"
temp_model <- h2o.loadModel(paste0(h2o_model_path, "your_model_here"))

# Convert for h2o
y <- "subdaySz"
x <- setdiff(names(matrix_trn), y)

# Convert to h2o object
train_h2o <- as.h2o(matrix_trn)
test_h2o <- as.h2o(matrix_tst)

{
  weight_absent_prompt <- as.numeric(readline(prompt = "Please select desired weight for H2o model when seizures are absent: "))
  weight_present_prompt <- as.numeric(readline(prompt = "Please select desired weight for H2o model when seizures are present: "))
}

if (!is.na(weight_absent_prompt) & !is.na(weight_present_prompt)) {

  message(glue("Using selected weight {weight_absent_prompt} and {weight_present_prompt} respectively for H2o model"))

  warning("If weights do not match loaded test model, results will vary from authors' initial results")

  w0 <- weight_absent_prompt
  w1 <- weight_present_prompt

  } else {

    print("These variables were not recognized, continuing analysis with default values")

    w0 <- 0.6068152
    w1 <- 2.840491

  }

weight_tst <- as.data.frame(test_h2o) %>% mutate(weights =
                                                   case_when(subdaySz == "0" ~ w0,
                                                             subdaySz == "1" ~ w1))
weight_tst <- as.h2o(weight_tst)

### FALSE POSITIVE RATE ----

# In addition to FPR, check some performance metrics,
# Including non-FPR metrics to double check re-created model results against previously created summary table results
# Note this is only works if "Recreate matrix to get H2o test data only" above is run first
temp_perf <- h2o.performance(temp_model, newdata = weight_tst)

h2o.confusionMatrix(temp_perf)
cohen_kappa(temp_perf)
h2o.F1(temp_perf, thresholds = model_threshold(temp_perf))[[1]]
h2o.precision(temp_perf, thresholds = model_threshold(temp_perf))[[1]]
h2o.recall(temp_perf, thresholds = model_threshold(temp_perf))[[1]]
h2o.aucpr(temp_perf)
h2o.fpr(temp_perf, thresholds = model_threshold(temp_perf))[[1]]

### Load non-H20 models for boostrap confidence intervals ----

source_section("seizure_predict_models.R",
               "Build non-H2o Models...",
               "Non-H2o Models built.")

source_section("seizure_predict_models.R",
               "Building weighted logistic regression models...",
               "...Weighted logistic regression models complete.")

### Create reference matrix of all models and their parameters for boostrap loop ----
model_matrix <- setNames(data.frame(matrix(ncol = 6, nrow = 1)),
                         c("model_name", "model", "package_type", "wt", "w0", "w1"))

model_matrix[1, ] <- c("log_regress_caret", "matrix_glm_mod", "caret", FALSE, NA, NA)
model_matrix[2, ] <- c("regresstree_caret", "matrix_tree_mod", "caret", FALSE, NA, NA)
model_matrix[3, ] <- c("random_model1", "model_1", "randomForest", FALSE, NA, NA)
model_matrix[4, ] <- c("random_model_mtry", "model_2", "randomForest", FALSE, NA, NA)
model_matrix[5, ] <- c("random_model_opt", "model_3", "randomForest", FALSE, NA, NA)
model_matrix[6, ] <- c("range_1", "range_1", "ranger", FALSE, NA, NA)
model_matrix[7, ] <- c("log_regress_caret_wb", "matrix_glm_mod_w", "caret", FALSE, NA, NA) # Note: log_regress models don't need explicit weighting for the for loop
model_matrix[8, ] <- c("log_regress_caret_w3", "matrix_glm_mod_w3", "caret", FALSE, NA, NA)
model_matrix[9, ] <- c("log_regress_caret_w5", "matrix_glm_mod_w5", "caret", FALSE, NA, NA)

## Load H2o Models based on all neonates OR HIE only into matrix
if (switch == "off")  {
  h2o_model_path <- "your_path_here"

  h2o_1 <- h2o.loadModel(paste0(h2o_model_path, "rf_h2o_1"))
  model_matrix[10, ] <- c("h2o_1", "h2o_1", "h2o", FALSE, NA, NA)

  h2o_balanced <- h2o.loadModel(paste0(h2o_model_path, "rf_h2o_balanced"))
  model_matrix[11, ] <- c("h2o_balanced", "h2o_balanced", "h2o", FALSE, NA, NA)

  h2o_custom_bal <- h2o.loadModel(paste0(h2o_model_path, "rf_h2o_custom_balanced"))
  model_matrix[12, ] <- c("h2o_custom_bal", "h2o_custom_bal", "h2o", FALSE, NA, NA)

  h2o_weighted_0.6068152_2.840491 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.6068152_2.840491_model_X"))
  model_matrix[13, ] <- c("h2o_weighted_0.6068152_2.840491", "h2o_weighted_0.6068152_2.840491", "h2o", TRUE, 0.6068152, 2.840491)

  h2o_weighted_0.5_1.5 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_1.5_model_X"))
  model_matrix[14, ] <- c("h2o_weighted_0.5_1.5", "h2o_weighted_0.5_1.5", "h2o", TRUE, 0.5, 1.5)

  h2o_weighted_0.5_2 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_2_model_X"))
  model_matrix[15, ] <- c("h2o_weighted_0.5_2", "h2o_weighted_0.5_2", "h2o", TRUE, 0.5, 2.0)

  h2o_weighted_0.5_3 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_3_model_X"))
  model_matrix[16, ] <- c("h2o_weighted_0.5_3", "h2o_weighted_0.5_3", "h2o", TRUE, 0.5, 3.0)

  h2o_weighted_0.5_4 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_4_model_X"))
  model_matrix[17, ] <- c("h2o_weighted_0.5_4", "h2o_weighted_0.5_4", "h2o", TRUE, 0.5, 4.0)

  h2o_weighted_0.5_5 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_5_model_X"))
  model_matrix[18, ] <- c("h2o_weighted_0.5_5", "h2o_weighted_0.5_5", "h2o", TRUE, 0.5, 5.0)

  h2o_weighted_0.5_10 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_10_model_X"))
  model_matrix[19, ] <- c("h2o_weighted_0.5_10", "h2o_weighted_0.5_10", "h2o", TRUE, 0.5, 10.0)

  h2o_weighted_0.5_15 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_15_model_X"))
  model_matrix[20, ] <- c("h2o_weighted_0.5_15", "h2o_weighted_0.5_15", "h2o", TRUE, 0.5, 15.0)

} else {
  h2o_model_path <- "your_path_here"

  h2o_1 <- h2o.loadModel(paste0(h2o_model_path, "rf_h2o_1"))
  model_matrix[10, ] <- c("h2o_1", "h2o_1", "h2o", FALSE, NA, NA)

  h2o_balanced <- h2o.loadModel(paste0(h2o_model_path, "rf_h2o_balanced"))
  model_matrix[11, ] <- c("h2o_balanced", "h2o_balanced", "h2o", FALSE, NA, NA)

  h2o_custom_bal <- h2o.loadModel(paste0(h2o_model_path, "rf_h2o_custom_balanced"))
  model_matrix[12, ] <- c("h2o_custom_bal", "h2o_custom_bal", "h2o", FALSE, NA, NA)

  h2o_weighted_0.6068152_2.840491 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.6068152_2.840491_model_X"))
  model_matrix[13, ] <- c("h2o_weighted_0.6068152_2.840491", "h2o_weighted_0.6068152_2.840491", "h2o", TRUE, 0.6068152, 2.840491)

  h2o_weighted_0.5_1.5 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_1.5_model_X"))
  model_matrix[14, ] <- c("h2o_weighted_0.5_1.5", "h2o_weighted_0.5_1.5", "h2o", TRUE, 0.5, 1.5)

  h2o_weighted_0.5_2 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_2_model_X"))
  model_matrix[15, ] <- c("h2o_weighted_0.5_2", "h2o_weighted_0.5_2", "h2o", TRUE, 0.5, 2.0)

  h2o_weighted_0.5_3 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_3_model_X"))
  model_matrix[16, ] <- c("h2o_weighted_0.5_3", "h2o_weighted_0.5_3", "h2o", TRUE, 0.5, 3.0)

  h2o_weighted_0.5_4 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_4_model_X"))
  model_matrix[17, ] <- c("h2o_weighted_0.5_4", "h2o_weighted_0.5_4", "h2o", TRUE, 0.5, 4.0)

  h2o_weighted_0.5_5 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_5_model_X"))
  model_matrix[18, ] <- c("h2o_weighted_0.5_5", "h2o_weighted_0.5_5", "h2o", TRUE, 0.5, 5.0)

  h2o_weighted_0.5_10 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_10_model_X"))
  model_matrix[19, ] <- c("h2o_weighted_0.5_10", "h2o_weighted_0.5_10", "h2o", TRUE, 0.5, 10.0)

  h2o_weighted_0.5_15 <- h2o.loadModel(paste0(h2o_model_path, "rf_grid_weighted_0.5_15_model_X"))
  model_matrix[20, ] <- c("h2o_weighted_0.5_15", "h2o_weighted_0.5_15", "h2o", TRUE, 0.5, 15.0)

}

### Iterate through models, generating and saving boostrap confidence intervals per desired metric ----
results_matrix <- setNames(data.frame(matrix(ncol = 11, nrow = 1)),
                         c("model_name", "metric_name", "metric_result",
                           "CI_norm_low", "CI_norm_up", "CI_basic_low", "CI_basic_up",
                           "CI_perc_low", "CI_perc_up", "CI_bca_low", "CI_bca_up")
                         )

metrics <- c("accuracy", "precision", "recall", "f1", "auc", "aucpr", "kappa", "fpr") # List of desired metrics

for (i in seq_len(NROW(model_matrix))) {

  print(paste0("Model: ", model_matrix[i, 1]))

  # Generate bootstrap for specific model
  set.seed(123)
  bootstrap_ml <- boot(data = matrix_tst,
                       statistic = model_perf, # From the "utils.R" script
                       R = 5000,
                       model = get(model_matrix[i, 2]),
                       package_type = model_matrix[i, 3],
                       wt = model_matrix[i, 4],
                       w0 = model_matrix[i, 5],
                       w1 = model_matrix[i, 6],
                       parallel = "multicore",
                       ncpus = 4
  )

  # Iterate through bootstrap results for each metric, saving confidence intervals (CIs)
  for (j in seq_len(NROW(metrics))) {

    metric <- metrics[j]
    print(paste0("Metric: ", metric))

    # Save each respective CI and use it to add to a results matrix, where each metric per model gets a row
    # If resampling values are too small and CI doesn't work, simply adds blank rows
    tryCatch({
        ci <- boot.ci(boot.out = bootstrap_ml, index = j)

        results_matrix <- rbind(results_matrix,
          c(model_matrix[i, 1],
           metric,
           ci$t0,
           ci$normal[, 2],
           ci$normal[, 3],
           ci$basic[, 4],
           ci$basic[, 5],
           ci$percent[, 4],
           ci$percent[, 5],
           ci$bca[, 4],
           ci$bca[, 5]
           )
        )

        # Save bootstrap results iteratively (in case loop breaks mid run)
        if (switch == "off") {
          write.csv("results/", results_matrix, glue("model_cis_all_neonates_{Sys.Date()}csv"))
        } else {
          write.csv("results/", results_matrix, glue("model_cis_hie_only_{Sys.Date()}csv"))
        }

    }, error = function(err) {

      message(glue("{metric} bootstrap does not have enough resampling data to produce results, please try again."))
      message(glue("Error: {err}"))

    }

    )

  }
}

# Save final bootstrap results
if (switch == "off") {
  write.csv("results/", results_matrix, glue("model_cis_all_neonates_{Sys.Date()}.csv"))
} else {
  write.csv("results/", results_matrix, glue("model_cis_hie_only_{Sys.Date()}.csv"))
}
