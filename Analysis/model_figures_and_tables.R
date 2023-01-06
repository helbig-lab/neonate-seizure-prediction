#' Seizure prediction model figures
#'
#' Visualize seizure prediction model performance
#' and review basic performance statistics


library(librarian)
librarian::shelf(ggplot2, tidyverse, scales, lubridate,
                 h2o, ROCR, bit64, caret, svglite)

source("../utils.R")

### Load and extract models ----
h2o.init()

model_path <- paste0("your_path_here") # Single Model

### Review single model results ----

# Load raw data to recreate test/train sets and load h2o model
temp_model <- h2o.loadModel(model_path)

message("If using an already created dataset, load the desired file as 'matrix_new' and skip to splitting train/test data.")
source("neonate_seizure_predict_matrix.R")

## Recreate matrix to get H2o test data
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
  dplyr::select(-c("END_EXAM_DTTM", "month", "year")) # Only needed for iterative model analyses

set.seed(123)
matrix_idx <- createDataPartition(matrix_new$subdaySz, p = 0.80, list = FALSE)
matrix_trn <- matrix_new[matrix_idx, ]
matrix_tst <- matrix_new[-matrix_idx, ]

# Convert for h2o
y <- "subdaySz"
x <- setdiff(names(matrix_trn), y)

train_h2o <- as.h2o(matrix_trn)
test_h2o <- as.h2o(matrix_tst)

{
  weight_absent_prompt <- as.numeric(readline(prompt = "Please select desired weight for H2o model when seizures are absent: "))
  weight_present_prompt <- as.numeric(readline(prompt = "Please select desired weight for H2o model when seizures are present: "))
}

if (!is.na(weight_absent_prompt) & !is.na(weight_present_prompt)) {
  
  message(glue("Using selected weight {weight_absent_prompt} and {weight_present_prompt} respectively for H2o model"))
  
  warning("If weights do not match loaded test model, results will vary from initial results")
  
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

## Check performance and save results
temp_perf <- h2o.performance(temp_model, newdata = weight_tst)

h2o.confusionMatrix(temp_perf)
cohen_kappa(temp_perf)
h2o.F1(temp_perf, thresholds = model_threshold(temp_perf))[[1]]
h2o.precision(temp_perf, thresholds = model_threshold(temp_perf))[[1]]
h2o.recall(temp_perf, thresholds = model_threshold(temp_perf))[[1]]

model <- as.data.frame(h2o.metric(temp_perf))
write.csv(model, "figures_and_tables/pr_data.csv")

## Plot results (change x and y based on desired metrics for review)
ggplot(model, aes(x = recall, y = precision)) +
  geom_ribbon(aes(x = recall, ymax = precision), ymin = 0, alpha = 0.2) +
  geom_line(col = "#238A8DFF", size = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  scale_fill_gradient2(low = "#DCE318FF", high = "#440154FF",
                        midpoint = median(model$recall)) +
  xlab("Recall") +
  ylab("Precision") +
  theme_classic()
ggsave("figures_and_tables/AUCPR_curve.svg", dpi = 300, width = 8, height = 8)

### Variable Importance ----

varimp <- h2o.varimp(temp_model)

ggplot(varimp, aes(y = reorder(variable, percentage), x = percentage)) +
  geom_bar(stat = "identity", fill = "#440154FF") +
  ylab("Variable") +
  xlab("Importance (%)") +
  theme_minimal() +
  scale_x_continuous(labels = scales::percent)
ggsave("figures_and_tables/importance.svg", dpi = 300, width = 12, height = 8)
