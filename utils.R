#' Utilities script
#'
#' Contains all necessary functions required for model generation, analysis, etc.

library(librarian)
librarian::shelf(h2o, caret)

# Custom way to calculate average accuracy of model
calc_acc <- function(actual, predicted) {
  mean(actual == predicted)
}

# Calculate kappa for h2o
cohen_kappa <- function(model) {

  po <- sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 2]]) / sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 1]] + h2o.confusionMatrix(model)[[1, 2]] + h2o.confusionMatrix(model)[[2, 2]])

  no <- sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[1, 2]]) / sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 1]] + h2o.confusionMatrix(model)[[1, 2]] + h2o.confusionMatrix(model)[[2, 2]]) * sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 1]]) / sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 1]] + h2o.confusionMatrix(model)[[1, 2]] + h2o.confusionMatrix(model)[[2, 2]])
  yes <- sum(h2o.confusionMatrix(model)[[2, 1]] + h2o.confusionMatrix(model)[[2, 2]]) / sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 1]] + h2o.confusionMatrix(model)[[1, 2]] + h2o.confusionMatrix(model)[[2, 2]]) * sum(h2o.confusionMatrix(model)[[1, 2]] + h2o.confusionMatrix(model)[[2, 2]]) / sum(h2o.confusionMatrix(model)[[1, 1]] + h2o.confusionMatrix(model)[[2, 1]] + h2o.confusionMatrix(model)[[1, 2]] + h2o.confusionMatrix(model)[[2, 2]])

  pe <- yes + no

  k <- (po - pe) / (1 - pe)
  return(k)
}

# Find threshold function for "best" F1 score (ignoring non-fininte values like NaN and NA)
model_threshold <- function(model) {
  model@metrics$thresholds_and_metric_scores$threshold[model@metrics$thresholds_and_metric_scores$f1 == max(model@metrics$thresholds_and_metric_scores$f1, na.rm = TRUE) & is.finite(model@metrics$thresholds_and_metric_scores$f1) == TRUE]
}

# Function for pulling in specific lines only from a source
source_section <- function(file, start_line, finish_line) {
  file_lines <- scan(file, what = character(), sep = "\n")
  start <- grep(start_line, file_lines)
  finish <- grep(finish_line, file_lines)
  source(textConnection(file_lines[(start + 1):(finish - 1)]))
}

# Generate machine learning model performance from a reference table
# Note: "data" here would be the test dataset and any added weight information ONLY applies to h2o models (log_regress can do this implicitly)
model_perf <- function(model, data, idx, package_type, ...) {

  dots <- list(...)

  if (package_type == "h2o") {

    if (dots$wt) {

      # Convert for h2o
      test_h2o <- as.h2o(data[idx, ])

      weight_tst <- as.data.frame(test_h2o) %>% mutate(weights =
                                                         case_when(subdaySz == "0" ~ as.numeric(dots$w0),
                                                                   subdaySz == "1" ~ as.numeric(dots$w1)))
      test_h2o <- as.h2o(weight_tst)
    } else {

      # Convert for h2o
      test_h2o <- as.h2o(data[idx, ])
    }

    data_index <- h2o.performance(model, newdata = test_h2o)

    precision <- h2o.precision(data_index, thresholds = model_threshold(data_index))[[1]]
    recall <- h2o.recall(data_index, thresholds = model_threshold(data_index))[[1]]
    f1 <- h2o.F1(data_index, thresholds = model_threshold(data_index))[[1]]
    auc <- h2o.auc(data_index)
    accuracy <- h2o.accuracy(data_index, thresholds = model_threshold(data_index))[[1]]
    aucpr <- h2o.aucpr(data_index)
    kappa <- cohen_kappa(data_index)
    fpr <- h2o.fpr(data_index, thresholds = model_threshold(data_index))[[1]]
  } else if (package_type == "ranger") {

    data_index <- data[idx, ]

    precision <- precision(predict(model, data_index)$predictions, data_index$subdaySz, relevant = "1")
    recall <- recall(predict(model, data_index)$predictions, data_index$subdaySz, relevant = "1")
    f1 <- F_meas(predict(model, data_index)$predictions, data_index$subdaySz, relevant = "1", beta = 1, na.rm = TRUE)
    pred <- prediction(as.numeric(predict(model, data_index)$predictions), as.numeric(data_index$subdaySz)) # Standardize data for evaluation
    auc <- as.numeric(performance(pred, "auc")@y.values)
    aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
    accuracy <- confusionMatrix(predict(model, data_index)$predictions, data_index$subdaySz)$overall["Accuracy"]
    kappa <- confusionMatrix(predict(model, data_index)$predictions, data_index$subdaySz)$overall["Kappa"]
    cfm <- confusionMatrix(predict(model, data_index)$predictions, data_index$subdaySz)
    fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])
  } else {

    data_index <- data[idx, ]

    precision <- precision(predict(model, data_index), data_index$subdaySz, relevant = "1")
    recall <- recall(predict(model, data_index), data$subdaySz, relevant = "1")
    f1 <- F_meas(predict(model, data_index), data_index$subdaySz, relevant = "1", beta = 1, na.rm = TRUE)
    pred <- prediction(as.numeric(predict(model, data_index)), as.numeric(data_index$subdaySz)) # Standardize data for evaluation
    auc <- as.numeric(performance(pred, "auc")@y.values)
    cfm <- confusionMatrix(predict(model, data_index), data_index$subdaySz)
    accuracy <-  cfm$overall["Accuracy"]
    aucpr <- as.numeric(performance(pred, "aucpr")@y.values)
    kappa <- cfm$overall[[2]]
    fpr <- cfm$table[2, 1] / (cfm$table[1, 1] + cfm$table[2, 1])
  }
  return(
    c(accuracy, precision, recall, f1, auc, aucpr, kappa, fpr)
  )
}
