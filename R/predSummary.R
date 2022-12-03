
#' Summarise prediction performance metrics
#' 
#' Quick function to calculate performance metrics: confusion matrix, accuracy
#' and balanced accuracy for classification; ROC AUC for binary classification;
#' RMSE and R^2 for regression.
#' 
#' @param output data.frame with columns `testy` containing observed response
#'   from test folds; `predy` predicted response; `predyp` (optional) predicted
#'   probabilities for classification to calculate ROC AUC
#' @return An object of class 'predSummary'. For classification a list is
#'   returned containing the confusion matrix table and a vector containing
#'   accuracy and balanced accuracy for classification, ROC AUC for binary
#'   classification. For regression a vector containing RMSE and R^2 is
#'   returned.
#' 
#' @export
predSummary <- function(output) {
  if (is.character(output$testy)) {
    output$testy <- factor(output$testy)
  }
  if (is.factor(output$testy)) {
    if (is.character(output$predy)) {
      output$predy <- factor(output$predy, levels=levels(output$testy))
    }
    cm <- table(output$predy, output$testy, dnn=c("Predicted", "Reference"))
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    if (nlevels(output$testy) == 2) {
      outputroc <- pROC::roc(output$testy, output$predyp, direction = "<", 
                             quiet = TRUE)
      auc <- outputroc$auc
      metrics <- setNames(c(auc, acc, b_acc), c("AUC", "Accuracy", "Balanced accuracy"))
    } else {
      metrics <- setNames(c(acc, b_acc), c("Accuracy", "Balanced accuracy"))
    }
    summary <- list(table = cm, metrics = metrics)
  } else {
    df <- data.frame(obs = output$testy, pred = output$predy)
    summary <- caret::defaultSummary(df)
  }
  class(summary) <- "predSummary"
  summary
}


#' @export
print.predSummary <- function(x, ...) {
  if (is.list(x)) {
    print(x$table)
    cat("\n")
    print(x$metrics, ...)
  } else print(unclass(x), ...)
}
