
#' Summarise prediction performance metrics
#'
#' Quick function to calculate performance metrics: confusion matrix, accuracy
#' and balanced accuracy for classification; ROC AUC for binary classification;
#' RMSE, R^2 and MAE for regression. Multi-class AUC is returned for multinomial
#' classification.
#'
#' @param output data.frame with columns `testy` containing observed response
#'   from test folds; `predy` predicted response; `predyp` (optional) predicted
#'   probabilities for classification to calculate ROC AUC. For multiclass
#'   `output`, columns 3 onwards contain probabilities for each class in
#'   columns.
#' @param family Optional character value to support specific glmnet models e.g.
#'   'mgaussian', 'cox'.
#' @return An object of class 'predSummary'. For classification a list is
#'   returned containing the confusion matrix table and a vector containing
#'   accuracy and balanced accuracy for classification, ROC AUC for
#'   classification. For regression a vector containing RMSE, R^2 and MAE is
#'   returned. For glmnet 'cox' models, Harrell's C-index is returned.
#'
#'   For glmnet 'mgaussian' models, an object of class 'predSummaryMulti' is
#'   returned which is a list of vectors with regression metrics (RMSE, R^2,
#'   MAE) for each response variable (i.e. each `y` column).
#' @details For multinomial classification, multi-class AUC as defined by Hand
#' and Till is calculated using [pROC::multiclass.roc()].
#'
#' Multi-class balanced accuracy is calculated as the mean of the Recall for
#' each class.
#'
#' R^2 is calculated as 1 - rss / tss, where rss = residual sum of squares, tss
#' = total sum of squares.
#' @seealso [metrics()]
#' @export
predSummary <- function(output, family = "") {
  if (family == "mgaussian") {
    nc <- ncol(output) /2
    summary <- lapply(1:nc, function(i) {
      df <- data.frame(testy = output[, i], predy = output[, i+nc])
      metrics_reg(df)
    })
    names(summary) <- colnames(output)[nc+ 1:nc]
    class(summary) <- "predSummaryMulti"
    return(summary)
  } else if (family == "cox") {
    C_ind <- glmnet::Cindex(output[, 3], output[, 1:2])
    summary <- setNames(C_ind, "C-index")
    class(summary) <- "predSummary"
    return(summary)
  }
  
  if (is.character(output$testy)) {
    output$testy <- factor(output$testy)
  }
  if (is.factor(output$testy)) {
    if (is.character(output$predy)) {
      output$predy <- factor(output$predy, levels=levels(output$testy))
    }
    cm <- table(output$predy, output$testy, dnn=c("Predicted", "Reference"))
    acc <- sum(diag(cm))/ sum(cm)
    b_acc <- multi_bal_acc(cm)
    ccm <- caret::confusionMatrix(cm)
    if (nlevels(output$testy) == 2) {
      outputroc <- pROC::roc(output$testy, output$predyp, direction = "<", 
                             quiet = TRUE)
      auc <- outputroc$auc
      metrics <- setNames(c(auc, acc, b_acc), c("AUC", "Accuracy", "Balanced accuracy"))
    } else {
      auc <- try(pROC::multiclass.roc(output$testy, output[, -c(1,2)])$auc,
                 silent = TRUE)
      if (inherits(auc, "try-error")) auc <- NA
      metrics <- setNames(c(auc, acc, b_acc), c("Multiclass AUC", "Accuracy",
                                                "Balanced accuracy"))
    }
    summary <- list(table = cm, metrics = metrics)
  } else {
    # regression
    summary <- metrics_reg(output)
  }
  class(summary) <- "predSummary"
  summary
}


#' @export
print.predSummary <- function(x, 
                              digits = max(3L, getOption("digits") - 3L),
                              ...) {
  if (is.list(x)) {
    print(x$table)
    cat("\n")
    print(x$metrics, digits = digits, print.gap = 3L)
  } else print(unclass(x), digits = digits, print.gap = 3L)
}


#' @export
print.predSummaryMulti <- function(x, 
                                   digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  print(unclass(x), digits = digits, print.gap = 3L)
}


metrics_reg <- function(output) {
  pred <- output$predy
  obs <- output$testy
  
  rmse <- sqrt(mean((pred - obs)^2))
  mae <- mean(abs(pred - obs))
  
  rss <- sum((pred - obs)^2)
  tss <- sum((obs - mean(obs))^2)
  Rsq <- 1 - rss/tss
  
  setNames(c(rmse, Rsq, mae), c("RMSE", "R.squared", "MAE"))
}

# multiclass balanced accuracy
# table with reference in columns, predicted in rows
multi_bal_acc <- function(tab) {
  recall <- diag(tab) / colSums(tab)
  mean(recall)
}
