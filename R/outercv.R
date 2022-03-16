# Outer CV with rf


#' Outer cross-validation with random forest
#' 
#' Outer cross-validation (CV) with random forest. This is designed to evaluate
#' performance of the model with fixed hyperparameters and no tuning. If tuning
#' of parameters on data is required, this will need full nested CV with inner
#' CV to tune model hyperparameters (see [nestcv.train]).
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter]. 
#' Any function can be provided and is passed `y` and `x`. Must return a 
#' character vector with names of filtered predictors.
#' @param filter_options List of additional arguments passed to the filter 
#' function specified by `filterFUN`. 
#' @param n_outer_folds Number of outer CV folds
#' @param cores Number of cores for parallel processing. Note this currently 
#' uses [parallel::mclapply].
#' @param ... Optional arguments passed to [randomForest].
#' @return An object with S3 class "outercv.rf"
#' @importFrom caret createFolds confusionMatrix defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#' @importFrom pROC roc
#' @importFrom randomForest randomForest
#' @importFrom stats predict setNames
#' @export
#' 
outercv.rf <- function(y, x,
                    filterFUN = NULL,
                    filter_options = NULL,
                    n_outer_folds = 10,
                    cores = 1,
                    ...) {
  reg <- !(is.factor(y) | is.character(y))  # y = regression
  outer_folds <- createFolds(y, k = n_outer_folds)
  outer_res <- mclapply(1:n_outer_folds, function(i) {
    test <- outer_folds[[i]]
    filtx <- if (is.null(filterFUN)) x else {
      args <- list(y = y[-test], x = x[-test, ])
      args <- append(args, filter_options)
      fset <- do.call(filterFUN, args)
      x[, fset]
    }
    fit <- randomForest(x = filtx[-test, ], y = y[-test], ...)
    # test on outer CV
    predy <- as.vector(predict(fit, newdata = filtx[test, ], type = "response"))
    preds <- data.frame(predy=predy, testy=y[test])
    # for AUC
    if (!reg & nlevels(y) == 2) {
      predyp <- predict(fit, newdata = filtx[test, ], type = "prob")
      predyp <- predyp[,2]
      preds$predyp <- predyp
    }
    rownames(preds) <- rownames(x[test, ])
    list(preds = preds,
         fit = fit,
         nfilter = ncol(filtx))
  }, mc.cores = cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  fit.roc <- NULL
  if (reg) {
    df <- data.frame(obs = output$testy, pred = output$predy)
    summary <- caret::defaultSummary(df)
  } else if (nlevels(y) == 2) {
    cm <- table(output$predy, output$testy)
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    fit.roc <- pROC::roc(output$testy, output$predyp, direction = "<")
    auc <- fit.roc$auc
    summary <- setNames(c(auc, acc, b_acc), c("AUC", "Accuracy", "Balanced accuracy"))
  } else {
    cm <- table(output$predy, output$testy)
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    summary <- setNames(c(acc, b_acc), c("Accuracy", "Balanced accuracy"))
  }
  
  # fit final model
  filtx <- if (is.null(filterFUN)) x else {
    args <- list(y = y, x = x)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    x[, fset]
  }
  fit <- randomForest(filtx, y, ...)
  out <- list(output = output,
              outer_result = outer_res,
              outer_folds = outer_folds,
              final_fit = fit,
              final_vars = colnames(filtx),
              roc = fit.roc,
              summary = summary)
  class(out) <- "outercv.rf"
  out
}


