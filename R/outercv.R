# outer CV


#' Outer cross-validation
#' 
#' Outer cross-validation (CV) with a specified model. This is designed to
#' evaluate performance of models for which no tuning of hyperparameters is
#' required or for which fixed hyperparameters are provided. If tuning of
#' parameters on data is required, this will need full nested CV with inner CV
#' to tune model hyperparameters (see [nestcv.train]).
#' 
#' @param model Model function to be fitted
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
#' @return An object with S3 class "outercv"
#' @importFrom caret createFolds confusionMatrix defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#' @importFrom pROC roc
#' @importFrom stats predict setNames
#' @export
#' 
outercv <- function(model, y, x,
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
    fit <- model(x = filtx[-test, ], y = y[-test], ...)
    # test on outer CV
    predy <- as.vector(predict(fit, newdata = filtx[test, ], type = "response"))
    # for AUC
    if (!reg & nlevels(y) == 2) {
      predyp <- predict(fit, newdata = filtx[test, ], type = "prob")
      nc <- ncol(predyp)
      if (!is.null(nc)) {
        if (nc == 2) predyp <- predyp[,2]
      }
      preds <- data.frame(predy=predy, predyp=predyp, testy=y[test])
    } else {
      preds <- data.frame(predy=predy, testy=y[test])
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
  fit <- model(filtx, y, ...)
  out <- list(output = output,
              outer_result = outer_res,
              outer_folds = outer_folds,
              final_fit = fit,
              final_vars = colnames(filtx),
              roc = fit.roc,
              summary = summary)
  class(out) <- "outercv"
  out
}



