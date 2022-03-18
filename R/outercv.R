# Outer CV with simple models


#' Outer cross-validation of selected models
#'
#' This is a convenience function designed to quickly evaluate performance of
#' specific models (random forest, naive Bayes) with fixed hyperparameters and
#' no tuning. If tuning of parameters on data is required, full nested CV with
#' inner CV is recommended to tune model hyperparameters (see [nestcv.train]).
#'
#' @param model String specifying function to fit. Currently supports
#'   randomForest and naive_bayes.
#' @param y Response vector
#' @param x Matrix of predictors
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter].
#'   Any function can be provided and is passed `y` and `x`. Must return a
#'   character vector with names of filtered predictors.
#' @param filter_options List of additional arguments passed to the filter
#'   function specified by `filterFUN`.
#' @param outer_method String of either `"cv"` or `"LOOCV"` specifying whether
#'   to do k-fold CV or leave one out CV (LOOCV) for the outer folds
#' @param n_outer_folds Number of outer CV folds
#' @param cores Number of cores for parallel processing. Note this currently
#'   uses [parallel::mclapply].
#' @param ... Optional arguments passed to the function specified by `model`.
#' @return An object with S3 class "outercv"
#' @details An alternative method of tuning a single model with fixed parameters
#'   is to use [nestcv.train] with `tuneGrid` set as a single row of a
#'   data.frame. The parameters which are needed for a specific model can be
#'   identified using [caret::modelLookup()].
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
                    outer_method = c("cv", "LOOCV"),
                    n_outer_folds = 10,
                    cores = 1,
                    ...) {
  reg <- !(is.factor(y) | is.character(y))  # y = regression
  outer_method <- match.arg(outer_method)
  outer_folds <- switch(outer_method,
                        cv = createFolds(y, k = n_outer_folds),
                        LOOCV = 1:length(y))
  outer_res <- mclapply(1:length(outer_folds), function(i) {
    test <- outer_folds[[i]]
    filtx <- if (is.null(filterFUN)) x else {
      args <- list(y = y[-test], x = x[-test, ])
      args <- append(args, filter_options)
      fset <- do.call(filterFUN, args)
      x[, fset]
    }
    # check if model uses formula
    if ("formula" %in% formalArgs(model)) {
      dat <- if (is.data.frame(filtx)) {filtx[-test, ]
      } else as.data.frame(filtx[-test, ], stringsAsFactors = TRUE)
      dat$.outcome <- y[-test]
      fit <- model(as.formula(".outcome ~ ."), data = dat, ...)
    } else {
      fit <- model(y = y[-test], x = filtx[-test, ], ...)
    }
    # test on outer CV
    predy <- predict(fit, newdata = filtx[test, ])
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
  # model.args <- list(x = filtx, y = y)
  # fit <- do.call(model, append(model.args, dot.args))
  if ("formula" %in% formalArgs(model)) {
    dat <- if (is.data.frame(filtx)) {filtx
    } else as.data.frame(filtx, stringsAsFactors = TRUE)
    dat$.outcome <- y
    fit <- model(as.formula(".outcome ~ ."), data = dat, ...)
  } else {
    fit <- model(y = y, x = filtx, ...)
  }
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


