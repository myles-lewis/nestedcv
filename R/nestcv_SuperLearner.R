
#' Outer cross-validation of SuperLearner model
#'
#' Provides a single loop of outer cross-validation to evaluate performance of
#' ensemble models from `SuperLearner` package.
#'
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter].
#'   Any function can be provided and is passed `y` and `x`. Must return a
#'   character vector with names of filtered predictors. Not available if
#'   `outercv` is called with a formula.
#' @param filter_options List of additional arguments passed to the filter
#'   function specified by `filterFUN`.
#' @param weights Weights applied to each sample for models which can use
#'   weights. Note `weights` and `balance` cannot be used at the same time.
#'   Weights are not applied in filters.
#' @param balance Specifies method for dealing with imbalanced class data.
#'   Current options are `"randomsample"` or `"smote"`. Not available if
#'   `outercv` is called with a formula. See [randomsample()] and [smote()]
#' @param balance_options List of additional arguments passed to the balancing
#'   function
#' @param outer_method String of either `"cv"` or `"LOOCV"` specifying whether
#'   to do k-fold CV or leave one out CV (LOOCV) for the outer folds
#' @param n_outer_folds Number of outer CV folds
#' @param outer_folds Optional list containing indices of test folds for outer
#'   CV. If supplied, `n_outer_folds` is ignored.
#' @param cv.cores Number of cores for parallel processing of the outer loops.
#'   NOTE: this uses `parallel::mclapply` on unix/mac and `parallel::parLapply`
#'   on windows.
#' @param na.option Character value specifying how `NA`s are dealt with.
#'   `"omit"` is equivalent to `na.action = na.omit`. `"omitcol"` removes cases
#'   if there are `NA` in 'y', but columns (predictors) containing `NA` are
#'   removed from 'x' to preserve cases. Any other value means that `NA` are
#'   ignored (a message is given).
#' @param ... Optional arguments passed to [SuperLearner::SuperLearner()]
#' @return An object with S3 class "nestcv.SuperLearner"
#'   \item{call}{the matched call}
#'   \item{output}{Predictions on the left-out outer folds}
#'   \item{outer_result}{List object of results from each outer fold containing
#'   predictions on left-out outer folds, model result and number of filtered
#'   predictors at each fold.}
#'   \item{dimx}{vector of number of observations and number of predictors}
#'   \item{y}{original response vector}
#'   \item{yfinal}{final response vector (post-balancing)}
#'   \item{outer_folds}{List of indices of outer test folds}
#'   \item{final_fit}{Final fitted model on whole data}
#'   \item{final_vars}{Column names of filtered predictors entering final model}
#'   \item{summary_vars}{Summary statistics of filtered predictors}
#'   \item{roc}{ROC AUC for binary classification where available.}
#'   \item{summary}{Overall performance summary. Accuracy and balanced accuracy
#'   for classification. ROC AUC for binary classification. RMSE for
#'   regression.}
#'   
#' @importFrom SuperLearner SuperLearner
#' @export

nestcv.SuperLearner <- function(y, x,
                                filterFUN = NULL,
                                filter_options = NULL,
                                weights = NULL,
                                balance = NULL,
                                balance_options = NULL,
                                outer_method = c("cv", "LOOCV"),
                                n_outer_folds = 10,
                                outer_folds = NULL,
                                cv.cores = 1,
                                na.option = "pass",
                                ...) {
  ncv.call <- match.call(expand.dots = TRUE)
  ok <- checkxy(y, x, na.option, weights)
  y <- y[ok$r]
  x <- x[ok$r, ok$c]
  weights <- weights[ok$r]
  if (!is.null(balance) & !is.null(weights)) {
    stop("`balance` and `weights` cannot be used at the same time")}
  if (is.character(y)) y <- factor(y)
  reg <- !is.factor(y)  # y = regression
  if (!is.null(balance) & reg) {
    stop("`balance` can only be used for classification")}
  
  outer_method <- match.arg(outer_method)
  if (is.null(outer_folds)) {
    outer_folds <- switch(outer_method,
                          cv = createFolds(y, k = n_outer_folds),
                          LOOCV = 1:length(y))
  }
  
  if (Sys.info()["sysname"] == "Windows" & cv.cores >= 2) {
    cl <- makeCluster(cv.cores)
    dots <- list(...)
    clusterExport(cl, varlist = c("outer_folds", "y", "x", "superLearner",
                                  "filterFUN", "filter_options",
                                  "weights", "balance", "balance_options",
                                  "nestSLcore", "dots"),
                  envir = environment())
    on.exit(stopCluster(cl))
    outer_res <- parLapply(cl = cl, outer_folds, function(test) {
      args <- c(list(test=test, y=y, x=x,
                     filterFUN=filterFUN, filter_options=filter_options,
                     weights=weights, balance=balance,
                     balance_options=balance_options), dots)
      do.call(nestSLcore, args)
    })
  } else {
    outer_res <- mclapply(outer_folds, function(test) {
      nestSLcore(test, y, x,
                 filterFUN, filter_options, weights,
                 balance, balance_options, ...)
    }, mc.cores = cv.cores)
  }
  
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  if (!is.null(rownames(x))) {
    rownames(output) <- unlist(lapply(predslist, rownames))}
  
  summary <- predSummary(output)
  if (!reg & nlevels(y) == 2) {
    fit.roc <- pROC::roc(output$testy, output$predyp, direction = "<",
                         quiet = TRUE)
  } else fit.roc <- NULL
  
  # fit final model
  dat <- nest_filt_bal(NULL, y, x, filterFUN, filter_options,
                       balance, balance_options)
  yfinal <- dat$ytrain
  filtx <- dat$filt_xtrain
  Y <- if (reg) yfinal else as.numeric(yfinal) -1
  
  fit <- SuperLearner(Y = Y, X = data.frame(filtx), obsWeights = weights, ...)
  
  out <- list(call = ncv.call,
              output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              outer_folds = outer_folds,
              dimx = dim(x),
              y = y,
              yfinal = yfinal,
              final_fit = fit,
              final_vars = colnames(filtx),
              summary_vars = summary_vars(filtx),
              roc = fit.roc,
              summary = summary)
  class(out) <- "nestcv.SuperLearner"
  out
}


nestSLcore <- function(test, y, x,
                       filterFUN, filter_options, weights,
                       balance, balance_options, ...) {
  dat <- nest_filt_bal(test, y, x, filterFUN, filter_options,
                       balance, balance_options)
  ytrain <- dat$ytrain
  ytest <- dat$ytest
  filt_xtrain <- data.frame(dat$filt_xtrain)
  filt_xtest <- data.frame(dat$filt_xtest)
  
  reg <- !is.factor(y)
  Y <- if (reg) ytrain else as.numeric(ytrain) -1
  fit <- SuperLearner(Y = Y, X = filt_xtrain,
                      obsWeights = weights[-test], ...)
  
  # test on outer CV
  predSL <- predict(fit, newdata = filt_xtest, onlySL = TRUE)
  if (reg) {
    predy <- c(predSL$pred)
  } else {
    # convert prob to class
    predy <- levels(y)[as.numeric(predSL$pred > 0.5) +1]
  }
  
  preds <- data.frame(predy=predy, testy=ytest)
  # for AUC
  if (!reg & nlevels(y) == 2) {
    preds$predyp <- c(predSL$pred)
  }
  rownames(preds) <- rownames(filt_xtest)
  list(preds = preds,
       fit = fit,
       nfilter = ncol(filt_xtest),
       ytrain = ytrain)
}


#' @importFrom matrixStats rowSds
#' @export
summary.nestcv.SuperLearner <- function(object, 
                                   digits = max(3L, getOption("digits") - 3L), 
                                   ...) {
  cat("Nested cross-validation with SuperLearner\n")
  cat("Outer loop: ", switch(object$outer_method,
                             cv = paste0(length(object$outer_folds), "-fold CV"),
                             LOOCV = "leave-one-out CV"))
  cat(paste0("\nInner loop:  ", object$final_fit$cvControl$V,
             "-fold CV (SuperLearner)\n"))
  balance <- object$call$balance
  if (!is.null(balance)) {
    cat("Balancing: ", balance, "\n")
  }
  cat(object$dimx[1], "observations,", object$dimx[2], "predictors\n")
  if (!is.numeric(object$y)) print(c(table(object$y)))
  cat("\nModel: SuperLearner\n")
  
  SLcoef <- lapply(object$outer_result, function(x) x$fit$coef)
  SLcoef <- do.call(cbind, SLcoef)
  SLrisk <- lapply(object$outer_result, function(x) x$fit$cvRisk)
  SLrisk <- do.call(cbind, SLrisk)
  
  if (!is.null(object$call$filterFUN)) {
    cat("Filter: ", object$call$filterFUN, "\n")
    nfilter <- unlist(lapply(object$outer_result, '[[', 'nfilter'))
    nfilter <- data.frame(n.filter = nfilter,
                          row.names = paste("Fold", seq_along(nfilter)))
    print(nfilter, digits = digits, print.gap = 2L)
  } else {
    nfilter <- NULL
    cat("No filter\n")
  }
  
  res <- data.frame(Risk = rowMeans(SLrisk),
                    `Risk SE` = rowSds(SLrisk)/sqrt(ncol(SLrisk)),
                    Coef = rowMeans(SLcoef),
                    `Coef SE` = rowSds(SLcoef)/sqrt(ncol(SLcoef)),
                    check.names = FALSE)
  print(res, digits = digits, print.gap = 2L)
  
  cat("\nFinal fit:")
  print(object$final_fit)
  cat("\nResult:\n")
  print(object$summary, digits = digits, print.gap = 2L)
  out <- list(dimx = object$dimx, nfilter = nfilter,
              result = object$summary)
  invisible(out)
}

