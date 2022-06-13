# Outer CV with simple models


#' Outer cross-validation of selected models
#'
#' This is a convenience function designed to use a single loop of
#' cross-validation to quickly evaluate performance of specific models (random
#' forest, naive Bayes, lm, glm) with fixed hyperparameters and no tuning. If
#' tuning of parameters on data is required, full nested CV with inner CV is
#' needed to tune model hyperparameters (see [nestcv.train]).
#'
#' @param y Response vector
#' @param x Matrix or data frame of predictors
#' @param formula A formula describing the model to be fitted
#' @param data A matrix or data frame containing variables in the model.
#' @param model Model function to be fitted.
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter].
#'   Any function can be provided and is passed `y` and `x`. Must return a
#'   character vector with names of filtered predictors. Not available if
#'   `outercv` is called with a formula.
#' @param filter_options List of additional arguments passed to the filter
#'   function specified by `filterFUN`.
#' @param outer_method String of either `"cv"` or `"LOOCV"` specifying whether
#'   to do k-fold CV or leave one out CV (LOOCV) for the outer folds
#' @param n_outer_folds Number of outer CV folds
#' @param cv.cores Number of cores for parallel processing. Note this currently
#'   uses [parallel::mclapply].
#' @param predict_type Only used with binary classification. Calculation of ROC
#'   AUC requires predicted class probabilities from fitted models. Most model
#'   functions use syntax of the form `predict(..., type = "prob")`. However,
#'   some models require a different `type` to be specified, which can be passed
#'   to `predict()` via `predict_type`.
#' @param na.option Character value specifying how `NA`s are dealt with.
#'   `"omit"` is equivalent to `na.action = na.omit`. `"omitcol"` removes cases
#'   if there are `NA` in 'y', but columns (predictors) containing `NA` are
#'   removed from 'x' to preserve cases. Any other value means that `NA` are
#'   ignored (a message is given).
#' @param na.action Formula S3 method only: a function to specify the action to
#'   be taken if NAs are found. The default action is for the procedure to fail.
#'   An alternative is `na.omit`, which leads to rejection of cases with missing
#'   values on any required variable. (NOTE: If given, this argument must be
#'   named.)
#' @param ... Optional arguments passed to the function specified by `model`.
#' @return An object with S3 class "outercv"
#'   \item{call}{the matched call}
#'   \item{output}{Predictions on the left-out outer folds}
#'   \item{outer_result}{List object of results from each outer fold containing
#'   predictions on left-out outer folds, model result and number of filtered
#'   predictors at each fold.}
#'   \item{dimx}{vector of number of observations and number of predictors}
#'   \item{outer_folds}{List of indices of outer training folds}
#'   \item{final_fit}{Final fitted model on whole data}
#'   \item{final_vars}{Column names of filtered predictors entering final model}
#'   \item{roc}{ROC AUC for binary classification where available.}
#'   \item{summary}{Overall performance summary. Accuracy and balanced accuracy
#'   for classification. ROC AUC for binary classification. RMSE for
#'   regression.}
#' @details 
#'   Some predictive model functions do not have an x & y interface. If the
#'   function specified by `model` requires a formula, `x` & `y` will be merged
#'   into a data.frame with `model()` called with a formula equivalent to 
#'   `y ~ .`.
#'   
#'   An alternative method of tuning a single model with fixed parameters
#'   is to use [nestcv.train] with `tuneGrid` set as a single row of a
#'   data.frame. The parameters which are needed for a specific model can be
#'   identified using [caret::modelLookup()].
#'
#'   The formula method for `outercv` is not recommended with large data sets -
#'   it is designed primarily for more basic models e.g. `lm()` and `glm()`.
#'   Filtering is not available if `outercv` is called with a formula - use `x`
#'   and `y` interface instead.
#'
#'   Note that in the case of `model = lm`, although additional arguments e.g.
#'   `subset`, `weights`, `offset` are passed into the model function via
#'   `"..."` the scoping is known to go wrong. Avoid using these arguments with
#'   `model = lm`.
#'   
#'   `NA` handling differs between the default S3 method and the formula S3
#'   method. The `na.option` argument takes a character string, while the more
#'   typical `na.action` argument takes a function.
#'   
#' @importFrom caret createFolds confusionMatrix defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom methods formalArgs
#' @importFrom parallel mclapply
#' @importFrom pROC roc
#' @importFrom stats as.formula predict setNames
#' @examples
#' 
#' ## Classification example
#' 
#' ## sigmoid function
#' sigmoid <- function(x) {1 / (1 + exp(-x))}
#' 
#' # load iris dataset and simulate a binary outcome
#' data(iris)
#' dt <- iris[, 1:4]
#' colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
#' dt <- as.data.frame(apply(dt, 2, scale))
#' x <- dt
#' y2 <- sigmoid(0.5 * dt$marker1 + 2 * dt$marker2) > runif(nrow(dt))
#' y2 <- factor(y2)
#' 
#' ## Random forest
#' library(randomForest)
#' cvfit <- outercv(y2, x, randomForest)
#' summary(cvfit)
#' plot(cvfit$roc)
#' 
#' ## Mixture discriminant analysis (MDA)
#' if (requireNamespace("mda", quietly = TRUE)) {
#'   library(mda)
#'   cvfit <- outercv(y2, x, mda, predict_type = "posterior")
#'   summary(cvfit)
#' }
#' 
#' 
#' ## Example with continuous outcome
#' y <- -3 + 0.5 * dt$marker1 + 2 * dt$marker2 + rnorm(nrow(dt), 0, 2)
#' dt$outcome <- y
#' 
#' ## simple linear model - formula interface
#' cvfit <- outercv(outcome ~ ., data = dt, model = lm)
#' summary(cvfit)
#' 
#' ## random forest for regression
#' cvfit <- outercv(y, x, randomForest)
#' summary(cvfit)
#' 
#' ## example with lm_filter() to reduce input predictors
#' cvfit <- outercv(y, x, randomForest, filterFUN = lm_filter,
#'                  filter_options = list(nfilter = 2))
#' summary(cvfit)
#' 
#' @export
#' 
outercv <- function(y, ...) {
  UseMethod("outercv")
}


#' @rdname outercv
#' @export
#' 
outercv.default <- function(y, x,
                            model,
                            filterFUN = NULL,
                            filter_options = NULL,
                            outer_method = c("cv", "LOOCV"),
                            n_outer_folds = 10,
                            cv.cores = 1,
                            predict_type = "prob",
                            na.option = "pass",
                            ...) {
  outercv.call <- match.call(expand.dots = TRUE)
  ok <- checkxy(y, x, na.option)
  y <- y[ok$r]
  x <- x[ok$r, ok$c]
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
      predyp <- predict(fit, newdata = filtx[test, ], type = predict_type)
      predyp <- predyp[,2]
      preds$predyp <- predyp
    }
    rownames(preds) <- rownames(x[test, , drop = FALSE])
    list(preds = preds,
         fit = fit,
         nfilter = ncol(filtx))
  }, mc.cores = cv.cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  if (!is.null(rownames(x))) {
    rownames(output) <- unlist(lapply(predslist, rownames))}
  summary <- predSummary(output)
  fit.roc <- NULL
  if (!reg & nlevels(y) == 2) {
    fit.roc <- pROC::roc(output$testy, output$predyp, direction = "<",
                         quiet = TRUE)
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
  out <- list(call = outercv.call,
              output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              outer_folds = outer_folds,
              dimx = dim(x),
              final_fit = fit,
              final_vars = colnames(filtx),
              roc = fit.roc,
              summary = summary)
  class(out) <- "outercv"
  out
}


#' @rdname outercv
#' @importFrom stats na.fail model.frame model.response reformulate terms
#' @export
#' 
outercv.formula <- function(formula, data,
                            model,
                            outer_method = c("cv", "LOOCV"),
                            n_outer_folds = 10,
                            cv.cores = 1, ...,
                            na.action = na.fail) {
  outercv.call <- match.call(expand.dots = TRUE)
  # if model does not use formula, then revert to outercv.default(x, y, ...)
  if (!"formula" %in% formalArgs(model)) {
    # sample code from randomForest.formula/ svm.formula/ train.formula
    m <- match.call(expand.dots = FALSE)
    names(m)[2] <- "formula"
    if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    y <- model.response(m)
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    attr(y, "na.action") <- attr(m, "na.action")
    m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)),
                     data.frame(m))
    out <- outercv.default(y, m, outer_method = outer_method, 
                           n_outer_folds = n_outer_folds, cv.cores = cv.cores, ...)
    return(out)
  }
  # for models designed for formula method
  outer_method <- match.arg(outer_method)
  y <- data[, all.vars(formula[[2]])]
  reg <- !(is.factor(y) | is.character(y))  # y = regression
  outer_folds <- switch(outer_method,
                        cv = createFolds(y, k = n_outer_folds),
                        LOOCV = 1:length(y))
  outer_res <- mclapply(1:length(outer_folds), function(i) {
    test <- outer_folds[[i]]
    fit <- model(formula = formula, data = data, ...)
    # test on outer CV
    predy <- predict(fit, newdata = data[test, ])
    preds <- data.frame(predy=predy, testy=y[test])
    # for AUC
    if (!reg & nlevels(y) == 2) {
      predyp <- predict(fit, newdata = data[test, ], type = "prob")
      predyp <- predyp[,2]
      preds$predyp <- predyp
    }
    rownames(preds) <- rownames(data[test, , drop = FALSE])
    list(preds = preds,
         fit = fit)
  }, mc.cores = cv.cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  rownames(output) <- unlist(lapply(predslist, rownames))
  summary <- predSummary(output)
  fit.roc <- NULL
  if (!reg & nlevels(y) == 2) {
    fit.roc <- pROC::roc(output$testy, output$predyp, direction = "<",
                         quiet = TRUE)
  }
  
  # fit final model
  fit <- model(formula = formula, data = data, ...)
  out <- list(call = outercv.call,
              output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              outer_folds = outer_folds,
              dimx = c(nrow(data), length(labels(terms(fit)))),
              final_fit = fit,
              roc = fit.roc,
              summary = summary)
  class(out) <- "outercv"
  out
}


#' @export
summary.outercv <- function(object, 
                                 digits = max(3L, getOption("digits") - 3L), 
                                 ...) {
  cat("Single cross-validation to measure performance\n")
  cat("Outer loop: ", switch(object$outer_method,
                             cv = paste0(length(object$outer_folds), "-fold CV"),
                             LOOCV = "leave-one-out CV"))
  cat("\nNo inner loop\n")
  cat(object$dimx[1], "observations,", object$dimx[2], "predictors\n")
  cat("Model: ", object$call$model, "\n")
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
  if (hasMethod2(object$final_fit, "print")) {
    cat("\nFinal fit:")
    print(object$final_fit)
  }
  cat("\nResult:\n")
  print(object$summary, digits = digits, print.gap = 2L)
  out <- list(dimx = object$dimx, nfilter = nfilter,
              result = object$summary)
  invisible(out)
}


#' @importFrom utils methods
#' 
hasMethod2 <- function(object, func) {
  func %in% unlist(lapply(class(object), function(x) {
        attr(methods(class = x), "info")$generic}))
}
