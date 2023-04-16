# nestedcv
# by Myles Lewis
# 03-03-2022


#' Nested cross-validation with glmnet
#'
#' This function enables nested cross-validation (CV) with glmnet including
#' tuning of elastic net alpha parameter. The function also allows the option of
#' embedded filtering of predictors for feature selection nested within the
#' outer loop of CV. Predictions on the outer test folds are brought back
#' together and error estimation/ accuracy determined. The default is 10x10
#' nested CV.
#'
#' @param y Response vector
#' @param x Matrix of predictors. Dataframes will be coerced to a matrix as
#'   is necessary for glmnet.
#' @param family Either a character string representing one of the built-in
#'   families, or else a `glm()` family object. Passed to [cv.glmnet] and
#'   [glmnet]
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter].
#'   Any function can be provided and is passed `y` and `x`. Must return a
#'   character vector with names of filtered predictors.
#' @param filter_options List of additional arguments passed to the filter
#'   function specified by `filterFUN`.
#' @param balance Specifies method for dealing with imbalanced class data.
#'   Current options are `"randomsample"` or `"smote"`. See [randomsample()] and
#'   [smote()]
#' @param balance_options List of additional arguments passed to the balancing
#'   function
#' @param outer_method String of either `"cv"` or `"LOOCV"` specifying whether
#'   to do k-fold CV or leave one out CV (LOOCV) for the outer folds
#' @param n_outer_folds Number of outer CV folds
#' @param n_inner_folds Number of inner CV folds
#' @param outer_folds Optional list containing indices of test folds for outer
#'   CV. If supplied, `n_outer_folds` is ignored.
#' @param alphaSet Vector of alphas to be tuned
#' @param min_1se Value from 0 to 1 specifying choice of optimal lambda from
#'   0=lambda.min to 1=lambda.1se
#' @param keep Logical indicating whether inner CV predictions are retained for
#'   calculating left-out inner CV fold accuracy etc. See argument `keep` in
#'   [cv.glmnet].
#' @param outer_train_predict Logical whether to save predictions on outer
#'   training folds to calculate performance on outer training folds.
#' @param weights Weights applied to each sample. Note `weights` and `balance`
#'   cannot be used at the same time. Weights are only applied in glmnet and not
#'   in filters.
#' @param penalty.factor Separate penalty factors can be applied to each
#'   coefficient. Can be 0 for some variables, which implies no shrinkage, and
#'   that variable is always included in the model. Default is 1 for all
#'   variables. See [glmnet]
#' @param cv.cores Number of cores for parallel processing of the outer loops.
#'   NOTE: this uses `parallel::mclapply` on unix/mac and `parallel::parLapply`
#'   on windows.
#' @param finalCV Logical whether to perform one last round of CV on the whole
#'   dataset to determine the final model parameters. If set to `FALSE`, the
#'   median of hyperparameters from outer CV folds are used for the final model.
#'   Performance metrics are independent of this last step. If set to `NA`,
#'   final model fitting is skipped altogether, which gives a useful speed boost
#'   if performance metrics are all that is needed.
#' @param na.option Character value specifying how `NA`s are dealt with.
#'   `"omit"` (the default) is equivalent to `na.action = na.omit`. `"omitcol"`
#'   removes cases if there are `NA` in 'y', but columns (predictors) containing
#'   `NA` are removed from 'x' to preserve cases. Any other value means that
#'   `NA` are ignored (a message is given).
#' @param ... Optional arguments passed to [cv.glmnet]
#' @return An object with S3 class "nestcv.glmnet"
#'   \item{call}{the matched call}
#'   \item{output}{Predictions on the left-out outer folds}
#'   \item{outer_result}{List object of results from each outer fold containing
#'   predictions on left-out outer folds, best lambda, best alpha, fitted glmnet
#'   coefficients, list object of inner fitted cv.glmnet and number of filtered
#'   predictors at each fold.}
#'   \item{outer_method}{the `outer_method` argument}
#'   \item{n_inner_folds}{number of inner folds}
#'   \item{outer_folds}{List of indices of outer test folds}
#'   \item{dimx}{dimensions of `x`}
#'   \item{xsub}{subset of `x` containing all predictors used in both outer CV
#'   folds and the final model}
#'   \item{y}{original response vector}
#'   \item{yfinal}{final response vector (post-balancing)}
#'   \item{final_param}{Final mean best lambda
#'   and alpha from each fold}
#'   \item{final_fit}{Final fitted glmnet model}
#'   \item{final_coef}{Final model coefficients and mean expression. Variables
#'   with coefficients shrunk to 0 are removed.}
#'   \item{final_vars}{Column names of filtered predictors entering final model.
#'   This is useful for subsetting new data for predictions.}
#'   \item{roc}{ROC AUC for binary classification where available.}
#'   \item{summary}{Overall performance summary. Accuracy and balanced accuracy
#'   for classification. ROC AUC for binary classification. RMSE for
#'   regression.}
#' @details
#' glmnet does not tolerate missing values, so `na.option = "omit"` is the
#' default.
#' @author Myles Lewis
#' @importFrom caret createFolds confusionMatrix defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom parallel mclapply makeCluster clusterExport stopCluster parLapply
#' @importFrom pROC roc
#' @importFrom Rfast colmeans
#' @importFrom stats predict setNames
#' @examples
#' \donttest{
#' ## Example binary classification problem with P >> n
#' x <- matrix(rnorm(150 * 2e+04), 150, 2e+04)  # predictors
#' y <- factor(rbinom(150, 1, 0.5))  # binary response
#' 
#' ## Partition data into 2/3 training set, 1/3 test set
#' trainSet <- caret::createDataPartition(y, p = 0.66, list = FALSE)
#' 
#' ## t-test filter using whole dataset
#' filt <- ttest_filter(y, x, nfilter = 100)
#' filx <- x[, filt]
#' 
#' ## Train glmnet on training set only using filtered predictor matrix
#' library(glmnet)
#' fit <- cv.glmnet(filx[trainSet, ], y[trainSet], family = "binomial")
#' plot(fit)
#' 
#' ## Predict response on test partition
#' predy <- predict(fit, newx = filx[-trainSet, ], s = "lambda.min", type = "class")
#' predy <- as.vector(predy)
#' predyp <- predict(fit, newx = filx[-trainSet, ], s = "lambda.min", type = "response")
#' predyp <- as.vector(predyp)
#' output <- data.frame(testy = y[-trainSet], predy = predy, predyp = predyp)
#' 
#' ## Results on test partition
#' ## shows bias since univariate filtering was applied to whole dataset
#' predSummary(output)
#' 
#' ## Nested CV
#' ## n_outer_folds reduced to speed up example
#' fit2 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1,
#'                       n_outer_folds = 3,
#'                       filterFUN = ttest_filter,
#'                       filter_options = list(nfilter = 100),
#'                       cv.cores = 2)
#' summary(fit2)
#' plot_lambdas(fit2, showLegend = "bottomright")
#' 
#' ## ROC plots
#' library(pROC)
#' testroc <- roc(output$testy, output$predyp, direction = "<")
#' inroc <- innercv_roc(fit2)
#' plot(fit2$roc)
#' lines(inroc, col = 'blue')
#' lines(testroc, col = 'red')
#' legend('bottomright', legend = c("Nested CV", "Left-out inner CV folds", 
#'                                  "Test partition, non-nested filtering"), 
#'        col = c("black", "blue", "red"), lty = 1, lwd = 2, bty = "n")
#' }
#' @export
#' 
nestcv.glmnet <- function(y, x,
                          family = c("gaussian", "binomial", "poisson", 
                                     "multinomial", "cox", "mgaussian"),
                          filterFUN = NULL,
                          filter_options = NULL,
                          balance = NULL,
                          balance_options = NULL,
                          outer_method = c("cv", "LOOCV"),
                          n_outer_folds = 10,
                          n_inner_folds = 10,
                          outer_folds = NULL,
                          alphaSet = seq(0, 1, 0.1),
                          min_1se = 0,
                          keep = TRUE,
                          outer_train_predict = FALSE,
                          weights = NULL,
                          penalty.factor = rep(1, ncol(x)),
                          cv.cores = 1,
                          finalCV = TRUE,
                          na.option = "omit",
                          ...) {
  family <- match.arg(family)
  nestcv.call <- match.call(expand.dots = TRUE)
  outer_method <- match.arg(outer_method)
  if (is.character(y)) y <- factor(y)
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
  ok <- checkxy(y, x, na.option, weights)
  y <- y[ok$r]
  x <- x[ok$r, ok$c]
  weights <- weights[ok$r]
  if (!is.null(balance) & !is.null(weights)) {
    stop("`balance` and `weights` cannot be used at the same time")}
  if (!is.null(balance) & is.numeric(y)) {
    stop("`balance` can only be used for classification")}
  
  if (is.null(outer_folds)) {
    outer_folds <- switch(outer_method,
                          cv = createFolds(y, k = n_outer_folds),
                          LOOCV = 1:length(y))
  }
  
  if (Sys.info()["sysname"] == "Windows" & cv.cores >= 2) {
    cl <- makeCluster(cv.cores)
    dots <- list(...)
    varlist = c("outer_folds", "y", "x", "filterFUN", "filter_options",
                "alphaSet", "min_1se",  "n_inner_folds", "keep", "family",
                "weights", "balance", "balance_options", "penalty.factor",
                "outer_train_predict", "nestcv.glmnetCore", "dots")
    clusterExport(cl, varlist = varlist, envir = environment())
    on.exit(stopCluster(cl))
    outer_res <- parLapply(cl = cl, outer_folds, function(test) {
      args <- c(list(test=test, y=y, x=x, filterFUN=filterFUN,
                     filter_options=filter_options,
                     balance=balance, balance_options=balance_options,
                     alphaSet=alphaSet, min_1se=min_1se,
                     n_inner_folds=n_inner_folds, keep=keep, family=family,
                     weights=weights, penalty.factor=penalty.factor,
                     outer_train_predict=outer_train_predict), dots)
      do.call(nestcv.glmnetCore, args)
    })
  } else {
    outer_res <- mclapply(outer_folds, function(test) {
      nestcv.glmnetCore(test, y, x, filterFUN, filter_options,
                        balance, balance_options,
                        alphaSet, min_1se, n_inner_folds, keep, family,
                        weights, penalty.factor, outer_train_predict, ...)
    }, mc.cores = cv.cores)
  }
  
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  if (!is.null(rownames(x))) {
    rownames(output) <- unlist(lapply(predslist, rownames))}
  
  summary <- predSummary(output)
  glmnet.roc <- NULL
  if (family == "binomial") {
    glmnet.roc <- pROC::roc(output$testy, output$predyp, direction = "<", 
                           quiet = TRUE)
  }
  
  if (is.na(finalCV)) {
    fit <- final_coef <- final_param <- yfinal <- final_vars <- xsub <- NA
  } else {
    dat <- nest_filt_bal(NULL, y, x, filterFUN, filter_options,
                         balance, balance_options,
                         penalty.factor = penalty.factor)
    yfinal <- dat$ytrain
    filtx <- dat$filt_xtrain
    filtpen.factor <- dat$filt_pen.factor
    
    if (finalCV) {
      # use CV on whole data to finalise parameters
      cvafit <- cva.glmnet(filtx, yfinal, alphaSet = alphaSet, family = family,
                           weights = weights, penalty.factor = filtpen.factor, ...)
      alphafit <- cvafit$fits[[cvafit$which_alpha]]
      s <- exp((log(alphafit$lambda.min) * (1-min_1se) + log(alphafit$lambda.1se) * min_1se))
      fit <- cvafit$fits[[cvafit$which_alpha]]
      final_param <- setNames(c(s, cvafit$best_alpha), c("lambda", "alpha"))
    } else {
      # use outer folds for final parameters
      lam <- exp(median(log(unlist(lapply(outer_res, '[[', 'lambda')))))
      alph <- median(unlist(lapply(outer_res, '[[', 'alpha')))
      final_param <- setNames(c(lam, alph), c("lambda", "alpha"))
      fit <- glmnet(filtx, yfinal, alpha = alph, family = family, 
                    weights = weights, penalty.factor = filtpen.factor, ...)
    }
    
    fin_coef <- glmnet_coefs(fit, s = final_param["lambda"])
    if (is.list(fin_coef) | length(fin_coef) == 1) {
      final_coef <- fin_coef  # multinomial
    } else {
      cfmean <- colmeans(x[, names(fin_coef)[-1], drop = FALSE])
      final_coef <- data.frame(coef = fin_coef, meanExp = c(NA, cfmean))
    }
    final_vars <- colnames(filtx)
    # collect all vars from outer_res
    all_vars <- unlist(lapply(outer_res, function(i) {
      cf <- i$coef
      if (!is.list(cf)) return(names(cf)[-1])
      # multinomial
      unlist(lapply(cf, function(j) {
        names(j)[-1]
      }))
    }))
    all_vars <- unique(c(all_vars, final_vars))
    xsub <- x[, all_vars]
  }
  out <- list(call = nestcv.call,
              output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              n_inner_folds = n_inner_folds,
              outer_folds = outer_folds,
              dimx = dim(x),
              xsub = xsub,
              y = y,
              yfinal = yfinal,
              final_param = final_param,
              final_fit = fit,
              final_coef = final_coef,
              final_vars = final_vars,
              roc = glmnet.roc,
              summary = summary)
  class(out) <- "nestcv.glmnet"
  out
}


nestcv.glmnetCore <- function(test, y, x, filterFUN, filter_options,
                              balance, balance_options,
                              alphaSet, min_1se, n_inner_folds, keep, family,
                              weights, penalty.factor,
                              outer_train_predict, ...) {
  dat <- nest_filt_bal(test, y, x, filterFUN, filter_options,
                       balance, balance_options, penalty.factor)
  ytrain <- dat$ytrain
  ytest <- dat$ytest
  filt_xtrain <- dat$filt_xtrain
  filt_xtest <- dat$filt_xtest
  filt_pen.factor <- dat$filt_pen.factor
  
  cvafit <- cva.glmnet(x = filt_xtrain, y = ytrain, 
                       alphaSet = alphaSet, nfolds = n_inner_folds,
                       keep = keep, family = family, weights = weights[-test],
                       penalty.factor = filt_pen.factor, ...)
  alphafit <- cvafit$fits[[cvafit$which_alpha]]
  s <- exp((log(alphafit$lambda.min) * (1-min_1se) + log(alphafit$lambda.1se) * min_1se))
  cf <- glmnet_coefs(alphafit, s = s)
  # test on outer CV
  predy <- as.vector(predict(alphafit, newx = filt_xtest, s = s, type = "class"))
  preds <- data.frame(testy=ytest, predy=predy)
  if (family == "binomial") {
    predyp <- as.vector(predict(alphafit, newx = filt_xtest, s = s))
    preds <- cbind(preds, predyp)
  } else if (family == "multinomial") {
    # glmnet generates 3d array
    predyp <- predict(alphafit, newx = filt_xtest, s = s)[,, 1]
    preds <- cbind(preds, predyp)
  }
  if (outer_train_predict) {
    train_predy <- as.vector(predict(alphafit, newx = filt_xtrain, s = s, type = "class"))
    train_preds <- data.frame(ytrain=ytrain, predy=train_predy)
    if (family == "binomial") {
      predyp <- as.vector(predict(alphafit, newx = filt_xtrain, s = s))
      train_preds <- cbind(train_preds, predyp)
    } else if (family == "multinomial") {
      # glmnet generates 3d array
      predyp <- predict(alphafit, newx = filt_xtrain, s = s)[,, 1]
      train_preds <- cbind(train_preds, predyp)
    }
  } else train_preds <- NULL
  rownames(preds) <- rownames(filt_xtest)
  ret <- list(preds = preds,
              train_preds = train_preds,
              lambda = s,
              alpha = cvafit$best_alpha,
              coef = cf,
              cvafit = cvafit,
              nfilter = ncol(filt_xtrain),
              ytrain = ytrain)
  # inner CV predictions
  if (keep) {
    ind <- alphafit$index["min", ]
    innerCV_preds <- if (family == "multinomial") {
      alphafit$fit.preval[, , ind]
    } else alphafit$fit.preval[, ind]
    ret <- append(ret, list(innerCV_preds = innerCV_preds))
  }
  ret
}


#' Cross-validation of alpha for glmnet
#' 
#' Performs k-fold cross-validation for glmnet, including alpha mixing parameter.
#' 
#' @param x Matrix of predictors
#' @param y Response vector
#' @param nfolds Number of folds (default 10)
#' @param alphaSet Sequence of alpha values to cross-validate
#' @param ... Other arguments passed to [cv.glmnet]
#' @return Object of S3 class "cva.glmnet", which is a list of the cv.glmnet 
#' objects for each value of alpha and `alphaSet`.
#' \item{fits}{List of fitted [cv.glmnet] objects}
#' \item{alphaSet}{Sequence of alpha values used}
#' \item{alpha_cvm}{The mean cross-validated error - a vector of length 
#' `length(alphaSet)`.}
#' \item{best_alpha}{Value of alpha giving lowest `alpha_cvm`.}
#' \item{which_alpha}{Index of `alphaSet` with lowest `alpha_cvm`}
#' @seealso [cv.glmnet], [glmnet]
#' @author Myles Lewis
#' @importFrom glmnet cv.glmnet
#' @importFrom utils tail
#' @export
#' 
cva.glmnet <- function(x, y, nfolds = 10, alphaSet = seq(0.1, 1, 0.1), ...) {
  foldid <- sample(rep(seq_len(nfolds), length = length(y)))
  fit1 <- cv.glmnet(x = x, y = y, 
                    alpha = tail(alphaSet, 1), foldid = foldid, ...)
  if (length(alphaSet) > 1) {
    fits <- lapply(alphaSet[1:(length(alphaSet)-1)], function(alpha) {
      cv.glmnet(x = x, y = y, 
                alpha = alpha, foldid = foldid, lambda = fit1$lambda, ...)
    })
    fits <- append(fits, list(fit1))
  } else fits <- list(fit1)
  alpha_cvm <- unlist(lapply(fits, function(i) min(i$cvm)))
  which_alpha <- which.min(alpha_cvm)
  best_alpha <- alphaSet[which_alpha]
  cvafit <- list(fits = fits,
                  alphaSet = alphaSet,
                  alpha_cvm = alpha_cvm,
                  best_alpha = best_alpha,
                  which_alpha = which_alpha)
  class(cvafit) <- "cva.glmnet"
  cvafit
}


#' Extract coefficients from a cva.glmnet object
#' 
#' Extracts model coefficients from a fitted [cva.glmnet()] object.
#' 
#' @param object Fitted `cva.glmnet` object.
#' @param ... Other arguments passed to `coef.glmnet()` e.g. `s` the value of
#'   lambda at which coefficients are required.
#' @export
coef.cva.glmnet <- function(object, ...) {
  coef(object$fits[[object$which_alpha]], ...)
}


#' glmnet coefficients
#' 
#' Convenience function for retrieving coefficients from a [cv.glmnet] model at 
#' a specified lambda. Sparsity is removed and non-intercept coefficients are 
#' ranked by absolute value.
#' 
#' @param fit A [cv.glmnet] fitted model object.
#' @param s Value of lambda. See [coef.glmnet] and [predict.cv.glmnet]
#' @param ... Other arguments passed to [coef.glmnet]
#' @return Vector or list of coefficients ordered with the intercept first, 
#' followed by highest absolute value to lowest.
#' @importFrom stats coef
#' @export
#' 
glmnet_coefs <- function(fit, s, ...) {
  if (length(fit) == 1 && is.na(fit)) return(NA)
  cf <- coef(fit, s = s, ...)
  if (is.list(cf)) {
    cf <- lapply(cf, function(i) {
      cf <- as.matrix(i)
      cf <- cf[cf != 0, ]
      cf2 <- cf[-1]
      cf2 <- cf2[order(abs(cf2), decreasing = TRUE)]
      c(cf[1], cf2)
    })
    return(cf)
  } 
  cf <- as.matrix(cf)
  cf <- cf[cf != 0, ]
  cf2 <- cf[-1]
  cf2 <- cf2[order(abs(cf2), decreasing = TRUE)]
  c(cf[1], cf2)  # keep intercept first
}


#' Extract coefficients from nestcv.glmnet object
#' 
#' Extracts coefficients from the final fit of a `"nestcv.glmnet"` object.
#' 
#' @param object Object of class `"nestcv.glmnet"`
#' @param s Value of penalty parameter lambda. Default is the mean of lambda 
#' values selected across each outer fold.
#' @param ... Other arguments passed to [coef.glmnet]
#' @return Vector or list of coefficients ordered with the intercept first, 
#' followed by highest absolute value to lowest.
#' @export
#'
coef.nestcv.glmnet <- function(object, s = object$final_param["lambda"], ...) {
  glmnet_coefs(object$final_fit, s = s, ...)
}


#' @export
print.nestcv.glmnet <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Nested cross-validation with glmnet\n")
  if (!is.null(x$call$filterFUN)) 
    cat("Filter: ", x$call$filterFUN, "\n") else cat("No filter\n")
  cat("\nFinal parameters:\n")
  if (length(x$final_param)==1 && is.na(x$final_param)) {
    cat("NA\n")
  } else print(x$final_param, digits = digits, print.gap = 2L)
  cat("\nFinal coefficients:\n")
  if (length(x$final_fit)==1 && is.na(x$final_fit)) {
    cat("NA\n")
  } else print(coef(x), digits = digits)
  cat("\nResult:\n")
  print(x$summary, digits = digits, print.gap = 2L)
}


#' @export
summary.nestcv.glmnet <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Nested cross-validation with glmnet\n")
  if (!is.null(object$call$filterFUN)) 
    cat("Filter: ", object$call$filterFUN, "\n") else cat("No filter\n")
  cat("Outer loop: ", switch(object$outer_method,
                         cv = paste0(length(object$outer_folds), "-fold CV"),
                         LOOCV = "leave-one-out CV"))
  cat("\nInner loop: ", paste0(object$n_inner_folds, "-fold CV\n"))
  balance <- object$call$balance
  if (!is.null(balance)) {
    cat("Balancing: ", balance, "\n")
  }
  cat(object$dimx[1], "observations,", object$dimx[2], "predictors\n")
  if (!is.numeric(object$y)) print(c(table(object$y)))
  cat("\n")
  
  alpha <- unlist(lapply(object$outer_result, '[[', 'alpha'))
  lambda <- unlist(lapply(object$outer_result, '[[', 'lambda'))
  nfilter <- unlist(lapply(object$outer_result, '[[', 'nfilter'))
  foldres <- data.frame(alpha = alpha, lambda = lambda, n.filter = nfilter,
                        row.names = paste("Fold", seq_along(alpha)))
  
  print(foldres, digits = digits)
  cat("\nFinal parameters:\n")
  if (length(object$final_param)==1 && is.na(object$final_param)) {
    cat("NA\n")
  } else print(object$final_param, digits = digits, print.gap = 2L)
  cat("\nFinal coefficients:\n")
  if (length(object$final_fit)==1 && is.na(object$final_fit)) {
    cat("NA\n")
  } else print(coef(object), digits = digits)
  cat("\nResult:\n")
  print(object$summary, digits = digits, print.gap = 3L)
  out <- list(dimx = object$dimx, folds = foldres,
              final_param = object$final_param,
              coef = coef(object), result = object$summary)
  invisible(out)
}

#' Predict method for nestcv.glmnet fits
#'
#' Obtains predictions from the final fitted model from a [nestcv.glmnet]
#' object.
#' 
#' @param object Fitted `nestcv.glmnet` object
#' @param newdata New data to predict outcome on
#' @param s Value of lambda for glmnet prediction
#' @param ... Other arguments passed to `predict.glmnet`.
#' @return Object returned depends on the `...` argument passed to predict
#'   method for `glmnet` objects.
#' @details Checks for missing predictors and if these are sparse (i.e. have
#'   zero coefficients) columns of 0 are automatically added to enable
#'   prediction to proceed.
#' @seealso [glmnet::glmnet]
#' @method predict nestcv.glmnet
#' @export
predict.nestcv.glmnet <- function(object, newdata,
                                  s = object$final_param["lambda"],
                                  ...) {
  newdata <- as.matrix(newdata)
  newx <- fix_cols(object$final_fit, newdata, s = s)
  predict(object$final_fit, newx = newx, s = unname(s), ...)
}


# fills in zero coefficent columns if missing
fix_cols <- function(x, newx, s) {
  cf <- coef(x, s = s)
  # check for multinomial
  if (is.list(cf)) {
    cf <- do.call(cbind, cf)
    cf <- as.matrix(cf)
  }
  final_vars <- rownames(cf)[-1]
  cf <- cf[-1,]
  # if full subset present
  if (all(final_vars %in% colnames(newx))) return(newx[, final_vars])
  # some cols missing
  nz <- if (is.vector(cf)) {
    cf != 0
  } else apply(cf, 1, function(i) any(i != 0))  # multinomial
  nonzeros <- final_vars[nz]
  if (!all(nonzeros %in% colnames(newx))) {
    stop("Some non-zero coefficient predictors are missing")}
  zeros <- final_vars[!nz]
  if (length(zeros) == 0) return(newx[, final_vars])
  miss <- zeros[!zeros %in% colnames(newx)] 
  if (length(miss) == 0) return(newx[, final_vars])
  mat <- matrix(0, nrow = nrow(newx), ncol = length(miss),
                dimnames = list(rownames(newx), miss))
  out <- cbind(newx, mat)
  out[, final_vars]  # col order seems to matter
}
