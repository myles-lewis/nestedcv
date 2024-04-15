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
#' @param y Response vector or matrix. Matrix is only used for 
#' `family = 'mgaussian'` or `'cox'`.
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
#' @param modifyX Character string specifying the name of a function to modify
#'   `x`. This can be an imputation function for replacing missing values, or a
#'   more complex function which alters or even adds columns to `x`. The
#'   required return value of this function depends on the `modifyX_useY`
#'   setting.
#' @param modifyX_useY Logical value whether the `x` modifying function makes
#'   use of response training data from `y`. If `FALSE` then the `modifyX`
#'   function simply needs to return a modified `x` object, which will be
#'   coerced to a matrix as required by `glmnet`. If `TRUE` then the `modifyX`
#'   function must return a model type object on which `predict()` can be
#'   called, so that train and test partitions of `x` can be modified
#'   independently.
#' @param modifyX_options List of additional arguments passed to the `x`
#'   modifying function
#' @param outer_method String of either `"cv"` or `"LOOCV"` specifying whether
#'   to do k-fold CV or leave one out CV (LOOCV) for the outer folds
#' @param n_outer_folds Number of outer CV folds
#' @param n_inner_folds Number of inner CV folds
#' @param outer_folds Optional list containing indices of test folds for outer
#'   CV. If supplied, `n_outer_folds` is ignored.
#' @param pass_outer_folds Logical indicating whether the same outer folds are
#'   used for fitting of the final model when final CV is applied. Note this can
#'   only be applied when `n_outer_folds` and `n_inner_folds` are the same and
#'   no balancing is applied.
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
#'   variables. See [glmnet]. Note this works separately from filtering. For
#'   some `nestedcv` filter functions you might need to set `force_vars` to
#'   avoid filtering out features.
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
#' @param verbose Logical whether to print messages and show progress
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
                          modifyX = NULL,
                          modifyX_useY = FALSE,
                          modifyX_options = NULL,
                          outer_method = c("cv", "LOOCV"),
                          n_outer_folds = 10,
                          n_inner_folds = 10,
                          outer_folds = NULL,
                          pass_outer_folds = FALSE,
                          alphaSet = seq(0.1, 1, 0.1),
                          min_1se = 0,
                          keep = TRUE,
                          outer_train_predict = FALSE,
                          weights = NULL,
                          penalty.factor = rep(1, ncol(x)),
                          cv.cores = 1,
                          finalCV = TRUE,
                          na.option = "omit",
                          verbose = FALSE,
                          ...) {
  start <- Sys.time()
  family <- match.arg(family)
  nestcv.call <- match.call(expand.dots = TRUE)
  outer_method <- match.arg(outer_method)
  if (is.character(y)) y <- factor(y)
  if (is.factor(y) && !family %in% c("binomial", "multinomial"))
    stop("`y` is not numeric: incorrect `family`")
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
  ok <- checkxy(y, x, na.option, weights)
  y <- if (is.matrix(y)) y[ok$r, ] else y[ok$r]
  x <- x[ok$r, ok$c, drop = FALSE]
  weights <- weights[ok$r]
  if (!is.null(balance) && !is.null(weights)) {
    stop("`balance` and `weights` cannot be used at the same time")}
  if (!is.null(balance) && is.numeric(y)) {
    stop("`balance` can only be used for classification")}
  
  if (is.null(outer_folds)) {
    y1 <- if (is.matrix(y)) y[,1] else y
    outer_folds <- switch(outer_method,
                          cv = createFolds(y1, k = n_outer_folds),
                          LOOCV = 1:NROW(y))
  } else {
    if ("n_outer_folds" %in% names(nestcv.call)) {
      if (n_outer_folds != length(outer_folds))
        stop("Mismatch between n_outer_folds and length(outer_folds)")
    }
    n_outer_folds <- length(outer_folds)
  }
  
  verbose <- as.numeric(verbose)
  if (verbose == 1) message("Performing ", n_outer_folds, "-fold outer CV, using ",
                       plural(cv.cores, "core(s)"))
  if (Sys.info()["sysname"] == "Windows" & cv.cores >= 2) {
    cl <- makeCluster(cv.cores)
    dots <- list(...)
    varlist = c("outer_folds", "y", "x", "filterFUN", "filter_options",
                "alphaSet", "min_1se",  "n_inner_folds", "keep", "family",
                "weights", "balance", "balance_options", "penalty.factor",
                "modifyX", "modifyX_useY", "modifyX_options",
                "outer_train_predict", "nestcv.glmnetCore", "dots")
    clusterExport(cl, varlist = varlist, envir = environment())
    on.exit(stopCluster(cl))
    if (verbose == 1) {
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("Package 'pbapply' must be installed", call. = FALSE)}
      outer_res <- pbapply::pblapply(seq_along(outer_folds), function(i) {
        args <- c(list(i=i, y=y, x=x, outer_folds=outer_folds,
                       filterFUN=filterFUN, filter_options=filter_options,
                       balance=balance, balance_options=balance_options,
                       modifyX=modifyX, modifyX_useY=modifyX_useY,
                       modifyX_options=modifyX_options,
                       alphaSet=alphaSet, min_1se=min_1se,
                       n_inner_folds=n_inner_folds, keep=keep, family=family,
                       weights=weights, penalty.factor=penalty.factor,
                       outer_train_predict=outer_train_predict), dots)
        do.call(nestcv.glmnetCore, args)
      }, cl = cl)
    } else {
      outer_res <- parLapply(cl = cl, seq_along(outer_folds), function(i) {
        args <- c(list(i=i, y=y, x=x, outer_folds=outer_folds,
                       filterFUN=filterFUN, filter_options=filter_options,
                       balance=balance, balance_options=balance_options,
                       modifyX=modifyX, modifyX_useY=modifyX_useY,
                       modifyX_options=modifyX_options,
                       alphaSet=alphaSet, min_1se=min_1se,
                       n_inner_folds=n_inner_folds, keep=keep, family=family,
                       weights=weights, penalty.factor=penalty.factor,
                       outer_train_predict=outer_train_predict), dots)
        do.call(nestcv.glmnetCore, args)
      })
    }
  } else {
    # linux/mac
    outer_res <- mclapply(seq_along(outer_folds), function(i) {
      nestcv.glmnetCore(i, y, x, outer_folds, filterFUN, filter_options,
                        balance, balance_options,
                        modifyX, modifyX_useY, modifyX_options,
                        alphaSet, min_1se, n_inner_folds, keep, family,
                        weights, penalty.factor, outer_train_predict,
                        verbose, ...)
    }, mc.cores = cv.cores)
  }
  
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  if (!is.null(rownames(x))) {
    rownames(output) <- unlist(lapply(predslist, rownames))}
  
  summary <- predSummary(output, family = family)
  glmnet.roc <- NULL
  if (family == "binomial") {
    glmnet.roc <- pROC::roc(output$testy, output$predyp, direction = "<", 
                           quiet = TRUE)
  }
  
  if (is.na(finalCV)) {
    fit <- final_coef <- final_param <- yfinal <- final_vars <- xsub <- filtx <- NA
  } else {
    dat <- nest_filt_bal(NULL, y, x, filterFUN, filter_options,
                         balance, balance_options,
                         modifyX, modifyX_useY, modifyX_options,
                         penalty.factor = penalty.factor)
    yfinal <- dat$ytrain
    filtx <- as.matrix(dat$filt_xtrain)
    filtpen.factor <- dat$filt_pen.factor
    
    if (finalCV) {
      # use CV on whole data to finalise parameters
      if (verbose == 1)
        message("Fitting final model using ", n_inner_folds, "-fold CV on whole data")
      foldid <- NULL
      if (pass_outer_folds) {
        if (n_outer_folds == n_inner_folds && is.null(balance)) {
          foldid <- rep(0, NROW(y))
          for (i in 1:length(outer_folds)) {
            foldid[outer_folds[[i]]] <- i
          }
        } else message("Cannot pass `outer_folds` to final CV")
      }
      cvafit <- cva.glmnet(filtx, yfinal, alphaSet = alphaSet, family = family,
                           weights = weights, penalty.factor = filtpen.factor, 
                           nfolds = n_inner_folds, foldid = foldid, ...)
      alphafit <- cvafit$fits[[cvafit$which_alpha]]
      s <- exp((log(alphafit$lambda.min) * (1-min_1se) + log(alphafit$lambda.1se) * min_1se))
      fit <- cvafit$fits[[cvafit$which_alpha]]
      final_param <- setNames(c(s, cvafit$best_alpha), c("lambda", "alpha"))
    } else {
      # use outer folds for final parameters
      if (verbose == 1) message("Fitting final model based on outer CV parameters")
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
      cfmean <- try(colMeans(x[, names(fin_coef)[-1], drop = FALSE], na.rm = TRUE),
                    silent = TRUE)
      if (inherits(cfmean, "try-error")) cfmean <- rep(NA, length(fin_coef) -1)
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
    all_vars <- all_vars[all_vars %in% colnames(x)]
    xsub <- x[, all_vars]
  }
  
  end <- Sys.time()
  if (verbose == 1) message("Duration: ", format(end - start))
  
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
  if (!is.na(finalCV) & !is.null(modifyX)) {
    out$xfinal <- filtx
    if (modifyX_useY) out$modify_fit <- dat$modify_fit
  }
  class(out) <- "nestcv.glmnet"
  out
}


nestcv.glmnetCore <- function(i, y, x, outer_folds, filterFUN, filter_options,
                              balance, balance_options,
                              modifyX, modifyX_useY, modifyX_options,
                              alphaSet, min_1se, n_inner_folds, keep, family,
                              weights, penalty.factor,
                              outer_train_predict, verbose = FALSE, ...) {
  start <- Sys.time()
  test <- outer_folds[[i]]
  dat <- nest_filt_bal(test, y, x, filterFUN, filter_options,
                       balance, balance_options,
                       modifyX, modifyX_useY, modifyX_options,
                       penalty.factor)
  ytrain <- dat$ytrain
  ytest <- dat$ytest
  filt_xtrain <- as.matrix(dat$filt_xtrain)
  filt_xtest <- as.matrix(dat$filt_xtest)
  filt_pen.factor <- dat$filt_pen.factor
  
  cvafit <- cva.glmnet(x = filt_xtrain, y = ytrain, 
                       alphaSet = alphaSet, nfolds = n_inner_folds,
                       keep = keep, family = family, weights = weights[-test],
                       penalty.factor = filt_pen.factor, ...)
  alphafit <- cvafit$fits[[cvafit$which_alpha]]
  s <- exp((log(alphafit$lambda.min) * (1-min_1se) + log(alphafit$lambda.1se) * min_1se))
  cf <- glmnet_coefs(alphafit, s = s)
  # test on outer CV
  if (family == "mgaussian") {
    # mgaussian
    predy <- predict(alphafit, newx = filt_xtest, s = s)[,, 1]
    preds <- as.data.frame(cbind(ytest, predy))
    colnames(preds)[1:ncol(y)] <- paste0("ytest.", colnames(ytest))
  } else if (family == "cox") {
    # cox
    predy <- as.vector(predict(alphafit, newx = filt_xtest, s = s))
    preds <- as.data.frame(cbind(ytest, predy))
  } else {
    # default
    predy <- as.vector(predict(alphafit, newx = filt_xtest, s = s, type = "class"))
    preds <- data.frame(testy=ytest, predy=predy)
  }
  if (family == "binomial") {
    predyp <- as.vector(predict(alphafit, newx = filt_xtest, s = s))
    preds <- cbind(preds, predyp)
  } else if (family == "multinomial") {
    # glmnet generates 3d array
    predyp <- predict(alphafit, newx = filt_xtest, s = s)[,, 1]
    preds <- cbind(preds, predyp)
  }
  if (outer_train_predict) {
    if (is.matrix(y)) {
      # mgaussian, cox
      train_predy <- predict(alphafit, newx = filt_xtrain, s = s)
      train_preds <- as.data.frame(cbind(ytrain, train_predy))
      if (family == "mgaussian") {
        colnames(train_preds)[1:ncol(y)] <- paste0("ytrain.", colnames(ytrain))
      }
    } else {
      train_predy <- as.vector(predict(alphafit, newx = filt_xtrain, s = s, type = "class"))
      train_preds <- data.frame(ytrain=ytrain, predy=train_predy)
    }
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
    innerCV_preds <- if (length(dim(alphafit$fit.preval)) == 3) {
      alphafit$fit.preval[, , ind]
    } else alphafit$fit.preval[, ind]
    ret <- append(ret, list(innerCV_preds = innerCV_preds))
  }
  if (verbose == 1) {
    end <- Sys.time()
    message_parallel("Fitted fold ", i, " (", format(end - start, digits = 3), ")")
  } else if (verbose == 2) cat_parallel("=")
  ret
}
