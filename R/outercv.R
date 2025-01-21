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
#' @param x Matrix or dataframe of predictors
#' @param formula A formula describing the model to be fitted
#' @param data A matrix or data frame containing variables in the model.
#' @param model Character value or function of the model to be fitted.
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter].
#'   Any function can be provided and is passed `y` and `x`. Ideally returns a
#'   numeric vector with indices of filtered predictors. The custom function can
#'   return a character vector of names of the filtered predictors, but this
#'   will not work with the `penalty.factor` argument in [nestcv.glmnet()]. Not
#'   available if `outercv` is called with a formula.
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
#' @param modifyX Character string specifying the name of a function to modify
#'   `x`. This can be an imputation function for replacing missing values, or a
#'   more complex function which alters or even adds columns to `x`. The
#'   required return value of this function depends on the `modifyX_useY`
#'   setting.
#' @param modifyX_useY Logical value whether the `x` modifying function makes
#'   use of response training data from `y`. If `FALSE` then the `modifyX`
#'   function simply needs to return a modified `x` object. If `TRUE` then the
#'   `modifyX` function must return a model type object on which `predict()` can
#'   be called, so that train and test partitions of `x` can be modified
#'   independently.
#' @param modifyX_options List of additional arguments passed to the `x`
#'   modifying function
#' @param outer_method String of either `"cv"` or `"LOOCV"` specifying whether
#'   to do k-fold CV or leave one out CV (LOOCV) for the outer folds
#' @param n_outer_folds Number of outer CV folds
#' @param outer_folds Optional list containing indices of test folds for outer
#'   CV. If supplied, `n_outer_folds` is ignored.
#' @param parallel_method parallelisation options "mclapply" (default) or "future".
#' @param cv.cores Number of cores for parallel processing of the outer loops.
#' If parallel_method="future" cv.cores will be ignored 
#'   for backward-compatibilty of [future::plan()].
#' @param multicore_fork Logical whether to use forked multicore parallel
#'   processing. Forked multicore processing uses `parallel::mclapply`. It is
#'   only available on unix/mac as windows does not allow forking. It is set to
#'   `FALSE` by default in windows and `TRUE` in unix/mac. Non-forked parallel
#'   processing is executed using `parallel::parLapply` or `pbapply::pblapply`
#'   if `verbose` is `TRUE`. If parallel_method="future" multicore_fork will be ignored 
#'   for backward-compatibilty of [future::plan()].

#' @param predict_type Only used with binary classification. Calculation of ROC
#'   AUC requires predicted class probabilities from fitted models. Most model
#'   functions use syntax of the form `predict(..., type = "prob")`. However,
#'   some models require a different `type` to be specified, which can be passed
#'   to `predict()` via `predict_type`.
#' @param outer_train_predict Logical whether to save predictions on outer
#'   training folds to calculate performance on outer training folds.
#' @param returnList Logical whether to return list of results after main outer
#'   CV loop without concatenating results. Useful for debugging.
#' @param final Logical whether to fit final model.
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
#' @param verbose Logical whether to print messages and show progress
#' @param suppressMsg Logical whether to suppress messages and printed output
#'   from model functions. This is necessary when using forked multicore 
#'   parallelisation.
#' @param ... Optional arguments passed to the function specified by `model`.
#' @return An object with S3 class "outercv"
#'   \item{call}{the matched call}
#'   \item{output}{Predictions on the left-out outer folds}
#'   \item{outer_result}{List object of results from each outer fold containing
#'   predictions on left-out outer folds, model result and number of filtered
#'   predictors at each fold.}
#'   \item{dimx}{vector of number of observations and number of predictors}
#'   \item{outer_folds}{List of indices of outer test folds}
#'   \item{final_fit}{Final fitted model on whole data}
#'   \item{final_vars}{Column names of filtered predictors entering final model}
#'   \item{roc}{ROC AUC for binary classification where available.}
#'   \item{summary}{Overall performance summary. Accuracy and balanced accuracy
#'   for classification. ROC AUC for binary classification. RMSE for
#'   regression.}
#' @details 
#'   Some predictive model functions do not have an x & y interface. If the
#'   function specified by `model` requires a formula, `x` & `y` will be merged
#'   into a dataframe with `model()` called with a formula equivalent to 
#'   `y ~ .`.
#'   
#'   The S3 formula method for `outercv` is not really recommended with large
#'   data sets - it is envisaged to be primarily used to compare
#'   performance of more basic models e.g. `lm()` specified by formulae for
#'   example incorporating interactions. NOTE: filtering is not available if
#'   `outercv` is called with a formula - use the `x-y` interface instead.
#'   
#'   An alternative method of tuning a single model with fixed parameters
#'   is to use [nestcv.train] with `tuneGrid` set as a single row of a
#'   data.frame. The parameters which are needed for a specific model can be
#'   identified using [caret::modelLookup()].
#'   
#'   Case weights can be passed to model function which accept these, however
#'   `outercv` assumes that these are passed to the model via an argument named
#'   `weights`.
#'
#'   Note that in the case of `model = "lm"`, although additional arguments e.g.
#'   `subset`, `weights`, `offset` are passed into the model function via
#'   `"..."` the scoping is known to go awry. Avoid using these arguments with
#'   `model = "lm"`.
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
#' cvfit <- outercv(y2, x, "randomForest")
#' summary(cvfit)
#' plot(cvfit$roc)
#' 
#' ## Mixture discriminant analysis (MDA)
#' if (requireNamespace("mda", quietly = TRUE)) {
#'   library(mda)
#'   cvfit <- outercv(y2, x, "mda", predict_type = "posterior")
#'   summary(cvfit)
#' }
#' 
#' 
#' ## Example with continuous outcome
#' y <- -3 + 0.5 * dt$marker1 + 2 * dt$marker2 + rnorm(nrow(dt), 0, 2)
#' dt$outcome <- y
#' 
#' ## simple linear model - formula interface
#' cvfit <- outercv(outcome ~ ., data = dt, model = "lm")
#' summary(cvfit)
#' 
#' ## random forest for regression
#' cvfit <- outercv(y, x, "randomForest")
#' summary(cvfit)
#' 
#' ## example with lm_filter() to reduce input predictors
#' cvfit <- outercv(y, x, "randomForest", filterFUN = lm_filter,
#'                  filter_options = list(nfilter = 2, p_cutoff = NULL))
#' summary(cvfit)
#' 
#' @export
#' 

outercv <- function(y, ...) {
  UseMethod("outercv")
}

#' @rdname outercv
#' @export
#' @importFrom future.apply future_lapply
#' 
outercv.default <- function(y, x,
                            model,
                            filterFUN = NULL,
                            filter_options = NULL,
                            weights = NULL,
                            balance = NULL,
                            balance_options = NULL,
                            modifyX = NULL,
                            modifyX_useY = FALSE,
                            modifyX_options = NULL,
                            outer_method = c("cv", "LOOCV"),
                            n_outer_folds = 10,
                            outer_folds = NULL,
                            parallel_method="mclapply",
                            cv.cores = 1,
                            multicore_fork = (Sys.info()["sysname"] != "Windows"),
                            predict_type = "prob",
                            outer_train_predict = FALSE,
                            returnList = FALSE,
                            final = TRUE,
                            na.option = "pass",
                            verbose = FALSE,
                            suppressMsg = verbose,
                            ...) {

   if ((!missing(cv.cores) | !missing(multicore_fork)) & parallel_method=="future") {
    warning("When parallel_method is future, cv.cores and multicore_fork arguments will be ignored for backward-compatibilty")
  }
  
  start <- Sys.time()
  outercv.call <- match.call(expand.dots = TRUE)
  if (!exists(substitute(model))) {
    stop("model '", substitute(model), "' not found", call. = FALSE)
  }
  ok <- checkxy(y, x, na.option, weights)
  y <- y[ok$r]
  x <- x[ok$r, ok$c, drop = FALSE]
  weights <- weights[ok$r]
  if (!is.null(balance) & !is.null(weights)) {
    stop("`balance` and `weights` cannot be used at the same time")}
  reg <- !(is.factor(y) | is.character(y))  # y = regression
  if (!is.null(balance) & reg) {
    stop("`balance` can only be used for classification")}
  if (!is.null(weights) & !checkArg("weights", model)) {
    warning("'weights' argument not found in model function")
  }
  outer_method <- match.arg(outer_method)
  if (is.null(outer_folds)) {
    outer_folds <- switch(outer_method,
                          cv = createFolds(y, k = n_outer_folds),
                          LOOCV = 1:length(y))
  }
  if (outercv.call$model == "glm") predict_type <- "response"
  if (outercv.call$model == "mda") predict_type <- "posterior"
  
  verbose <- as.numeric(verbose)
  if(parallel_method=="mclapply"){
  if (verbose == 1 && (!multicore_fork || Sys.getenv("RSTUDIO") == "1")) {
    message("Performing ", n_outer_folds, "-fold outer CV, using ",
            plural(cv.cores, "core(s)"))}
  
  dots <- list(...)
  if (!multicore_fork && cv.cores >= 2) {
    cl <- makeCluster(cv.cores)
    varlist <- c("outer_folds", "y", "x", "model", "reg","filterFUN", 
                 "filter_options", "weights", "balance", "balance_options",
                 "modifyX", "modifyX_useY", "modifyX_options",
                 "predict_type", "outer_train_predict", "outercvCore", 
                 "suppressMsg", "dots")
    clusterExport(cl, varlist = varlist, envir = environment())
    on.exit(stopCluster(cl))
    if (verbose) {
      # pblapply
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("Package 'pbapply' must be installed", call. = FALSE)}
      outer_res <- pbapply::pblapply(seq_along(outer_folds), function(i) {
        args <- c(list(i=i, y=y, x=x, outer_folds=outer_folds, model=model, reg=reg,
                       filterFUN=filterFUN, filter_options=filter_options,
                       weights=weights, balance=balance,
                       balance_options=balance_options,
                       modifyX=modifyX, modifyX_useY=modifyX_useY,
                       modifyX_options=modifyX_options,
                       predict_type=predict_type,
                       outer_train_predict=outer_train_predict,
                       suppressMsg=suppressMsg), dots)
        do.call(outercvCore, args)
      }, cl = cl)
    } else {
      # parLapply
      outer_res <- parLapply(cl = cl, seq_along(outer_folds), function(i) {
        args <- c(list(i=i, y=y, x=x, outer_folds=outer_folds, model=model, reg=reg,
                       filterFUN=filterFUN, filter_options=filter_options,
                       weights=weights, balance=balance,
                       balance_options=balance_options,
                       modifyX=modifyX, modifyX_useY=modifyX_useY,
                       modifyX_options=modifyX_options,
                       predict_type=predict_type,
                       outer_train_predict=outer_train_predict,
                       suppressMsg=suppressMsg), dots)
        do.call(outercvCore, args)
      })
    }
  } else {
    # linux/mac, forked
    outer_res <- mclapply(seq_along(outer_folds), function(i) {
      outercvCore(i, y, x, outer_folds, model, reg,
                  filterFUN, filter_options, weights,
                  balance, balance_options,
                  modifyX, modifyX_useY, modifyX_options,
                  predict_type,
                  outer_train_predict, verbose, suppressMsg, ...)
    }, mc.cores = cv.cores)
  }
}else{
     if (verbose == 1) {
    message("Performing ", n_outer_folds, "-fold outer CV")
  }

  dots <- list(...)
  outer_res <- future_lapply(seq_along(outer_folds), function(i) {
      outercvCore(i, y, x, outer_folds, model, reg,
                  filterFUN, filter_options, weights,
                  balance, balance_options,
                  modifyX, modifyX_useY, modifyX_options,
                  predict_type,
                  outer_train_predict, verbose, suppressMsg, ...)
  }, future.seed = TRUE)
}

  if (returnList) return(outer_res)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  if (!is.null(rownames(x))) {
    rownames(output) <- unlist(lapply(predslist, rownames))}
  if (outercv.call$model == "glm") {
    output$predy <- levels(output$testy)[as.numeric(output$predyp > 0.5) +1]
  }
  summary <- predSummary(output)
  fit.roc <- NULL
  if (!reg & nlevels(y) == 2) {
    fit.roc <- pROC::roc(output$testy, output$predyp, direction = "<",
                         quiet = TRUE)
  }
  
  # fit final model
  if (final) {
    if (verbose == 1) message("Fitting final model on whole data")
    dat <- nest_filt_bal(NULL, y, x, filterFUN, filter_options,
                         balance, balance_options,
                         modifyX, modifyX_useY, modifyX_options)
    yfinal <- dat$ytrain
    filtx <- dat$filt_xtrain
    
    if ("formula" %in% formalArgs(model)) {
      dat <- if (is.data.frame(filtx)) {filtx
      } else as.data.frame(filtx, stringsAsFactors = TRUE)
      dat$.outcome <- yfinal
      if (is.null(weights)) {
        args <- c(list(as.formula(".outcome ~ ."), data = quote(dat)), dots)
      } else {
        args <- c(list(as.formula(".outcome ~ ."), data = quote(dat),
                       weights = quote(weights)), dots)
      }
    } else {
      if (is.null(weights)) {
        args <- c(alist(y = yfinal, x = filtx), dots)  # prevents evaluation
      } else {
        args <- c(alist(y = yfinal, x = filtx, weights = weights), dots)
      }
    }
    if (suppressMsg) {
      printlog <- capture.output({ fit <- do.call(model, args) })
    } else fit <- do.call(model, args)
  } else yfinal <- fit <- filtx <- NULL
  
  end <- Sys.time()
  if (verbose == 1) message("Duration: ", format(end - start))
  out <- list(call = outercv.call,
              output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              outer_folds = outer_folds,
              dimx = dim(x),
              y = y,
              yfinal = yfinal,
              final_fit = fit,
              final_vars = colnames(filtx),
              roc = fit.roc,
              summary = summary)
  if (final && !is.null(modifyX)) {
    out$xfinal <- filtx
    if (modifyX_useY) out$modify_fit <- dat$modify_fit
  }
  class(out) <- "outercv"
  out
}


outercvCore <- function(i, y, x, outer_folds, model, reg,
                        filterFUN, filter_options, weights,
                        balance, balance_options,
                        modifyX, modifyX_useY, modifyX_options,
                        predict_type,
                        outer_train_predict, verbose = FALSE,
                        suppressMsg = FALSE,
                        ...) {
  start <- Sys.time()
  test <- outer_folds[[i]]
  dat <- nest_filt_bal(test, y, x, filterFUN, filter_options,
                       balance, balance_options,
                       modifyX, modifyX_useY, modifyX_options)
  ytrain <- dat$ytrain
  ytest <- dat$ytest
  filt_xtrain <- dat$filt_xtrain
  filt_xtest <- dat$filt_xtest
  dots <- list(...)
  
  # check if model uses formula
  if ("formula" %in% formalArgs(model)) {
    dat <- if (is.data.frame(filt_xtrain)) {filt_xtrain
    } else as.data.frame(filt_xtrain, stringsAsFactors = TRUE)
    dat$.outcome <- ytrain
    if (is.null(weights)) {
      args <- c(list(as.formula(".outcome ~ ."), data = quote(dat)), dots)
    } else {
      args <- c(list(as.formula(".outcome ~ ."), data = quote(dat),
                   weights = quote(weights[-test])), dots)
    }
  } else {
    if (is.null(weights)) {
      args <- c(alist(y = ytrain, x = filt_xtrain), dots)
    } else {
      args <- c(alist(y = ytrain, x = filt_xtrain, weights = weights[-test]),
                dots)
    }
    
  }
  
  if (suppressMsg) {
    printlog <- capture.output({ fit <- do.call(model, args) })
  } else fit <- do.call(model, args)
  
  # test on outer CV
  predy <- predict(fit, newdata = filt_xtest)
  preds <- data.frame(predy=predy, testy=ytest)
  # for AUC
  if (!reg & nlevels(y) == 2) {
    predyp <- predict(fit, newdata = filt_xtest, type = predict_type)
    if (!is.vector(predyp)) predyp <- predyp[,2]
    preds$predyp <- predyp
  }
  rownames(preds) <- rownames(filt_xtest)
  if (outer_train_predict) {
    predy <- predict(fit, newdata = filt_xtrain)
    train_preds <- data.frame(ytrain=ytrain, predy=predy)
    rownames(train_preds) <- rownames(filt_xtrain)
    # for AUC
    if (!reg & nlevels(y) == 2) {
      predyp <- predict(fit, newdata = filt_xtrain, type = predict_type)
      if (!is.vector(predyp)) predyp <- predyp[,2]
      train_preds <- cbind(train_preds, predyp)
    }
  } else train_preds <- NULL
  
  if (verbose == 1) {
    end <- Sys.time()
    message_parallel("Fitted fold ", i, " (", format(end - start, digits = 3), ")")
  } else if (verbose == 2) cat_parallel("=")
  
  list(preds = preds,
       train_preds = train_preds,
       fit = fit,
       nfilter = ncol(filt_xtest),
       ytrain = ytrain)
}


#' @rdname outercv
#' @importFrom stats na.fail model.frame model.response reformulate terms
#' @export
#' 
outercv.formula <- function(formula, data,
                            model,
                            outer_method = c("cv", "LOOCV"),
                            n_outer_folds = 10,
                            outer_folds = NULL,
                            parallel_method="mclapply",
                            cv.cores = 1,
                            multicore_fork = (Sys.info()["sysname"] != "Windows"),
                            predict_type = "prob",
                            outer_train_predict = FALSE,
                            verbose = FALSE,
                            suppressMsg = verbose,
                            ..., na.action = na.fail) {

   if ((!missing(cv.cores) | !missing(multicore_fork)) & parallel_method=="future") {
    warning("When parallel_method is future, cv.cores and multicore_fork arguments will be ignored for backward-compatibilty")
  }
 
  start <- Sys.time()
  outercv.call <- match.call(expand.dots = TRUE)
  if (!exists(substitute(model))) {
    stop("model '", substitute(model), "' not found", call. = FALSE)
  }
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
    if(parallel_method=="mclapply"){
      out <- outercv.default(y, m, outer_method = outer_method,n_outer_folds = n_outer_folds, cv.cores = cv.cores, ...)
     }else{
      out <- outercv.default(y, m, outer_method = outer_method, n_outer_folds = n_outer_folds, ...)
     }
    return(out)
  }
  # for models designed for formula method
  verbose <- as.numeric(verbose)
  if(parallel_method=="mclapply"){
    if(verbose == 1 && (!multicore_fork || Sys.getenv("RSTUDIO") == "1")) {
        message("Performing ", n_outer_folds, "-fold outer CV, using ",plural(cv.cores, "core(s)"))
    }
    }else{
       if (verbose == 1) { message("Performing ", n_outer_folds, "-fold outer CV")}
    }

  outer_method <- match.arg(outer_method)
  y <- data[, all.vars(formula[[2]])]
  reg <- !(is.factor(y) | is.character(y))  # y = regression
  if (is.null(outer_folds)) {
    outer_folds <- switch(outer_method,
                          cv = createFolds(y, k = n_outer_folds),
                          LOOCV = 1:length(y))
  }
  if (outercv.call$model == "glm") predict_type <- "response"
  
  dots <- list(...)
  if(parallel_method=="mclapply"){
  if (Sys.info()["sysname"] == "Windows" & cv.cores >= 2) {
    cl <- makeCluster(cv.cores)
    clusterExport(cl, varlist = c("outer_folds", "formula", "data", "y", 
                                  "model", "reg", "predict_type",
                                  "outer_train_predict", "outercvFormulaCore",
                                  "suppressMsg", "dots"),
                  envir = environment())
    on.exit(stopCluster(cl))
    outer_res <- parLapply(cl = cl, seq_along(outer_folds), function(i) {
      args <- c(list(i=i, outer_folds=outer_folds, formula=formula, data=data, 
                     y=y, model=model, reg=reg, predict_type=predict_type,
                     suppressMsg=suppressMsg), 
                dots)
      do.call(outercvFormulaCore, args)
    })
  } else {
    outer_res <- mclapply(seq_along(outer_folds), function(i) {
      outercvFormulaCore(i, outer_folds, formula, data, y, model,
                         reg, predict_type, outer_train_predict, verbose,
                         suppressMsg, ...)
    }, mc.cores = cv.cores)
  }
}else{
    outer_res <- future_lapply(seq_along(outer_folds), function(i) {
      outercvFormulaCore(i, outer_folds, formula, data, y, model,
                         reg, predict_type, outer_train_predict, verbose,
                         suppressMsg, ...)
  }, future.seed = TRUE)
}
  
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  rownames(output) <- unlist(lapply(predslist, rownames))
  if (outercv.call$model == "glm") {
    output$predy <- levels(output$testy)[as.numeric(output$predyp > 0.5) +1]
  }
  summary <- predSummary(output)
  fit.roc <- NULL
  if (!reg & nlevels(y) == 2) {
    fit.roc <- pROC::roc(output$testy, output$predyp, direction = "<",
                         quiet = TRUE)
  }
  
  # fit final model
  if (verbose == 1) message("Fitting final model on whole data")
  args <- c(list(formula = formula, data = quote(data)), dots)
  if (suppressMsg) {
    printlog <- capture.output({ fit <- do.call(model, args) })
  } else fit <- do.call(model, args)
  
  end <- Sys.time()
  if (verbose == 1) message("Duration: ", format(end - start))
  
  out <- list(call = outercv.call,
              output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              outer_folds = outer_folds,
              dimx = c(nrow(data), length(labels(terms(fit)))),
              final_fit = fit,
              y = y,
              roc = fit.roc,
              summary = summary)
  class(out) <- "outercv"
  out
}

  
outercvFormulaCore <- function(i, outer_folds, formula, data, y, model,
                               reg, predict_type, outer_train_predict,
                               verbose = FALSE, suppressMsg = FALSE, ...) {
  start <- Sys.time()
  test <- outer_folds[[i]]
  dots <- list(...)
  args <- c(list(formula = formula, data = quote(data[-test, ])), dots)
  if (suppressMsg) {
    printlog <- capture.output({ fit <- do.call(model, args) })
  } else fit <- do.call(model, args)
  
  # test on outer CV
  predy <- predict(fit, newdata = data[test, ])
  preds <- data.frame(predy=predy, testy=y[test])
  # for AUC
  if (!reg & nlevels(y) == 2) {
    predyp <- predict(fit, newdata = data[test, ], type = predict_type)
    if (!is.vector(predyp)) predyp <- predyp[,2]
    preds$predyp <- predyp
  }
  rownames(preds) <- rownames(data)[test]
  if (outer_train_predict) {
    predy <- predict(fit, newdata = data[test, ])
    train_preds <- data.frame(ytrain=y[-test], predy=predy)
    # for AUC
    if (!reg & nlevels(y) == 2) {
      predyp <- predict(fit, newdata = data[-test, ], type = predict_type)
      if (!is.vector(predyp)) predyp <- predyp[,2]
      train_preds <- cbind(train_preds, predyp)
    }
    rownames(train_preds) <- rownames(data)[-test]
  } else train_preds <- NULL
  
  if (verbose == 1) {
    end <- Sys.time()
    message_parallel("Fitted fold ", i, " (", format(end - start, digits = 3), ")")
  } else if (verbose == 2) cat_parallel("=")
  
  list(preds = preds,
       train_preds = train_preds,
       fit = fit)
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
  if (!is.null(object$call$modifyX))
    cat("Modifier: ", object$call$modifyX, "\n")
  balance <- object$call$balance
  if (!is.null(balance)) {
    cat("Balancing: ", balance, "\n")
  }
  cat(object$dimx[1], "observations,", object$dimx[2], "predictors\n")
  if (!is.numeric(object$y)) print(c(table(object$y)))
  cat("\n")
  
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
  print(object$summary, digits = digits, print.gap = 3L)
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
