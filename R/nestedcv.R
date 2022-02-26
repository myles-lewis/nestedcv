# nestedcv
# by Myles Lewis
# 24-02-2022

#' Univariate filter
#' 
#' Simple univariate filter using t-test or anova.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param return Type of vector returned. Default "names" returns predictor 
#' names, "full" returns a named vector of p values.
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a named vector of p values is returned.
#' @importFrom Rfast ttests
#' @export
#' 
uni_filter <- function(y,
                       x,
                       nfilter = NULL,
                       p_cutoff = 0.05,
                       return = "names") {
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  res <- Rfast::ttests(x[indx1, ], x[indx2, ])
  rownames(res) <- colnames(x)
  if (return == "full") return(res[, 'stat'])
  out <- res[, 'pvalue']
  out <- sort(out[out < p_cutoff])
  if (!is.null(nfilter)) out <- out[1:nfilter]
  names(out)
}

#' glmnet coefficients
#' 
#' Convenience function for retrieving coefficients from a [cv.glmnet] model at 
#' a specified lambda. Sparsity is removed and non-intercept coefficients are 
#' ranked by absolute value.
#' 
#' @param fit A [cv.glmnet] fitted model object.
#' @param s Value of lambda. See [coef.glmnet] and [predict.cv.glmnet]
#' @return Vector of coefficients ordered with the intercept first, followed 
#' by highest absolute value to lowest.
#' @importFrom stats coef
#' @export
#' 
glmnet_coefs <- function(fit, s) {
  cf <- coef(fit, s = s)
  cf <- as.matrix(cf)
  cf <- cf[cf != 0, ]
  cf2 <- cf[-1]
  cf2 <- cf2[order(abs(cf2), decreasing = TRUE)]
  c(cf[1], cf2)  # keep intercept first
}

#' Random forest filter
#' 
#' Fits a random forest model and ranks variables by variable importance. 
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param return Type of vector returned. Default "names" returns predictor 
#' names, "full" returns a named vector of variable importance.
#' @param ntree Number of trees to grow. See [randomForest].
#' @param mtry Number of predictors randomly sampled as candidates at each 
#' split. See [randomForest].
#' @param ... Optional arguments passed to [randomForest].
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a named vector of variable importance is returned.
#' @details
#' This filter uses the [randomForest] function from the randomForest package.
#' Variable importance is calculated using the [importance] function, specifying
#' type 1 = mean decrease in accuracy. See [importance].
#' @importFrom randomForest randomForest importance
#' @export
#' 
rf_filter <- function(y, x, nfilter = NULL, return = "names",
                      ntree = 1000,
                      mtry = ncol(x) * 0.2,
                      ...) {
  fit <- randomForest::randomForest(x, y, importance = TRUE,
                                    ntree = ntree, mtry = mtry, ...)
  vi <- as.vector(importance(fit, type = 1))
  names(vi) <- colnames(x)
  if (return == "full") return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) vi <- vi[1:min(nfilter, length(vi))]
  names(vi)
}

#' ReliefF filter
#' 
#' Uses ReliefF algorithm from the CORElearn package to rank predictors in order 
#' of importance.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param estimator Type of algorithm used, see [CORElearn::attrEval]
#' @param return Type of vector returned. Default "names" returns predictor 
#' names, "full" returns a named vector of variable importance.
#' @param ... Other arguments passed to [CORElearn::attrEval]
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a named vector of variable importance is returned.
#' @importFrom CORElearn attrEval
#' @export
#'
relieff_filter <- function(y, x, nfilter = NULL, 
                           estimator = "ReliefFequalK",
                           return = "names", ...) {
  df <- as.data.frame(x)
  df$y <- y
  ref <- CORElearn::attrEval('y', df, estimator = estimator, ...)
  names(ref) <- colnames(x)
  if (return == "full") return(ref)
  ref <- sort(ref, decreasing = TRUE)
  if (!is.null(nfilter)) ref <- ref[1:min(nfilter, length(ref))]
  names(ref)
}

#' Combo filter
#' 
#' Filter combining univariate filtering and reliefF filtering in equal measure.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return, using 1/2 from `uni_filter` 
#' and 1/2 from `relieff_filter`. Since `unique` is applied, the final number 
#' returned may be less than `nfilter`.
#' @param return Type of output returned. Default "names" returns predictor 
#' names, "full" returns full output.
#' @param ... Optional arguments passed via [relieff_filter] to 
#' [CORElearn::attrEval]
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a list of 2 vectors of full output from both [uni_filter] 
#' and [relieff_filter] is returned.
#' @export
#' 
combo_filter <- function(y, x, nfilter, return = "names", ...) {
  uni_set <- uni_filter(y, x, nfilter, return = return)
  relf_set <- relieff_filter(y, x, nfilter, return = return, ...)
  if (return == "full") {
    return(list(unifilt = uni_set, relieff_filter = relf_set))
  }
  n <- round(nfilter / 2)
  unique(c(uni_set[1:n], relf_set[1:n]))
}

#' Nested cross-validation with glmnet
#' 
#' Nested cross-validation (CV) with glmnet including tuning of elastic net 
#' alpha parameter and embedding of a filter function within the nested CV.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param family Either a character string representing one of the built-in 
#' families, or else a `glm()` family object. Passed to [cv.glmnet] and [glmnet]
#' @param filterFUN Filter function, e.g. [uni_filter] or [relieff_filter]. 
#' Any function can be provided and is passed `y` and `x`. Must return a 
#' character vector with names of filtered predictors.
#' @param n_outer_folds Number of outer CV folds
#' @param n_inner_folds Number of inner CV folds
#' @param alphaSet Vector of alphas to be tuned
#' @param min_1se Value from 0 to 1 specifying choice of optimal lambda from 
#' 0=lambda.min to 1=lambda.1se
#' @param keep Logical indicating whether inner CV predictions are 
#' retained for calculating left-out inner CV fold accuracy etc. See argument 
#' `keep` in [cv.glmnet].
#' @param cores Number of cores for parallel processing. Note this currently 
#' uses [parallel::mclapply].
#' @param ... Optional arguments passed both to [cv.glmnet] as well as the 
#' filter function specified by `filterFUN`
#' @return An object with S3 class "nestcv.glmnet"
#' @importFrom caret createFolds confusionMatrix defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom parallel mclapply
#' @importFrom stats predict setNames
#' @export
#' 
nestcv.glmnet <- function(y, x,
                          family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
                          filterFUN = NULL,
                          n_outer_folds = 10,
                          n_inner_folds = 10,
                          alphaSet = seq(0, 1, 0.1),
                          min_1se = 0,
                          keep = TRUE,
                          cores = 1,
                          ...) {
  
  outer_folds <- createFolds(y, k = n_outer_folds, returnTrain = TRUE)
  outer_res <- mclapply(1:n_outer_folds, function(i) {
    trainIndex <- outer_folds[[i]]
    foldid <- sample(rep(seq_len(n_inner_folds), length = length(trainIndex)))
    # expand data with interactions
    filtx <- if (is.null(filterFUN)) x else {
      fset <- filterFUN(y[trainIndex], x[trainIndex, ], ...)
      x[, fset]
    }
    fit <- lapply(alphaSet, function(alpha) {
      cv.glmnet(x = filtx[trainIndex, ], y = y[trainIndex], 
                alpha = alpha, nfolds = n_inner_folds, foldid = foldid, 
                keep = keep, family = family, ...)
    })
    alphas <- unlist(lapply(fit, function(fitx) {
      w <- which.min(fitx$cvm)
      fitx$cvm[w]
    }))
    alpha.x <- which.min(alphas)
    s <- exp((log(fit[[alpha.x]]$lambda.min) * (1-min_1se) + log(fit[[alpha.x]]$lambda.1se) * min_1se))
    cf <- as.matrix(coef(fit[[alpha.x]], s = s))
    cf <- cf[cf != 0, ]
    # test on outer CV
    predy <- as.vector(predict(fit[[alpha.x]], newx = filtx[-trainIndex, ], s = s, type = "class"))
    predyp <- as.vector(predict(fit[[alpha.x]], newx = filtx[-trainIndex, ], s = s))
    preds <- data.frame(predy=predy, predyp=predyp, testy=y[-trainIndex])
    rownames(preds) <- rownames(x[-trainIndex, ])
    ret <- list(preds = preds,
                alpha = alphaSet[[alpha.x]], lambda = s, coef = cf, cvfit = fit[[alpha.x]],
                cv_alpha = alphas,
                nfilter = ncol(filtx))
    # inner CV predictions
    if (keep) {
      ind <- fit[[alpha.x]]$index["min", ]
      innerCV_preds <- fit[[alpha.x]]$fit.preval[, ind]
      ytrain <- y[trainIndex]
      ret <- append(ret, list(ytrain = ytrain, innerCV_preds = innerCV_preds))
    }
    ret
  }, mc.cores = cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  
  if (family == "binomial") {
    cm <- table(output$predy, output$testy)
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    glmnet.roc <- pROC::roc(output$testy, output[, 2], direction = "<")
    auc <- glmnet.roc$auc
    summary <- setNames(c(auc, acc, b_acc), c("AUC", "Accuracy", "Bal_accuracy"))
  } else if (family == "multinomial") {
    cm <- table(output$predy, output$testy)
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    summary <- setNames(c(acc, b_acc), c("Accuracy", "Bal_accuracy"))
  } else {
    df <- data.frame(obs = output$testy, pred = output$predy)
    summary <- caret::defaultSummary(df)
  }
  
  # fit final glmnet
  lam <- mean(unlist(lapply(outer_res, '[[', 'lambda')))
  alph <- mean(unlist(lapply(outer_res, '[[', 'alpha')))
  filtx <- if (is.null(filterFUN)) x else {
    fset <- filterFUN(y, x, ...)
    x[, fset]
  }
  fit <- glmnet(filtx, y, alpha = alph, family = family, ...)
  out <- list(output = output,
              outer_result = outer_res,
              mean_lambda = lam,
              mean_alpha = alph,
              final_fit = fit,
              roc = glmnet.roc,
              alphaSet = alphaSet,
              summary = summary)
  class(out) <- "nestcv.glmnet"
  out
}

#' Plot cross-validated glmnet alpha
#' 
#' Plot of cross-validated glmnet alpha parameter against deviance.
#' 
#' @param x Fitted "nestcv.glmnet" object
#' @param col Optional vector of line colours for each fold
#' @param ... other arguments passed to plot
#' @return No return value
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @export
#' 
plot_alphas <- function(x,
                        col = NULL,
                        ...) {
  cv_alpha <- lapply(x$outer_result, '[[', 'cv_alpha')
  n <- length(cv_alpha)
  if (is.null(col)) {
    col <- rainbow(n)
  } else {
    col <- rep_len(col, n)
  }
  new.args <- list(...)
  plot.args <- list(y = cv_alpha[[1]], x = x$alphaSet,
       type = 'l',
       ylim = range(unlist(cv_alpha)),
       xlab = 'Alpha',
       ylab = x$outer_result[[1]]$cvfit$name,
       col = col[1])
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call(plot, plot.args)
  for (i in 2:n) {
    lines(cv_alpha[[i]], x = x$alphaSet, col = col[i])
  }
}

#' Build ROC curve from left-out folds from inner CV
#' 
#' Build ROC (receiver operating characteristic) curve from left-out folds 
#' from inner CV. Object can be plotted or passed to functions [auc] etc.
#' 
#' @param cva Fitted `"nestcv.glmnet"` object 
#' @param direction Passed to [pROC::roc]
#' @param ... Other arguments passed to [pROC::roc]
#' @return `"roc"` object, see [pROC::roc]
#' @importFrom pROC roc auc
#' @export
#' 
innercv_roc <- function(cva, direction = "<", ...) {
  innerpreds <- unlist(lapply(cva$outer_result, '[[', 'innerCV_preds'))
  ytrain <- unlist(lapply(cva$outer_result, '[[', 'ytrain'))
  pROC::roc(ytrain, innerpreds, direction = direction, ...)
}

