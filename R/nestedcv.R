# nestedcv
# by Myles Lewis
# 03-03-2022


#' Nested cross-validation with glmnet
#' 
#' Nested cross-validation (CV) with glmnet including tuning of elastic net 
#' alpha parameter and embedding of a filter function within the nested CV.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param family Either a character string representing one of the built-in 
#' families, or else a `glm()` family object. Passed to [cv.glmnet] and [glmnet]
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter]. 
#' Any function can be provided and is passed `y` and `x`. Must return a 
#' character vector with names of filtered predictors.
#' @param filter_options List of additional arguments passed to the filter 
#' function specified by `filterFUN`.
#' @param outer_method String of either `"cv"` or `"loocv"` specifying whether to 
#' do k-fold CV or leave one out CV (LOOCV) for the outer folds
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
#' @param ... Optional arguments passed to [cv.glmnet]
#' @return An object with S3 class "nestcv.glmnet"
#' \item{output}{Predictions on the left-out outer folds}
#' \item{outer_result}{List object of results from each outer fold containing 
#' predictions on left-out outer folds, best lambda, best alpha, fitted glmnet 
#' coefficients, list object of inner fitted cv.glmnet and number of filtered 
#' predictors at each fold.}
#' \item{outer_folds}{List of indices of outer training folds}
#' \item{mean_lambda}{Final mean best lambda from each fold}
#' \item{mean_alpha}{Final mean best alpha from each fold}
#' \item{final_fit}{Final fitted glmnet model}
#' \item{roc}{ROC AUC for binary classification where available.}
#' \item{summary}{Overall performance summary. Accuracy and balanced accuracy 
#' for classification. ROC AUC for binary classification. RMSE for regression.}
#' @importFrom caret createFolds confusionMatrix defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom parallel mclapply
#' @importFrom pROC roc
#' @importFrom stats predict setNames
#' @export
#' 
nestcv.glmnet <- function(y, x,
                          family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
                          filterFUN = NULL,
                          filter_options = NULL,
                          outer_method = c("cv", "loocv"),
                          n_outer_folds = 10,
                          n_inner_folds = 10,
                          alphaSet = seq(0, 1, 0.1),
                          min_1se = 0,
                          keep = TRUE,
                          cores = 1,
                          ...) {
  family <- match.arg(family)
  outer_method <- match.arg(outer_method)
  outer_folds <- switch(outer_method,
                        cv = createFolds(y, k = n_outer_folds),
                        loocv = 1:length(y))
  outer_res <- mclapply(1:length(outer_folds), function(i) {
    test <- outer_folds[[i]]
    # expand data with interactions
    filtx <- if (is.null(filterFUN)) x else {
      args <- list(y = y[-test], x = x[-test, ])
      args <- append(args, filter_options)
      fset <- do.call(filterFUN, args)
      x[, fset]
    }
    cvafit <- cva.glmnet(x = filtx[-test, ], y = y[-test], 
                alphaSet = alphaSet, nfolds = n_inner_folds,
                keep = keep, family = family, ...)
    alphafit <- cvafit$fits[[cvafit$which_alpha]]
    s <- exp((log(alphafit$lambda.min) * (1-min_1se) + log(alphafit$lambda.1se) * min_1se))
    cf <- as.matrix(coef(alphafit, s = s))
    cf <- cf[cf != 0, ]
    # test on outer CV
    predy <- as.vector(predict(alphafit, newx = filtx[test, ], s = s, type = "class"))
    predyp <- as.vector(predict(alphafit, newx = filtx[test, ], s = s))
    preds <- data.frame(predy=predy, predyp=predyp, testy=y[test])
    rownames(preds) <- rownames(x[test, ])
    ret <- list(preds = preds,
                lambda = s,
                alpha = cvafit$best_alpha,
                coef = cf,
                cvafit = cvafit,
                nfilter = ncol(filtx))
    # inner CV predictions
    if (keep) {
      ind <- alphafit$index["min", ]
      innerCV_preds <- alphafit$fit.preval[, ind]
      ytrain <- y[-test]
      ret <- append(ret, list(ytrain = ytrain, innerCV_preds = innerCV_preds))
    }
    ret
  }, mc.cores = cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  
  glmnet.roc <- NULL
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
    args <- list(y = y, x = x)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    x[, fset]
  }
  fit <- glmnet(filtx, y, alpha = alph, family = family, ...)
  out <- list(output = output,
              outer_result = outer_res,
              outer_method = outer_method,
              outer_folds = outer_folds,
              mean_lambda = lam,
              mean_alpha = alph,
              final_fit = fit,
              roc = glmnet.roc,
              summary = summary)
  class(out) <- "nestcv.glmnet"
  out
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
#' @importFrom glmnet cv.glmnet
#' @export
#' 
cva.glmnet <- function(x, y, nfolds = 10, alphaSet = seq(0.1, 1, 0.1), ...) {
  foldid <- sample(rep(seq_len(nfolds), length = length(y)))
  fits <- lapply(alphaSet, function(alpha) {
    cv.glmnet(x = x, y = y, 
              alpha = alpha, foldid = foldid, ...)
  })
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


#' glmnet coefficients
#' 
#' Convenience function for retrieving coefficients from a [cv.glmnet] model at 
#' a specified lambda. Sparsity is removed and non-intercept coefficients are 
#' ranked by absolute value.
#' 
#' @param fit A [cv.glmnet] fitted model object.
#' @param s Value of lambda. See [coef.glmnet] and [predict.cv.glmnet]
#' @return Vector or list of coefficients ordered with the intercept first, 
#' followed by highest absolute value to lowest.
#' @importFrom stats coef
#' @export
#' 
glmnet_coefs <- function(fit, s) {
  cf <- coef(fit, s = s)
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


#' Outer cross-validation with randomForest
#' 
#' Outer cross-validation (CV) with randomForest. Note, no tuning of 
#' parameters is performed. If tuning of parameters is required, this will 
#' require full nested CV.
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
  outer_folds <- createFolds(y, k = n_outer_folds, returnTrain = TRUE)
  outer_res <- mclapply(1:n_outer_folds, function(i) {
    trainIndex <- outer_folds[[i]]
    filtx <- if (is.null(filterFUN)) x else {
      args <- list(y = y[trainIndex], x = x[trainIndex, ])
      args <- append(args, filter_options)
      fset <- do.call(filterFUN, args)
      x[, fset]
    }
    fit <- randomForest(x = filtx[trainIndex, ], y = y[trainIndex], ...)
    # test on outer CV
    predy <- as.vector(predict(fit, newdata = filtx[-trainIndex, ], type = "response"))
    predyp <- as.vector(predict(fit, newdata = filtx[-trainIndex, ], type = "prob")[,2])
    preds <- data.frame(predy=predy, predyp=predyp, testy=y[-trainIndex])
    rownames(preds) <- rownames(x[-trainIndex, ])
    list(preds = preds,
                rf = fit,
                nfilter = ncol(filtx))
  }, mc.cores = cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  
  if (reg) {
    df <- data.frame(obs = output$testy, pred = output$predy)
    summary <- caret::defaultSummary(df)
  } else if (nlevels(y) == 2) {
    cm <- table(output$predy, output$testy)
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    rf.roc <- pROC::roc(output$testy, output[, 2], direction = "<")
    auc <- rf.roc$auc
    summary <- setNames(c(auc, acc, b_acc), c("AUC", "Accuracy", "Bal_accuracy"))
  } else {
    cm <- table(output$predy, output$testy)
    acc <- sum(diag(cm))/ sum(cm)
    ccm <- caret::confusionMatrix(cm)
    b_acc <- ccm$byClass[11]
    summary <- setNames(c(acc, b_acc), c("Accuracy", "Bal_accuracy"))
  }
  
  # fit final rf
  filtx <- if (is.null(filterFUN)) x else {
    args <- list(y = y, x = x)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    x[, fset]
  }
  fit <- randomForest(filtx, y, ...)
  out <- list(output = output,
              outer_result = outer_res,
              final_fit = fit,
              roc = rf.roc,
              summary = summary)
  class(out) <- "outercv.rf"
  out
}

#' Inner CV for random forest
#' 
#' Using caret.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfolds Number of folds for cross-validation
#' @param repeats Number of repeats for repeated CV
#' @param mtry Vector of possible values of `mtry` to tune. `mtry` is the 
#' number of predictors tried at each split, see [randomForest]
#' @param ... other arguments passed to [caret::train] or [randomForest]
#' @importFrom caret trainControl train
#' @importFrom randomForest randomForest
#' @export
#' 
cv.rf <- function(y, x,
                  nfolds = 10, repeats = 1,
                  mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
                  ...) {
  if (length(mtry) == 1) {
    return(randomForest(y = y, x = x, mtry = mtry, ...))
  }
  trCtrl <- caret::trainControl(method = "repeatedcv", number = nfolds,
                                repeats = repeats,
                                classProbs = TRUE,
                                summaryFunction = mnLogLoss)
  tuneGrid <- expand.grid(mtry = mtry)
  fit <- caret::train(x = x, y = y,
                      method = "rf",
                      metric = "logLoss",
                      trControl = trCtrl,
                      tuneGrid = tuneGrid, ...)
  fit
}


#' Nested CV for caret
#' 
#' Wrapper function for applying nested CV and predictor filtering followed by 
#' training using `caret`.
#' 
#' @param y Response vector. For classification this should be a factor.
#' @param x Matrix of predictors
#' @param filterFUN Filter function, e.g. [ttest_filter] or [relieff_filter]. 
#' Any function can be provided and is passed `y` and `x`. Must return a 
#' character vector with names of filtered predictors.
#' @param filter_options List of additional arguments passed to the filter 
#' function specified by `filterFUN`.
#' @param n_outer_folds Number of outer CV folds
#' @param metric A string that specifies what summary metric will be used to 
#' select the optimal model. By default, "logLoss" is used for classification 
#' and "RMSE" is used for regression. Note this differs from the default setting 
#' in caret which uses "Accuracy" for classification. See details.
#' @param trControl A list of values generated by the `caret` function 
#' [trainControl]. This defines how inner CV training through `caret` is 
#' performed. See http://topepo.github.io/caret/using-your-own-model-in-train.html.
#' @param tuneGrid Data frame of tuning values, see [caret::train]
#' @param savePredictions Indicates whether hold-out predictions for 
#' each inner CV fold should be saved for ROC curves, accuracy etc
#' see [caret::trainControl].Set to `"final"` to capture predictions for inner 
#' CV ROC.
#' @param cores Number of cores for parallel processing. Note this currently 
#' uses `parallel::mclapply`.
#' @param ... Arguments passed to [caret::train]
#' @details Parallelisation is performed on the outer folds using `mclapply`. 
#' For classification `metric` defaults to using 'logLoss' with the `trControl` 
#' arguments `classProbs = TRUE, summaryFunction = mnLogLoss`, rather than 
#' 'Accuracy' which is the default classification metric in `caret`. See 
#' [trainControl]. LogLoss is arguably more consistent than Accuracy for tuning 
#' parameters in datasets with small sample size.
#' @importFrom caret createFolds train trainControl mnLogLoss confusionMatrix
#' defaultSummary
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#' @importFrom pROC roc
#' @importFrom stats predict setNames
#' @export
#' 
nestcv.train <- function(y, x,
                           filterFUN = NULL,
                           filter_options = NULL, 
                           n_outer_folds = 10,
                           cores = 1,
                           metric = ifelse(is.factor(y), "logLoss", "RMSE"),
                           trControl = NULL,
                           tuneGrid = NULL,
                           savePredictions = FALSE,
                           ...) {
  if (is.null(trControl)) {
    trControl <- if (is.factor(y)) {
      trainControl(method = "repeatedcv", 
                   number = 10, repeats = 1,
                   classProbs = TRUE,
                   savePredictions = savePredictions,
                   summaryFunction = mnLogLoss)
    } else trainControl()
  }
  outer_folds <- createFolds(y, k = n_outer_folds, returnTrain = TRUE)
  outer_res <- mclapply(1:n_outer_folds, function(i) {
    trainIndex <- outer_folds[[i]]
    filtx <- if (is.null(filterFUN)) x else {
      args <- list(y = y[trainIndex], x = x[trainIndex, ])
      args <- append(args, filter_options)
      fset <- do.call(filterFUN, args)
      x[, fset]
    }
    fit <- caret::train(x = filtx[trainIndex, ], y = y[trainIndex],
                        metric = metric,
                        trControl = trControl,
                        tuneGrid = tuneGrid, ...)
    predy <- predict(fit, newdata = filtx[-trainIndex, ], type = "raw")
    predyp <- predict(fit, newdata = filtx[-trainIndex, ], type = "prob")
    # note predyp has 2 columns
    
    preds <- data.frame(predy=predy, predyp=predyp[,2], testy=y[-trainIndex])
    rownames(preds) <- rownames(x[-trainIndex, ])
    ret <- list(preds = preds,
                fit = fit,
                nfilter = ncol(filtx))
    ret
  }, mc.cores = cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  caret.roc <- NULL
  if (is.factor(y)) {
    if (nlevels(y) == 2) {
      cm <- table(output$predy, output$testy)
      acc <- sum(diag(cm))/ sum(cm)
      ccm <- caret::confusionMatrix(cm)
      b_acc <- ccm$byClass[11]
      caret.roc <- pROC::roc(output$testy, output[, 2], direction = "<")
      auc <- caret.roc$auc
      summary <- setNames(c(auc, acc, b_acc), c("AUC", "Accuracy", "Bal_accuracy"))
    } else {
      # multinomial class
      cm <- table(output$predy, output$testy)
      acc <- sum(diag(cm))/ sum(cm)
      ccm <- caret::confusionMatrix(cm)
      b_acc <- ccm$byClass[11]
      summary <- setNames(c(acc, b_acc), c("Accuracy", "Bal_accuracy"))
    }
  } else {
    # regression
    df <- data.frame(obs = output$testy, pred = output$predy)
    summary <- caret::defaultSummary(df)
  }
  bestTunes <- lapply(outer_res, function(i) i$fit$bestTune)
  bestTunes <- as.data.frame(data.table::rbindlist(bestTunes))
  rownames(bestTunes) <- paste0('Fold', 1:n_outer_folds)
  finalTune <- colMeans(bestTunes)
  finalTune <- data.frame(as.list(finalTune))
  filtx <- if (is.null(filterFUN)) x else {
    args <- list(y = y, x = x)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    x[, fset]
  }
  fitControl <- trainControl(method = "none", classProbs = TRUE)
  final_fit <- caret::train(x = filtx, y = y, 
                            trControl = fitControl,
                            tuneGrid = finalTune, ...)
  
  out <- list(output = output,
              outer_result = outer_res,
              final_fit = final_fit,
              roc = caret.roc,
              bestTunes = bestTunes,
              finalTune = finalTune,
              summary = summary)
  class(out) <- "nestcv.train"
  out
}

