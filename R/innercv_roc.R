# Extract inner CV predictions or outer training fold predictions


#' Build ROC curve from left-out folds from inner CV
#' 
#' Build ROC (receiver operating characteristic) curve from left-out folds 
#' from inner CV. Object can be plotted using `plot()` or passed to functions 
#' [auc()] etc.
#' 
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @param direction Set ROC directionality [pROC::roc]
#' @param ... Other arguments passed to [pROC::roc]
#' @return `"roc"` object, see [pROC::roc]
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
#' fit2 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1,
#'                       filterFUN = ttest_filter,
#'                       filter_options = list(nfilter = 100),
#'                       n_outer_folds = 3)
#' summary(fit2)
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
innercv_roc <- function(x, direction = "<", ...) {
  innerpreds <- innercv_preds(x)
  pROC::roc(innerpreds$testy, innerpreds$predyp, direction = direction, quiet = TRUE, ...)
}


#' Inner CV predictions
#' 
#' Obtain predictions on held-out test inner CV folds
#' 
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @return Dataframe with columns `testy` and `predy`, and for binomial and
#'   multinomial models additional columns containing probabilities or log
#'   likelihood values.
#' @export
#' 
innercv_preds <- function(x) {
  UseMethod("innercv_preds")
}


#' @rdname innercv_preds
#' @export
innercv_preds.nestcv.glmnet <- function(x) {
  ytrain <- unlist(lapply(x$outer_result, '[[', 'ytrain'))
  if (is.character(ytrain)) ytrain <- factor(ytrain)
  fam <- x$call$family
  if (is.null(fam)) fam <- ""
  if (fam == "binomial") {
    # binomial
    innerpreds <- unlist(lapply(x$outer_result, '[[', 'innerCV_preds'))
    predy <- levels(ytrain)[(innerpreds > 0) + 1]
    out <- data.frame(testy = ytrain, predy = predy, predyp = innerpreds)
  } else if (fam == "multinomial") {
    # multinomial
    innerpreds <- lapply(x$outer_result, '[[', 'innerCV_preds')
    innerpreds <- do.call(rbind, innerpreds)
    predy <- colnames(innerpreds)[max.col(innerpreds)]
    predy <- factor(predy, levels = levels(ytrain))
    out <- data.frame(testy = ytrain, predy = predy)
    out <- cbind(out, innerpreds)
  } else {
    # regression
    innerpreds <- unlist(lapply(x$outer_result, '[[', 'innerCV_preds'))
    out <- data.frame(testy = ytrain, predy = innerpreds)
  }
  if (!is.null(rownames(x$outer_result[[1]]$innerCV_preds))) {
    rn <- unlist(lapply(1:length(x$outer_result), function(i) {
      paste(names(x$outer_result)[i], rownames(x$outer_result[[i]]$innerCV_preds),
            sep=".")
    }))
    rownames(out) <- rn
  }
  out
}


#' @rdname innercv_preds
#' @export
innercv_preds.nestcv.train <- function(x) {
  names(x$outer_result) <- paste0(names(x$outer_result), ".")
  ytrain <- unlist(lapply(x$outer_result, function(i) i$fit$pred$obs))
  predy <- unlist(lapply(x$outer_result, function(i) i$fit$pred$pred))
  res <- data.frame(testy = ytrain, predy = predy)
  if (x$outer_result[[1]]$fit$modelType == "Classification") {
    if (length(x$outer_result[[1]]$fit$levels) == 2) {
      # binomial
      predyp <- unlist(lapply(x$outer_result, function(i) i$fit$pred[, i$fit$levels[2]]))
    } else {
      # multinomial
      predyp <- lapply(x$outer_result, function(i) i$fit$pred[, i$fit$levels])
      predyp <- do.call(rbind, predyp)
    }
    res <- cbind(res, predyp)
  }
  res
}


#' Summarise performance on inner CV test folds
#' 
#' Calculates performance metrics on inner CV held-out test folds: confusion
#' matrix, accuracy and balanced accuracy for classification; ROC AUC for binary
#' classification; RMSE, R^2 and mean absolute error (MAE) for regression. 
#' 
#' @param x a `nestcv.glmnet` or `nestcv.train` object
#' @return Returns performance metrics from outer training folds, see
#'   [predSummary].
#' @seealso [predSummary]
#' @examples
#' data(iris)
#' x <- iris[, 1:4]
#' y <- iris[, 5]
#' 
#' fit <- nestcv.glmnet(y, x,
#'                      family = "multinomial",
#'                      alpha = 1,
#'                      n_outer_folds = 3)
#' summary(fit)
#' innercv_summary(fit)
#' 
#' @export
innercv_summary <- function(x) {
  innerpreds <- innercv_preds(x)
  predSummary(innerpreds)
}


#' Outer training fold predictions
#' 
#' Obtain predictions on outer training folds which can be used for performance
#' metrics and ROC curves.
#' 
#' @param x a `nestcv.glmnet`, `nestcv.train` or `outercv` fitted object
#' @return Dataframe with columns `ytrain` and `predy` containing observed and
#'   predicted values from training folds. For binomial and multinomial models
#'   additional columns are added with class probabilities or log likelihood
#'   values.
#' @details Note: the argument `outer_train_predict` must be set to `TRUE` in
#'   the original call to either `nestcv.glmnet`, `nestcv.train` or `outercv`.
#' @export
train_preds <- function(x) {
  trainpreds <- lapply(x$outer_result, '[[', 'train_preds')
  trainpreds <- do.call(rbind, trainpreds)
  trainpreds
}


#' Summarise performance on outer training folds
#' 
#' Calculates performance metrics on outer training folds: confusion matrix,
#' accuracy and balanced accuracy for classification; ROC AUC for binary
#' classification; RMSE, R^2 and mean absolute error (MAE) for regression.
#' 
#' @param x a `nestcv.glmnet`, `nestcv.train` or `outercv` object
#' @return Returns performance metrics from outer training folds, see
#'   [predSummary]
#' @details Note: the argument `outer_train_predict` must be set to `TRUE` in
#'   the original call to either `nestcv.glmnet`, `nestcv.train` or `outercv`.
#' @seealso [predSummary]
#' @examples
#' \donttest{
#' data(iris)
#' x <- iris[, 1:4]
#' y <- iris[, 5]
#' 
#' fit <- nestcv.glmnet(y, x,
#'                      family = "multinomial",
#'                      alpha = 1,
#'                      outer_train_predict = TRUE,
#'                      n_outer_folds = 3)
#' summary(fit)
#' innercv_summary(fit)
#' train_summary(fit)
#' 
#' fit2 <- nestcv.train(y, x,
#'                     model="svm",
#'                     outer_train_predict = TRUE,
#'                     n_outer_folds = 3,
#'                     cv.cores = 2)
#' summary(fit2)
#' innercv_summary(fit2)
#' train_summary(fit2)
#' }
#' @export
train_summary <- function(x) {
  trainpreds <- train_preds(x)
  if (is.null(trainpreds)) stop("No saved training prediction data")
  df <- trainpreds
  colnames(df)[colnames(df) == "ytrain"] <- "testy"
  predSummary(df)
}


#' Build ROC curve from outer CV training folds
#' 
#' Build ROC (receiver operating characteristic) curve from outer training
#' folds. Object can be plotted using `plot()` or passed to functions [auc()]
#' etc.
#' 
#' @param x a `nestcv.glmnet`, `nestcv.train` or `outercv` object
#' @param direction Set ROC directionality [pROC::roc]
#' @param ... Other arguments passed to [pROC::roc]
#' @return `"roc"` object, see [pROC::roc]
#' @details Note: the argument `outer_train_predict` must be set to `TRUE` in
#'   the original call to either `nestcv.glmnet`, `nestcv.train` or `outercv`.
#' @export
train_roc <- function(x, direction = "<", ...) {
  trainpreds <- train_preds(x)
  pROC::roc(trainpreds$ytrain, trainpreds$predyp, direction = direction,
            quiet = TRUE, ...)
}
