# Extract inner CV ROC curve


#' Build ROC curve from left-out folds from inner CV
#' 
#' Build ROC (receiver operating characteristic) curve from left-out folds 
#' from inner CV. Object can be plotted using `plot()` or passed to functions 
#' [auc()] etc.
#' 
#' @param x Fitted `nestedcv` object 
#' @param direction Set ROC directionality [pROC::roc]
#' @param ... Other arguments passed to [pROC::roc]
#' @return `"roc"` object, see [pROC::roc]
#' @examples
#' 
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
#'                       filter_options = list(nfilter = 100))
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
#'
#' @export innercv_roc
#' 
innercv_roc <- function(x, direction = "<", ...) {
  innerpreds <- innercv_preds(x)
  pROC::roc(innerpreds$testy, innerpreds$predyp, direction = direction, quiet = TRUE, ...)
}


#' Inner CV predictions
#' 
#' Obtain predictions on held-out test inner CV folds
#' 
#' @param x Fitted `nestedcv` object
#' @return Dataframe with 2 columns `obs` and `pred`
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
  if (x$call$family == "binomial") {
    # binomial
    innerpreds <- unlist(lapply(x$outer_result, '[[', 'innerCV_preds'))
    predy <- levels(ytrain)[(innerpreds > 0) + 1]
    out <- data.frame(testy = ytrain, predy = predy, predyp = innerpreds)
  } else if (x$call$family == "multinomial") {
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
  innerpreds <- unlist(lapply(x$outer_result, function(i) i$fit$pred[, i$fit$levels[2]]))
  ytrain <- unlist(lapply(x$outer_result, function(i) i$fit$pred$obs))
  data.frame(testy = ytrain, predyp = innerpreds)
}
