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
#' @export innercv_roc
#' 
innercv_roc <- function(x, ...) {
  UseMethod("innercv_roc")
}


#' @rdname innercv_roc
#' @importFrom pROC roc auc
#' @export
#' 
innercv_roc.nestcv.glmnet <- function(x, direction = "<", ...) {
  innerpreds <- unlist(lapply(x$outer_result, '[[', 'innerCV_preds'))
  ytrain <- unlist(lapply(x$outer_result, '[[', 'ytrain'))
  pROC::roc(ytrain, innerpreds, direction = direction, ...)
}


#' @rdname innercv_roc
#' @importFrom pROC roc auc
#' @export
#' 
innercv_roc.nestcv.train <- function(x, direction = "<", ...) {
  innerpreds <- unlist(lapply(x$outer_result, function(i) i$fit$pred[, i$fit$levels[2]]))
  ytrain <- unlist(lapply(x$outer_result, function(i) i$fit$pred$obs))
  pROC::roc(ytrain, innerpreds, direction = direction, ...)
}
