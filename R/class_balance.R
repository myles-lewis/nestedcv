
#' Check class balance in training folds
#' 
#' @param object Object of class `nestedcv.glmnet`, `nestcv.train` or `outercv`
#' @return Invisibly a table of the response classes in the training folds
#' @export
#' 
class_balance <- function(object) {
  UseMethod("class_balance")
}


#' @rdname class_balance
#' @export
#' 
class_balance.nestcv.glmnet <- function(object) {
  ytrain <- unlist(lapply(object$outer_result, '[[', 'ytrain'))
  if (is.numeric(ytrain)) stop("Not classification", call. = FALSE)
  tab <- table(ytrain)
  cat("Training folds:\n")
  print(c(tab))
  invisible(tab)
}


#' @rdname class_balance
#' @export
#' 
class_balance.nestcv.train <- function(object) {
  ytrain <- unlist(lapply(object$outer_result, function(i) i$fit$pred$obs))
  if (is.numeric(ytrain)) stop("Not classification", call. = FALSE)
  tab <- table(ytrain)
  cat("Training folds:\n")
  print(c(tab))
  invisible(tab)
}


#' @rdname class_balance
#' @export
#' 
class_balance.outercv <- function(object) {
  class_balance.nestcv.glmnet(object)
}
