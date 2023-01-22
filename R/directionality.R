
var_dir <- function(y, x) {
  tt <- ttest_filter(y, x, p_cutoff = NULL, type = "full")
  sign(tt[, "stat"])
}


#' Variable directionality
#'
#' Determines directionality of final predictors for binary models.
#' 
#' @param object a `nestcv.glmnet` or `nestcv.train` fitted model
#' @return named vector showing sign of final predictors
#' @export
var_direction <- function(object) {
  UseMethod("var_direction")
}

#' @rdname var_direction
#' @export
var_direction.default <- function(object) {
  y <- object$final_data$y
  x <- object$final_data$x
  var_dir(y, x)
}

#' @rdname var_direction
#' @export
var_direction.nestcv.train <- function(object) {
  y <- object$final_fit$trainingData$.outcome
  x <- object$final_fit$trainingData
  x$.outcome <- NULL
  x <- as.matrix(x)
  var_dir(y, x)
}
