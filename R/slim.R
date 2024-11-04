
#' Slim nestedcv models
#' 
#' Slims nestedcv objects to only the models in the outer CV folds.
#' 
#' @param x A 'nestedcv' fitted model object.
#' @returns A list object of the same class but only containing `outer_result`.
#' @seealso [nestcv.glmnet()] [nestcv.train()] [outercv()] 
#' @export
slim <- function(x) {
  UseMethod("slim")  
}

#' @export
slim.nestcv.glmnet <- function(x) {
  outer_result <- lapply(x$outer_result, function(i) {
    list(coef = i$coef, cvafit = i$cvafit)
  })
  ret <- list(outer_result = outer_result)
  class(ret) <- "nestcv.glmnet"
  ret
}

#' @export
slim.default <- function(x) {
  outer_result <- lapply(x$outer_result, function(i) list(fit = i$fit))
  ret <- list(outer_result = outer_result)
  class(ret) <- class(x)
  ret
}
