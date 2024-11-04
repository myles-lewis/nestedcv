
#' Slim nestedcv models
#' 
#' Slims nestedcv objects to only the models in the outer CV folds.
#' 
#' @param x A 'nestedcv' or 'cva.glmnet' fitted model object.
#' @returns For 'nestedcv' objects, a list object of the same class but only
#'   containing `outer_result`. For 'cva.glmnet' models, only the cv.glmnet
#'   model with the best alpha value is kept. Models for all other values of
#'   alpha are discarded.
#' @seealso [nestcv.glmnet()] [nestcv.train()] [outercv()] [cva.glmnet()]
#' @export
slim <- function(x) {
  UseMethod("slim")  
}

#' @export
slim.nestcv.glmnet <- function(x) {
  outer_result <- lapply(x$outer_result, function(i) {
    list(coef = i$coef, cvafit = slim.cva.glmnet(i$cvafit))
  })
  ret <- list(outer_result = outer_result)
  class(ret) <- "nestcv.glmnet"
  ret
}

#' @export
slim.cva.glmnet <- function(x) {
  fits <- x$fits[x$which_alpha]
  alphaSet <- x$alphaSet[x$which_alpha]
  alpha_cvm <- x$alpha_cvm[x$which_alpha]
  ret <- list(fits = fits, alphaSet = alphaSet, alpha_cvm = alpha_cvm,
              best_alpha = x$best_alpha, which_alpha = 1)
  structure(ret, class = "cva.glmnet")
}

#' @export
slim.default <- function(x) {
  outer_result <- lapply(x$outer_result, function(i) list(fit = i$fit))
  ret <- list(outer_result = outer_result)
  class(ret) <- class(x)
  ret
}
