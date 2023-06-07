
var_dir <- function(y, x) {
  if (is.numeric(y)) {
    # regression
    cc <- correl_filter(y, x, p_cutoff = NULL, type = "full")
    return(sign(cc[, "r"]))
  }
  # classification
  if (nlevels(y) != 2) return(NULL)
  x <- data.matrix(x)
  tt <- ttest_filter(y, x, p_cutoff = NULL, type = "full")
  -sign(tt[, "stat"])
}


#' Variable directionality
#'
#' Determines directionality of final predictors for binary or regression
#' models, using the sign of the t-statistic or correlation coefficient
#' respectively for each variable compared to the outcomes.
#'
#' @param object a `nestcv.glmnet` or `nestcv.train` fitted model
#' @return named vector showing the directionality of final predictors. If the
#'   response vector is multinomial `NULL` is returned.
#' @details 
#' Categorical features with >2 levels are assumed to have a meaningful order
#' for the purposes of directionality. Factors are coerced to ordinal using
#' `data.matrix()`. If factors are multiclass then directionality results should
#' be ignored.
#' @export

var_direction <- function(object) {
  y <- object$y
  x <- object$xsub
  var_dir(y, x)
}

