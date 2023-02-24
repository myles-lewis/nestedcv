
var_dir <- function(y, x) {
  if (is.numeric(y)) {
    # regression
    cc <- correl_filter(y, x, p_cutoff = NULL, type = "full")
    return(sign(cc[, "r"]))
  }
  # classification
  if (nlevels(y) != 2) return(NULL)
  tt <- ttest_filter(y, x, p_cutoff = NULL, type = "full")
  -sign(tt[, "stat"])
}


#' Variable directionality
#'
#' Determines directionality of final predictors for binary or regression
#' models.
#'
#' @param object a `nestcv.glmnet` or `nestcv.train` fitted model
#' @return named vector showing sign of final predictors. If the response vector
#'   is multinomial `NULL` is returned.
#' @export

var_direction <- function(object) {
  y <- object$y
  x <- object$xsub
  var_dir(y, x)
}

