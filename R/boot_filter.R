
#' Bootstrap for filter functions
#'
#' Randomly samples predictors and averages the ranking to give an ensemble
#' measure of predictor variable importance.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param filterFUN Filter function, e.g. [ttest_filter()].
#' @param B Number of times to bootstrap
#' @param type Type of vector returned. Default "index" returns indices, "full"
#'   returns full output.
#' @param ... Optional arguments passed to the function specified by `filterFUN`
#' @return Integer vector of indices of filtered parameters (`type = "index"`)
#'   or if `type = "full"` a matrix of rankings from each bootstrap is returned.
#' @seealso [boot_ttest()]
#' @export
#' 
boot_filter <- function(y, x, filterFUN, B = 50, type = "index", ...) {
  ranks <- sapply(1:B, function(i) {
    ind <- sample.int(length(y), replace = TRUE)
    out <- filterFUN(y[ind], x[ind, ], type = "full", ...)
    rank(out[, grep("p.?val", colnames(out))])
  })
  meanRank <- rowMeans(ranks)
  names(meanRank) <- colnames(x)
  if (type == "full") return(sort(meanRank))
  order(meanRank)
}


#' Bootstrap t-test filter
#'
#' Randomly samples predictors and averages the ranking from [ttest_filter()] to
#' give an ensemble measure of best predictors by bootstrapped t-test.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param B Number of times to bootstrap
#' @param ... Optional arguments passed to [ttest_filter()]
#' @return Integer vector of indices of filtered parameters (`type = "index"`)
#'   or if `type = "full"` a matrix of rankings from each bootstrap is returned.
#' @seealso [boot_filter()]
#' @export
#' 
boot_ttest <- function(y, x, B = 50, ...) {
  boot_filter(y, x, ttest_filter, B=B, ...)
}

