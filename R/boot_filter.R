
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
boot_filter <- function(y, x, filterFUN, B = 50,
                        nfilter = NULL, type = "index", ...) {
  ranks <- sapply(1:B, function(i) {
    ind <- sample.int(length(y), replace = TRUE)
    out <- filterFUN(y[ind], x[ind, ], type = "full", ...)
    rank(out[, grep("p.?val", colnames(out))])
  })
  meanRank <- rowMeans(ranks)
  names(meanRank) <- colnames(x)
  if (type == "full") return(sort(meanRank))
  out <- order(meanRank)
  if (!is.null(nfilter)) out <- out[1:nfilter]
  out
}


#' Bootstrap univariate filters
#'
#' Randomly samples predictors and averages the ranking from filtering functions
#' including [ttest_filter()], [wilcoxon_filter()], [anova_filter()],
#' [correl_filter()] and [lm_filter()] to give an ensemble measure of best
#' predictors by repeated random sampling subjected to a statistical test.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param B Number of times to bootstrap
#' @param ... Optional arguments passed to the filter function
#' @return Integer vector of indices of filtered parameters (`type = "index"`),
#'   or if `type = "full"`, a matrix of rankings from each bootstrap is
#'   returned.
#' @seealso [ttest_filter()], [wilcoxon_filter()], [anova_filter()],
#'   [correl_filter()], [lm_filter()] and [boot_filter()]
#' @export
#' 
boot_ttest <- function(y, x, B = 50, ...) {
  boot_filter(y, x, ttest_filter, B=B, ...)
}


#' @rdname boot_ttest
#' @export
#' 
boot_wilcoxon <- function(y, x, B = 50, ...) {
  boot_filter(y, x, wilcoxon_filter, B=B, ...)
}


#' @rdname boot_ttest
#' @export
#' 
boot_anova <- function(y, x, B = 50, ...) {
  boot_filter(y, x, anova_filter, B=B, ...)
}


#' @rdname boot_ttest
#' @export
#' 
boot_correl <- function(y, x, B = 50, ...) {
  boot_filter(y, x, correl_filter, B=B, ...)
}


#' @rdname boot_ttest
#' @export
#' 
boot_lm <- function(y, x, B = 50, ...) {
  boot_filter(y, x, lm_filter, B=B, ...)
}
