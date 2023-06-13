
#' Univariate filter for binary classification
#'
#' Univariate statistic filter for dataframes of predictors with mixed numeric
#' and categorical datatypes. Different statistical tests are used depending on
#' the data type of predictors: t-test for continuous data, chi-squared test for
#' categorical data.
#'
#' @param y Response vector
#' @param x Matrix or dataframe of predictors
#' @param force_vars Vector of column names within `x` which are always retained
#'   in the model (i.e. not filtered). Default `NULL` means all predictors will
#'   be passed to `filterFUN`.
#' @param nfilter Number of predictors to return. If `NULL` all predictors with
#'   p-values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param rsq_cutoff r^2 cutoff for removing predictors due to collinearity.
#'   Default `NULL` means no collinearity filtering. Predictors are ranked based
#'   on t-test. If 2 or more predictors are collinear, the first ranked
#'   predictor by t-test is retained, while the other collinear predictors are
#'   removed. See [collinear()].
#' @param type Type of vector returned. Default "index" returns indices, "names"
#'   returns predictor names, "full" returns a list of p values.
#' @param ... optional arguments, e.g. `rsq_method`: see [collinear()].
#'
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters in order
#'   of t-test/chi-squared test p-value. If `type` is `"full"` full output is
#'   returned containing a matrix of results from [Rfast::ttests] and a matrix
#'   of [chisq.test()] results.
#' @importFrom stats chisq.test
#' @export
#'
bin_stat_filter <- function(y,
                            x,
                            force_vars = NULL,
                            nfilter = NULL,
                            p_cutoff = 0.05,
                            rsq_cutoff = NULL,
                            type = c("index", "names", "full"),
                            ...) {
  type <- match.arg(type)
  if (length(unique(y)) != 2) stop("y is not binary")
  factor_ind <- index_factor(x)
  if (sum(factor_ind) == 0) {
    return(ttest_filter(y, x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, ...))
  }
  if (is.null(colnames(x))) colnames(x) <- seq_len(ncol(x))
  x1 <- as.matrix(x[, !factor_ind])
  x2 <- x[, factor_ind]
  yi1 <- as.numeric(y) == 1
  yi2 <- as.numeric(y) == 2
  res1 <- NULL
  if (ncol(x1) > 0) {
    res1 <- Rfast::ttests(x1[yi1, ], x1[yi2, ])
    rownames(res1) <- colnames(x1)
  }
  res2 <- chisq.tests(x2, y)
  if (type == "full") return(list("t.test" = res1, "chisq.test" = res2))
  rescomb <- matrix(nrow = ncol(x), ncol = 3,
                    dimnames = list(colnames(x), colnames(res2)))
  rescomb[!factor_ind, ] <- res1
  rescomb[factor_ind, ] <- res2
  filter_end(rescomb[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff,
             type, keep_factors = FALSE, factor_ind = NULL, ...)
}


chisq.tests <- function(x, y) {
  res <- t(apply(x, 2, function(j) {
    chsq <- suppressWarnings(chisq.test(j, y))
    unlist(chsq[c(1, 3, 2)])
  }))
  colnames(res) <- c("stat", "pvalue", "dof")
  res
}


index_factor <- function(x) {
  if (is.matrix(x)) return(NULL)
  !unlist(lapply(x, is.numeric))
}

