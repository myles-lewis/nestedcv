
#' Univariate filter for binary classification with mixed predictor datatypes
#'
#' Univariate statistic filter for dataframes of predictors with mixed numeric
#' and categorical datatypes. Different statistical tests are used depending on
#' the data type of response vector and predictors:
#' \describe{
#'   \item{Binary class response: `bin_stat_filter()`}{t-test for continuous 
#'   data, chi-squared test for categorical data}
#'   \item{Multiclass response: `class_stat_filter()`}{one-way ANOVA for 
#'   continuous data, chi-squared test for categorical data}
#'   \item{Continuous response: `lm_stat_filter()`}{linear regression for 
#'   continuous data and binary data, one-way ANOVA for categorical data}
#' }
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
#'   returns predictor names, "full" returns a dataframe of statistics, "list"
#'   returns a list of 2 matrices of statistics, one for continuous predictors,
#'   one for categorical predictors.
#' @param ... optional arguments, e.g. `rsq_method`: see [collinear()].
#'
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters in order
#'   of test p-value. If `type` is `"full"` full output is
#'   returned containing a dataframe of statistical results. If `type` is
#'   `"list"` the output is returned as a list of 2 matrices containing
#'   statistical results separated by continuous and categorical predictors.
#' @details
#' `stat_filter()` is a wrapper which calls `bin_stat_filter()`,
#' `class_stat_filter()` or `lm_stat_filter()` depending on whether `y` is
#' binary, multiclass or continuous respectively.
#' 
#' @export
#' 
stat_filter <- function(y, x, ...) {
  if (is.factor(y) || is.character(y)) {
    if (nlevels(factor(y)) == 2) {
      bin_stat_filter(y, x, ...)
    } else class_stat_filter(y, x, ...)
  } else {
    lm_stat_filter(y, x, ...)
  }
}


#' @rdname stat_filter
#' @importFrom stats chisq.test
#' @export
#'
bin_stat_filter <- function(y,
                            x,
                            force_vars = NULL,
                            nfilter = NULL,
                            p_cutoff = 0.05,
                            rsq_cutoff = NULL,
                            type = c("index", "names", "full", "list"),
                            ...) {
  type <- match.arg(type)
  if (length(unique(y)) != 2) stop("y is not binary")
  factor_ind <- index_factor(x)
  if (sum(factor_ind) == 0) {
    return(ttest_filter(y, x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, ...))
  }
  if (is.null(colnames(x))) colnames(x) <- seq_len(ncol(x))
  x1 <- as.matrix(x[, !factor_ind])
  x2 <- x[, factor_ind, drop = FALSE]
  yi1 <- as.numeric(y) == 1
  yi2 <- as.numeric(y) == 2
  res1 <- NULL
  if (ncol(x1) > 0) {
    res1 <- Rfast::ttests(x1[yi1, ], x1[yi2, ])
    rownames(res1) <- colnames(x1)
  }
  res2 <- chisq.tests(x2, y)
  if (type == "list") return(list("t.test" = res1, "chisq.test" = res2))
  rescomb <- matrix(nrow = ncol(x), ncol = 3,
                    dimnames = list(colnames(x), colnames(res2)))
  rescomb[!factor_ind, ] <- res1
  rescomb[factor_ind, ] <- res2
  if (type == "full") {
    test <- ifelse(factor_ind, "chi-squared", "t-test")
    test <- data.frame(test)
    return(cbind(test, rescomb))
  }
  filter_end(rescomb[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff,
             type, keep_factors = FALSE, factor_ind = NULL, ...)
}

#' @rdname stat_filter
#' @export
#'
class_stat_filter <- function(y,
                              x,
                              force_vars = NULL,
                              nfilter = NULL,
                              p_cutoff = 0.05,
                              rsq_cutoff = NULL,
                              type = c("index", "names", "full", "list"),
                              ...) {
  type <- match.arg(type)
  factor_ind <- index_factor(x)
  if (sum(factor_ind) == 0) {
    return(anova_filter(y, x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, ...))
  }
  if (is.null(colnames(x))) colnames(x) <- seq_len(ncol(x))
  x1 <- as.matrix(x[, !factor_ind])
  x2 <- x[, factor_ind, drop = FALSE]
  y <- factor(y)
  res1 <- NULL
  if (ncol(x1) > 0) {
      res1 <- matrixTests::col_oneway_welch(x1, y)
      rownames(res1) <- colnames(x1)
  }
  res2 <- chisq.tests(x2, y)[, 1:2]
  if (type == "list") return(list("F-test" = res1, "chisq.test" = res2))
  rescomb <- matrix(nrow = ncol(x), ncol = 2,
                    dimnames = list(colnames(x), colnames(res2)))
  rescomb[!factor_ind, ] <- res1
  rescomb[factor_ind, ] <- res2
  if (type == "full") {
    test <- ifelse(factor_ind, "chi-squared", "F-test")
    test <- data.frame(test)
    return(cbind(test, rescomb))
  }
  filter_end(rescomb[, "pvalue"],
             x, force_vars, nfilter, p_cutoff, rsq_cutoff,
             type, keep_factors = FALSE, factor_ind = NULL, ...)
}


#' @rdname stat_filter
#' @export
#'
lm_stat_filter <- function(y,
                           x,
                           force_vars = NULL,
                           nfilter = NULL,
                           p_cutoff = 0.05,
                           rsq_cutoff = NULL,
                           type = c("index", "names", "full", "list"),
                           ...) {
  type <- match.arg(type)
  factor_ind <- index_factor(x, convert_bin = TRUE)
  if (sum(factor_ind) == 0)
    return(lm_filter(y, x, force_vars, nfilter, p_cutoff, rsq_cutoff, type, ...))
  if (is.null(colnames(x)))
    colnames(x) <- seq_len(ncol(x))
  x1 <- data.matrix(x[, !factor_ind])
  x2 <- x[, factor_ind, drop = FALSE]
  res1 <- NULL
  if (ncol(x1) > 0) {
    res1 <- lm_filter(y, x1, type = "full")
  }
  res2 <- oneway.tests.g(y, x2)
  if (type == "list") return(list("lm" = res1, "oneway.test" = res2))
  rescomb <- matrix(nrow = ncol(x), ncol = 2,
                    dimnames = list(colnames(x), c("stat", "pvalue")))
  rescomb[!factor_ind, ] <- res1[, 3:4]
  rescomb[factor_ind, ] <- res2[, c(1, 4)]
  if (type == "full") {
    test <- ifelse(factor_ind, "F-test", "lm.t-test")
    test <- data.frame(test)
    return(cbind(test, rescomb))
  }
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


index_factor <- function(x, convert_bin = FALSE) {
  if (is.matrix(x)) return(NULL)
  if (convert_bin) {
    num_ind <- unlist(lapply(x, function(i) {
      is.numeric(i) || nlevels(factor(i)) <= 2
    }))
  } else {
    num_ind <- unlist(lapply(x, is.numeric))
  }
  !num_ind
}


which_factor <- function(x) {
  if (is.matrix(x)) return(NULL)
  which(index_factor(x, convert_bin = TRUE))
}


# apply over g
#' @importFrom stats oneway.test
oneway.tests.g <- function(x, g) {
  res <- t(apply(g, 2, function(gi) {
    fit <- suppressWarnings(oneway.test(x ~ gi))
    unlist(fit[1:3])
  }))
  colnames(res) <- c("stat", "num df", "denom df", "pvalue")
  res
}
