# Filters to reduce number of predictors in nested CV

#' t-test filter
#' 
#' Simple univariate filter using t-test using the Rfast package for speed.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a matrix of p values.
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` full output from [Rfast::ttests] is returned.
#' @importFrom Rfast ttests
#' @export
#' 
ttest_filter <- function(y,
                       x,
                       nfilter = NULL,
                       p_cutoff = 0.05,
                       type = c("index", "names", "full")) {
  type <- match.arg(type)
  y <- factor(y)
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  res <- Rfast::ttests(x[indx1, ], x[indx2, ])
  rownames(res) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(res)
  out <- res[, "pvalue"]
  out <- sort(out)
  if (!is.null(nfilter)) out <- out[1:min(nfilter, length(out))]
  if (!is.null(p_cutoff)) out <- out[out < p_cutoff]
  out <- names(out)
  if (length(out) == 0) stop("No predictors selected")
  if (type == "index") out <- as.integer(out)
  out
}


#' ANOVA filter
#' 
#' Simple univariate filter using anova (Welch's F-test) using the 
#' Rfast package for speed.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a matrix of p values.
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` full output from [Rfast::ftests] is returned.
#' @importFrom Rfast ftests
#' @export
#' 
anova_filter <- function(y,
                       x,
                       nfilter = NULL,
                       p_cutoff = 0.05,
                       type = c("index", "names", "full")) {
  type <- match.arg(type)
  y <- factor(y)
  res <- Rfast::ftests(x, y)
  rownames(res) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(res)
  out <- res[, "pval"]
  out <- sort(out)
  if (!is.null(nfilter)) out <- out[1:min(nfilter, length(out))]
  if (!is.null(p_cutoff)) out <- out[out < p_cutoff]
  out <- names(out)
  if (length(out) == 0) stop("No predictors selected")
  if (type == "index") out <- as.integer(out)
  out
}


#' Wilcoxon test filter
#' 
#' Simple univariate filter using Wilcoxon (Mann-Whitney) test using the 
#' matrixTests package.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a matrix of p-values.
#' @param exact Logical whether exact or approximate p-value is calculated. 
#' Default is `FALSE` for speed.
#' @param ... Further arguments passed to [matrixTests::row_wilcoxon_twosample]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` full output from [matrixTests::row_wilcoxon_twosample] 
#' is returned.
#' @importFrom matrixTests row_wilcoxon_twosample
#' @export
#' 
wilcoxon_filter <- function(y,
                         x,
                         nfilter = NULL,
                         p_cutoff = 0.05,
                         type = c("index", "names", "full"),
                         exact = FALSE,
                         ...) {
  type <- match.arg(type)
  y <- factor(y)
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  res <- suppressWarnings(
    matrixTests::row_wilcoxon_twosample(t(x[indx1, ]), t(x[indx2, ]),
                                        exact = FALSE, ...)
  )
  if (type == "full") return(res)
  out <- res[, "pvalue"]
  names(out) <- if (type == "index") 1:ncol(x) else colnames(x)
  out <- sort(out)
  if (!is.null(nfilter)) out <- out[1:min(nfilter, length(out))]
  if (!is.null(p_cutoff)) out <- out[out < p_cutoff]
  out <- names(out)
  if (length(out) == 0) stop("No predictors selected")
  if (type == "index") out <- as.integer(out)
  out
}


#' Random forest filter
#' 
#' Fits a random forest model and ranks variables by variable importance. 
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ntree Number of trees to grow. See [randomForest].
#' @param mtry Number of predictors randomly sampled as candidates at each 
#' split. See [randomForest].
#' @param ... Optional arguments passed to [randomForest].
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` a named vector of variable importance is returned.
#' @details
#' This filter uses the [randomForest] function from the randomForest package.
#' Variable importance is calculated using the [importance] function, specifying
#' type 1 = mean decrease in accuracy. See [importance].
#' @importFrom randomForest randomForest importance
#' @export
#' 
rf_filter <- function(y, x, nfilter = NULL,
                      type = c("index", "names", "full"),
                      ntree = 1000,
                      mtry = ncol(x) * 0.2,
                      ...) {
  type <- match.arg(type)
  fit <- randomForest::randomForest(x, y, importance = TRUE,
                                    ntree = ntree, mtry = mtry, ...)
  vi <- as.vector(importance(fit, type = 1))
  names(vi) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) vi <- vi[1:min(nfilter, length(vi))]
  out <- names(vi)
  if (type == "index") out <- as.integer(out)
  out
}

#' ReliefF filter
#' 
#' Uses ReliefF algorithm from the CORElearn package to rank predictors in order 
#' of importance.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param estimator Type of algorithm used, see [CORElearn::attrEval]
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ... Other arguments passed to [CORElearn::attrEval]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` a named vector of variable importance is returned.
#' @importFrom CORElearn attrEval
#' @export
#'
relieff_filter <- function(y, x, nfilter = NULL, 
                           estimator = "ReliefFequalK",
                           type = c("index", "names", "full"), ...) {
  type <- match.arg(type)
  df <- as.data.frame(x)
  df$y <- y
  ref <- CORElearn::attrEval('y', df, estimator = estimator, ...)
  names(ref) <- if (type == "index") 1:ncol(x) else colnames(x)
  if (type == "full") return(ref)
  ref <- sort(ref, decreasing = TRUE)
  if (!is.null(nfilter)) ref <- ref[1:min(nfilter, length(ref))]
  out <- names(ref)
  if (type == "index") out <- as.integer(out)
  out
}

#' Combo filter
#' 
#' Filter combining univariate (t-test or anova) filtering and reliefF filtering 
#' in equal measure.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return, using 1/2 from `ttest_filter` 
#' or `anova_filter` and 1/2 from `relieff_filter`. Since `unique` is applied, 
#' the final number returned may be less than `nfilter`.
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns full output.
#' @param ... Optional arguments passed via [relieff_filter] to 
#' [CORElearn::attrEval]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If `type` 
#' is `"full"` a list containing full outputs from either [ttest_filter] or 
#' [anova_filter] and [relieff_filter] is returned.
#' @export
#' 
combo_filter <- function(y, x,
                         nfilter,
                         type = c("index", "names", "full"), ...) {
  uni_set <- if (nlevels(y) == 2) {ttest_filter(y, x, nfilter, type = type)
  } else anova_filter(y, x, nfilter, type = type)
  relf_set <- relieff_filter(y, x, nfilter, type = type, ...)
  if (type == "full") {
    return(list(unifilt = uni_set, relieff_filter = relf_set))
  }
  n <- round(nfilter / 2)
  unique(c(uni_set[1:min(n, length(uni_set))],
           relf_set[1:min(n, length(relf_set))]))
}
