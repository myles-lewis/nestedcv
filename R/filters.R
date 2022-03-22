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


#' Correlation between a vector and a matrix
#'
#' Fast pearson/spearman correlation where y is vector, x is matrix, adapted 
#' from [stats::cor.test].
#' 
#' @param y Numerical vector
#' @param x Matrix
#' @param method Type of correlation, either "pearson" or "spearman".
#' @param use Optional character string giving a method for computing 
#' covariances in the presence of missing values. See [cor]
#' @details For speed, p-values for Spearman's test are computed by 
#' asymptotic t approximation, equivalent to [cor.test] with `exact = FALSE`.
#' @importFrom Rfast Rank colRanks
#' @importFrom stats complete.cases cor pt
#' @export
#' 
correls2 <- function(y, x,
                    method = "pearson",
                    use = "complete.obs") {
  ok <- complete.cases(x, y)
  x <- x[ok, ]
  y <- y[ok]
  n <- length(y)
  df <- n - 2
  if (method == "spearman") {
    r <- as.vector( cor(Rfast::Rank(y), Rfast::colRanks(x)) )
    q <- (n^3 - n) * (1 - r) / 6  # AS 89 algorithm
    den <- (n*(n^2-1))/6
    r <- 1 - q/den
    pval <- 2 * pt(abs(r) / sqrt((1 - r^2)/df), df, lower.tail=FALSE)
    rname <- 'rho'
  } else {
    r <- suppressWarnings(as.vector(cor(y, x, use = use)))
    STATISTIC <- sqrt(df) * r / sqrt(1 - r^2)
    pval <- 2 * pt(abs(STATISTIC), df, lower.tail = FALSE)  # two-tailed
    rname <- 'r'
  }
  res <- cbind(r, pval)
  colnames(res) <- c(rname, 'p-value')
  rownames(res) <- colnames(x)
  res
}


#' Correlation filter
#' 
#' Filter using correlation (Pearson or Spearman) for ranking variables.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a matrix of p-values.
#' @param method Type of correlation, either "pearson" or "spearman".
#' @param ... Further arguments passed to [correls]
#' @return Integer vector of indices of filtered parameters (type = "index") or 
#' character vector of names (type = "names") of filtered parameters. If 
#' `type` is `"full"` full output from [correls] is returned.
#' @export
#' 
correl_filter <- function(y,
                      x,
                      nfilter = NULL,
                      p_cutoff = 0.05,
                      method = "pearson",
                      type = c("index", "names", "full"),
                      ...) {
  type <- match.arg(type)
  res <- correls2(y, x, method = method, ...)
  if (type == "full") return(res)
  out <- res[, "p-value"]
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
#' @seealso [CORElearn::attrEval]
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


#' glmnet filter
#'
#' Filter using properties of elastic net regression using glmnet to calculate
#' variable importance.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return
#' @param method String indicating method of determining variable importance.
#'   "mean" (the default) uses the mean absolute coefficients across the range
#'   of lambdas; "nonzero" counts the number of times variables are retained in
#'   the model across all values of lambda.
#' @param type Type of vector returned. Default "index" returns indices, "names"
#'   returns predictor names, "full" returns full output.
#' @param ... Other arguments passed to [glmnet]
#' @details The glmnet elastic net mixing parameter alpha can be varied to
#'   include a larger number of predictors. Default alpha = 1 is pure LASSO,
#'   resulting in greatest sparsity, while alpha = 0 is pure ridge regression,
#'   retaining all predictors in the regression model. Note, the `family`
#'   argument is commonly needed, see [glmnet].
#' @seealso [glmnet]
#' @importFrom glmnet glmnet
#' @export
#' 
glmnet_filter <- function(y,
                          x,
                          nfilter = NULL,
                          method = c("mean", "nonzero"),
                          type = c("index", "names", "full"),
                          ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  fit <- glmnet(x, y, ...)
  cf <- as.matrix(coef(fit))
  if (method == "mean") {
    cf <- abs(cf)
    out <- rowMeans(cf)  # mean abs coefs
  } else {
    cf <- cf != 0
    out <- rowSums(cf)  # number of nonzero coefs
  }
  out <- out[-1]  # remove intercept
  if (type == "full") return(out)
  if (type == "index") names(out) <- 1:ncol(x)
  out <- out[out != 0]
  out <- sort(out, decreasing = TRUE)
  if (!is.null(nfilter)) out <- out[1:min(nfilter, length(out))]
  out <- names(out)
  if (length(out) == 0) stop("No predictors selected")
  if (type == "index") out <- as.integer(out)
  out
}

# Code modified from caret::findCorrelation
# https://github.com/topepo/caret/blob/master/pkg/caret/R/findCorrelation.R

#' Filter to reduce collinearity in predictors
#'
#' This function identifies predictors with r^2 above a given cut-off and
#' produces an index of predictors to be removed. The function takes a matrix or
#' data.frame of predictors, and the columns need to be ordered in terms of
#' importance - first column of any pair that are correlated is retained and
#' subsequent columns which correlate above the cut-off are flagged for removal.
#'
#' @param x A matrix or data.frame of values. The order of columns is used to
#'   determine which columns to retain, so the columns in `x` should be sorted
#'   with the most important columns first.
#' @param rsq_cutoff Value of cut-off for r-squared
#' @param verbose Boolean whether to print details
#' @return Integer vector of the indices of columns in `x` to remove due to
#'   collinearity
#' @export
#' 
collinear <- function(x, rsq_cutoff = 0.9, verbose = FALSE) {
  rsq <- cor(x)^2
  rsq[lower.tri(rsq, diag = TRUE)] <- NA
  combsAboveCutoff <- which(rsq > rsq_cutoff)
  colsToCheck <- ceiling(combsAboveCutoff / nrow(rsq))
  rowsToCheck <- combsAboveCutoff %% nrow(rsq)
  colsToDiscard <- colsToCheck > rowsToCheck
  rowsToDiscard <- !colsToDiscard
  if (any(rowsToDiscard)) warning("Unexpected rows to discard")
  if (verbose) {
    df <- data.frame(keep = rowsToCheck, remove = colsToCheck)
    df <- df[order(df$keep), ]
    remd <- NULL
    for (i in unique(df$keep)) {
      rem <- df$remove[df$keep %in% i]
      rem <- rem[!rem %in% remd]
      if (length(rem) > 0) {
        if (!i %in% df$remove) cat("Keep ") else cat ("Removed ")
        cat(paste0(colnames(x)[i], ", remove "))
        cat(paste(colnames(x)[rem], collapse = ", "))
        remd <- c(remd, rem)
        cat("\n")
      }
    }
  }
  deletecol <- colsToCheck
  deletecol <- unique(deletecol)
  deletecol
}
