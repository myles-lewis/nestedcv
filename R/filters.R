# Filters to reduce number of predictors in nested CV

#' Univariate filter
#' 
#' Simple univariate filter using t-test or anova (Welch's F-test) using the 
#' Rfast package for speed.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors with 
#' p values < `p_cutoff` are returned.
#' @param p_cutoff p value cut-off
#' @param return Type of vector returned. Default "names" returns predictor 
#' names, "full" returns a named vector of p values.
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a named vector of p values is returned.
#' @importFrom Rfast ttests ftests
#' @export
#' 
uni_filter <- function(y,
                       x,
                       nfilter = NULL,
                       p_cutoff = 0.05,
                       return = "names") {
  y <- factor(y)
  if (nlevels(y) == 2) {
    indx1 <- as.numeric(y) == 1
    indx2 <- as.numeric(y) == 2
    res <- Rfast::ttests(x[indx1, ], x[indx2, ])
  } else {
    res <- Rfast::ftests(x, y)
  }
  rownames(res) <- colnames(x)
  if (return == "full") return(res)
  out <- res[, grep("pval", colnames(res))]
  out <- sort(out[out < p_cutoff])
  if (!is.null(nfilter)) out <- out[1:nfilter]
  names(out)
}

#' Random forest filter
#' 
#' Fits a random forest model and ranks variables by variable importance. 
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return. If `NULL` all predictors are 
#' returned.
#' @param return Type of vector returned. Default "names" returns predictor 
#' names, "full" returns a named vector of variable importance.
#' @param ntree Number of trees to grow. See [randomForest].
#' @param mtry Number of predictors randomly sampled as candidates at each 
#' split. See [randomForest].
#' @param ... Optional arguments passed to [randomForest].
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a named vector of variable importance is returned.
#' @details
#' This filter uses the [randomForest] function from the randomForest package.
#' Variable importance is calculated using the [importance] function, specifying
#' type 1 = mean decrease in accuracy. See [importance].
#' @importFrom randomForest randomForest importance
#' @export
#' 
rf_filter <- function(y, x, nfilter = NULL, return = "names",
                      ntree = 1000,
                      mtry = ncol(x) * 0.2,
                      ...) {
  fit <- randomForest::randomForest(x, y, importance = TRUE,
                                    ntree = ntree, mtry = mtry, ...)
  vi <- as.vector(importance(fit, type = 1))
  names(vi) <- colnames(x)
  if (return == "full") return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) vi <- vi[1:min(nfilter, length(vi))]
  names(vi)
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
#' @param return Type of vector returned. Default "names" returns predictor 
#' names, "full" returns a named vector of variable importance.
#' @param ... Other arguments passed to [CORElearn::attrEval]
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a named vector of variable importance is returned.
#' @importFrom CORElearn attrEval
#' @export
#'
relieff_filter <- function(y, x, nfilter = NULL, 
                           estimator = "ReliefFequalK",
                           return = "names", ...) {
  df <- as.data.frame(x)
  df$y <- y
  ref <- CORElearn::attrEval('y', df, estimator = estimator, ...)
  names(ref) <- colnames(x)
  if (return == "full") return(ref)
  ref <- sort(ref, decreasing = TRUE)
  if (!is.null(nfilter)) ref <- ref[1:min(nfilter, length(ref))]
  names(ref)
}

#' Combo filter
#' 
#' Filter combining univariate filtering and reliefF filtering in equal measure.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Number of predictors to return, using 1/2 from `uni_filter` 
#' and 1/2 from `relieff_filter`. Since `unique` is applied, the final number 
#' returned may be less than `nfilter`.
#' @param return Type of output returned. Default "names" returns predictor 
#' names, "full" returns full output.
#' @param ... Optional arguments passed via [relieff_filter] to 
#' [CORElearn::attrEval]
#' @return Ordered vector of names of filtered parameters. If `return` is 
#' `"full"` a list of 2 vectors of full output from both [uni_filter] 
#' and [relieff_filter] is returned.
#' @export
#' 
combo_filter <- function(y, x, nfilter, return = "names", ...) {
  uni_set <- uni_filter(y, x, nfilter, return = return)
  relf_set <- relieff_filter(y, x, nfilter, return = return, ...)
  if (return == "full") {
    return(list(unifilt = uni_set, relieff_filter = relf_set))
  }
  n <- round(nfilter / 2)
  unique(c(uni_set[1:n], relf_set[1:n]))
}
