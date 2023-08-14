
#' One-hot encode
#' 
#' Fast one-hot encoding of all factor and character columns in a dataframe to
#' convert it into a numeric matrix by creating dummy (binary) columns.
#' 
#' @param x A dataframe, matrix or tibble. Matrices are returned untouched.
#' @param all_levels Logical, whether to create dummy variables for all levels
#'   of each factor. Default is `FALSE` to avoid issues with regression models.
#' @param rename_binary Logical, whether to rename binary factors by appending
#'   the 2nd level of the factor to aid interpretation of encoded factor levels
#'   and to allow consistency with naming.
#' @param sep Character for separating factor variable names and levels for 
#'   encoded columns.
#' @details
#'  Binary factor columns and logical columns are converted to integers (0 or
#'  1). Multi-level unordered factors are converted to multiple columns of 0/1
#'  (dummy variables): if `all_levels` is set to `FALSE` (the default), then the
#'  first level is assumed to be a reference level and additional columns are
#'  created for each additional level; if `all_levels` is set to `TRUE` one
#'  column is used for each level. Unused levels are dropped. Character columns
#'  are first converted to factors and then encoded. Ordered factors are
#'  replaced by their internal codes. Numeric or integer columns are left
#'  untouched.
#'  
#'  Having dummy variables for all levels of a factor can cause problems with
#'  multicollinearity in regression (the dummy variable trap), so `all_levels`
#'  is set to `FALSE` by default which is necessary for regression models such
#'  as `glmnet` (equivalent to full rank parameterisation). However, setting
#'  `all_levels` to `TRUE` can aid with interpretability (e.g. with SHAP
#'  values), and in some cases filtering might result in some dummy variables
#'  being excluded. Note this function is designed to quickly generate dummy
#'  variables for more general machine learning purposes. To create a proper
#'  design matrix object for regression models, use [model.matrix()].
#' @return A numeric matrix with the same number of rows as the input data.
#'   Dummy variable columns replace the input factor or character columns.
#'   Numeric columns are left intact.
#' @seealso [caret::dummyVars()], [model.matrix()]
#' @examples
#' data(iris)
#' x <- iris
#' x2 <- one_hot(x)
#' head(x2)  # 3 columns for Species
#' 
#' x2 <- one_hot(x, all_levels = FALSE)
#' head(x2)  # 2 columns for Species
#' 
#' @export
#' 
one_hot <- function(x, all_levels = FALSE, rename_binary = TRUE, sep = ".") {
  if (is.matrix(x)) return(x)
  x <- as.data.frame(x)
  x <- droplevels(x)
  factor_ind <- index_factor(x, convert_bin = TRUE)
  bin_ind <- unlist(lapply(x, function(i) {
    !is.numeric(i) && !is.logical(i) && nlevels(factor(i)) == 2
  }))
  if (sum(factor_ind) == 0) {
    out <- fix_binary(x, bin_ind, sep, rename_binary)
    return(out)
  }
  x0 <- fix_binary(x, bin_ind, sep, rename_binary)
  x1 <- x0[, !factor_ind, drop = FALSE]
  f <- which(factor_ind)
  if (all_levels) {
    expand_x2 <- lapply(f, function(i) {
      lev <- levels(factor(x[, i]))
      title <- colnames(x)[i]
      cn <- paste(title, lev, sep = sep)
      m <- matrix(0L, nrow = nrow(x), ncol = length(lev),
                  dimnames = list(NULL, cn))
      m[is.na(x[, i]), ] <- NA
      for (j in seq_along(lev)) {
        m[x[, i] == lev[j], j] <- 1L
      }
      m
    })
  } else {
    expand_x2 <- lapply(f, function(i) {
      lev <- levels(factor(x[, i]))[-1]
      title <- colnames(x)[i]
      cn <- paste(title, lev, sep = sep)
      m <- matrix(0L, nrow = nrow(x), ncol = length(lev),
                  dimnames = list(NULL, cn))
      m[is.na(x[, i]), ] <- NA
      for (j in seq_along(lev)) {
        m[x[, i] == lev[j], j] <- 1L
      }
      m
    })
  }
  expand_x2 <- do.call(cbind, expand_x2)
  out <- cbind(x1, expand_x2)
  out
}

fix_binary <- function(x, bin_ind, sep, rename_binary) {
  if (rename_binary) {
    which_bin <- which(bin_ind)
    lev2 <- unlist(lapply(which_bin, function(i) levels(factor(x[, i]))[2]))
    cn <- paste(colnames(x)[which_bin], lev2, sep = sep)
    colnames(x)[which_bin] <- cn
  }
  x <- data.matrix(x)
  x[, bin_ind] <- x[, bin_ind] - 1
  x
}
