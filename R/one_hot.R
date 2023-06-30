
#' One-hot encode
#' 
#' One-hot encoding of factors in dataframes. Binary factors are each converted
#' to a single column of 0 and 1. Multi-level unordered factors are converted to
#' multiple columns, one for each level. Unused levels are dropped. Numeric or
#' integer columns are left untouched.
#' 
#' @param x A dataframe or matrix. Matrices are returned untouched.
#' @param sep Character for separating factor variable names and levels for 
#'   encoded columns.
#' @param rename_binary Logical, whether to rename binary factors by appending
#'   the 2nd level of the factor to aid interpretation of encoded factor levels
#'   and to allow consistency with naming.
#' @return A numeric matrix with multi-level factors converted to one-hot
#'   encoded extra columns with 0 and 1. Binary factors are each converted to
#'   single columns of 0 and 1. Ordered factors are converted to integer levels.
#' @seealso [caret::dummyVars()]
#' @examples
#' data(iris)
#' x <- iris
#' x2 <- one_hot(x)
#' head(x2)
#' 
#' @export
#' 
one_hot <- function(x, sep = ".", rename_binary = TRUE) {
  if (is.matrix(x)) return(x)
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
  expand_x2 <- lapply(f, function(i) {
    lev <- levels(factor(x[, i]))
    title <- colnames(x)[i]
    cn <- paste(title, lev, sep = sep)
    m <- matrix(0, nrow = nrow(x), ncol = length(lev))
    m[is.na(x[, i]), ] <- NA
    for (j in seq_along(lev)) {
      m[x[, i] == lev[j], j] <- 1
    }
    colnames(m) <- cn
    m
  })
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
