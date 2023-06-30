
#' One-hot encode
#' 
#' One-hot encoding of factors in dataframes. Binary factors are converted to a
#' single column of 0 and 1. Multi-level unordered factors are converted to
#' multiple columns, one for each level. Unused levels are dropped. Numeric or
#' integer columns are left untouched.
#' 
#' @param x A dataframe or matrix. Matrices are returned untouched.
#' @param sep Character for separating factor variable names and levels when new
#'   columns are generated.
#' @return A numeric matrix with multi-level factors converted to one-hot
#'   encoded extra columns with 0 and 1. Binary factors are converted to single
#'   columns of 0 and 1. Ordered factors are converted to integer levels.
#' @seealso [caret::dummyVars()]
#' @examples
#' data(iris)
#' x <- iris
#' one_hot(x)
#' 
#' @export
#' 
one_hot <- function(x, sep = ".") {
  if (is.matrix(x)) return(x)
  x <- droplevels(x)
  factor_ind <- index_factor(x, convert_bin = TRUE)
  bin_ind <- unlist(lapply(x, function(i) {
    !is.numeric(i) && nlevels(factor(i)) == 2
  }))
  if (sum(factor_ind) == 0) {
    out <- data.matrix(x)
    out[, bin_ind] <- out[, bin_ind] - 1
    return(out)
  }
  x0 <- data.matrix(x)
  x0[, bin_ind] <- x0[, bin_ind] - 1
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
