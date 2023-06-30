

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
  x1 <- x0[, !factor_ind]
  f <- which(factor_ind)
  res_x2 <- lapply(f, function(i) {
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
  res_x2 <- do.call(cbind, res_x2)
  out <- cbind(x1, res_x2)
  out
}
