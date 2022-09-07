
#' SMOTE
#' 
#' Synthetic Minority Oversampling Technique (SMOTE) algorithm for imbalanced
#' classification data.
#' 
#' @param y Vector of response outcome as a factor
#' @param x Matrix of predictors
#' @param k 
#' @param over Amount of oversampling of the minority class. If set to `NULL`
#'   then all classes will be oversampled up to the number of samples in the
#'   majority class. To turn off oversampling set `minor = 1`.
#' @param yminor Optional character value specifying the level in `y` which is
#'   to be oversampled. If `NULL`, this is set automatically to the class with
#'   the smallest sample size.
#' @return Matrix of synthesised data
#' @importFrom stats dist runif
#' @export

smote <- function(y, x, k = 5, over = NULL, yminor = NULL) {
  stopifnot(k >= 1)
  ytab <- table(y)
  if (is.null(yminor)) yminor <- names(ytab)[which.min(ytab)]
  xminor <- x[y == yminor, ]
  nmin <- nrow(xminor)
  
  if (is.null(over)) {
    ymajor <- names(ytab)[which.max(ytab)]
    nymajor <- round(max(ytab))  # number to sample up to
    n <- nymajor - nmin
  } else n <- round(nmin * (over -1))
  
  d <- as.matrix(dist(xminor))
  d_order <- apply(d, 1, order)
  knn <- t(d_order[2:(k+1), , drop = FALSE])  # samples in rows
  
  s1 <- sample.int(nmin, n, replace = (n > nmin))
  s2 <- sample.int(k, n, replace = TRUE)
  out <- vapply(1:n, function(i) {
    x1 <- xminor[s1[i], ]
    x2 <- xminor[knn[s1[i], s2[i]], ]
    r <- runif(1)
    x1 * r + x2 * (1-r)
  }, numeric(ncol(x)))
  colnames(out) <- paste0(rownames(xminor)[s1], ".sm")
  t(out)
}

