
#' SMOTE
#' 
#' Synthetic Minority Oversampling Technique (SMOTE) algorithm for imbalanced
#' classification data.
#' 
#' @param y Vector of response outcome as a factor
#' @param x Matrix of predictors
#' @param k Range of KNN to consider for generation of new data
#' @param over Amount of oversampling of the minority class. If set to `NULL`
#'   then all classes will be oversampled up to the number of samples in the
#'   majority class.
#' @param yminor Optional character value specifying the level in `y` which is
#'   to be oversampled. If `NULL`, this is set automatically to the class with
#'   the smallest sample size.
#' @return List containing extended matrix `x` of synthesised data and extended
#'   response vector `y`
#' @importFrom stats dist runif
#' @export

smote <- function(y, x, k = 5, over = NULL, yminor = NULL) {
  stopifnot(k >= 1)
  y <- droplevels(y)
  ytab <- table(y)
  
  if (is.null(over)) {
    # equalise
    ymajor <- names(ytab)[which.max(ytab)]
    n_ymajor <- round(max(ytab))  # number to sample up to
    yset <- names(ytab)[!names(ytab) %in% ymajor]
    newx <- lapply(yset, function(i) {
      xminor <- x[y == i, ]
      n <- n_ymajor - nrow(xminor)
      if (n == 0) return(NULL)
      smoteN(xminor, k, n)
    })
    newx <- do.call(rbind, newx)
    n_yset <- n_ymajor - ytab[names(ytab) != ymajor]
    newy <- unlist(lapply(seq_along(yset), function(i) {
      rep(yset[i], n_yset[i])
    }))
  } else {
    # single round
    if (is.null(yminor)) yminor <- names(ytab)[which.min(ytab)]
    xminor <- x[y == yminor, ]
    n <- round(nrow(xminor) * (over -1))
    newx <- smoteN(xminor, k, n)
    newy <- rep(yminor, n)
  }
  
  x <- rbind(x, newx)
  y <- c(y, factor(newy))
  list(y = y, x = x)
}
  

smoteN <- function(xminor, k, n) {
  nmin <- nrow(xminor)
  d <- as.matrix(dist(xminor))
  d_order <- apply(d, 1, order)
  knn <- t(d_order[2:(k+1), , drop = FALSE])  # samples in rows
  
  s1 <- smote_sample(nmin, n)
  s2 <- sample.int(k, n, replace = TRUE)
  out <- vapply(1:n, function(i) {
    x1 <- xminor[s1[i], ]
    x2 <- xminor[knn[s1[i], s2[i]], ]
    r <- runif(1)
    x1 * r + x2 * (1-r)
  }, numeric(ncol(xminor)))
  colnames(out) <- paste0(rownames(xminor)[s1], ".sm")
  t(out)
}


smote_sample <- function(n, size) {
  nfull <- size %/% n
  rem <- size %% n
  sam <- if (nfull > 0) rep(seq_len(n), nfull) else c()
  if (rem > 0) sam <- c(sam, sample.int(n, rem, replace = FALSE))
  sam
}

