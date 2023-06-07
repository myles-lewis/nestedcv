
checkxy <- function(y, x, na.option = "pass", weights = NULL) {
  if (length(y) != nrow(x))
    stop("Mismatch in length of 'y' and number of rows in 'x'", call. = FALSE)
  if (!is.null(weights)) {
    if (length(y) != length(weights)) stop("Mismatch in 'y' and 'weights'")
    if (any(weights < 0)) stop("Weights must be positive")
  }
  nay <- is.na(y)
  if (any(nay)) message("'y' contains ", sum(nay), " NA")
  naxr <- !complete.cases(x)
  naxc <- !complete.cases(t(x))
  okr <- rep.int(TRUE, length(y))
  okc <- rep.int(TRUE, ncol(x))
  if (na.option == "omit") {
    if (sum(naxr) == 1) {message("1 row in 'x' contains NA")
    } else if (sum(naxr) > 1) message(sum(naxr), " rows in 'x' contain NA")
    okr <- !nay & !naxr
  } else {
    if (sum(naxc) == 1) {message("1 column in 'x' contains NA")
    } else if (sum(naxc) > 1) message(sum(naxc), " columns in 'x' contain NA")
  }
  if (na.option == "omitcol") {
    okr <- !nay
    okc <- !naxc
  }
  x <- data.matrix(x)
  sc <- scale(x)
  sds <- attr(sc, "scaled:scale")  # Rfast & matrixStats have bugs if var=0
  var0 <- sds == 0
  if (any(var0)) {
    message(sum(var0), " predictor(s) have var=0")
    okc <- okc & !var0
  }
  list(r = okr, c = okc)
}
