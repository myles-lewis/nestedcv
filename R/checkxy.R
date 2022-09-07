
checkxy <- function(y, x, na.option) {
  if (length(y) != nrow(x))
    stop("Mismatch in length of 'y' and number of rows in 'x'", call. = FALSE)
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
  list(r = okr, c = okc)
}
