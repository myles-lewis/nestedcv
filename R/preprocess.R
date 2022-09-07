
#' Preprocess data
#' 
#' Preprocessing of datasets including checking and optionally omitting `NA` and
#' handling of class imbalance.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param balance Specifies method for dealing with imbalanced data. Current
#'   options are `"none"` or `"randomsample"`. See [randomsample()]
#' @param balance_options List of additional arguments passed to the balancing
#'   function
#' @param na.option Character value specifying how `NA`s are dealt with.
#'   `"omit"` (the default) is equivalent to `na.action = na.omit`. `"omitcol"`
#'   removes cases if there are `NA` in 'y', but columns (predictors) containing
#'   `NA` are removed from 'x' to preserve cases. Any other value means that
#'   `NA` are ignored (a message is given).
#' @return A list containing `y`, `x` and `pre` which is a list containing
#'   balancing information and original values of `y`
#' @export

preprocess <- function(y, x,
                       balance = c("none", "randomsample", "smote"),
                       balance_options = list(minor = NULL, major = 1),
                       na.option = "pass") {
  balance <- match.arg(balance)
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
  y <- y[okr]
  x <- x[okr, okc]
  if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
  
  orig_y <- y
  pre <- list(balance = balance)
  if (!is.numeric(y)) {
    pre <- append(pre, list(balance_options = balance_options, orig_y = orig_y))
    if (balance == "randomsample") {
      args <- list(y = y)
      args <- append(args, balance_options)
      samples <- do.call(randomsample, args)
      y <- y[samples]
      x <- x[samples,]
      rownames(x) <- make.names(rownames(x), unique = TRUE)
      pre <- append(pre, samples)
    } else if (balance == "smote") {
      args <- list(y = y, x = x)
      args <- append(args, balance_options)
      smotedata <- do.call(smote, args)
      y <- smotedata$y
      x <- smotedata$x
    }
  }
  
  list(y = y, x = x, pre = pre)
}

