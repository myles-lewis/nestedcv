
# nested filtering and oversampling for balancing

nest_filt_bal <- function(test, y, x,
                          filterFUN, filter_options,
                          balance = NULL, balance_options,
                          penalty.factor = NULL) {
  if (is.null(test)) {
    ytrain <- y
    xtrain <- x
    ytest <- NULL
    xtest <- NULL
  } else {
    ytrain <- y[-test]
    xtrain <- x[-test, , drop = FALSE]
    ytest <- y[test]
    xtest <- x[test, , drop = FALSE]
  }
  
  # if (!is.null(balance)) {
  #   args <- list(y = ytrain, x = xtrain)
  #   args <- append(args, balance_options)
  #   bal_dat <- do.call(balance, args)
  #   ytrain <- bal_dat$y
  #   xtrain <- bal_dat$x
  # }
  
  if (is.null(filterFUN)) {
    filt_xtrain <- xtrain
    filt_xtest <- xtest
    filt_pen.factor <- penalty.factor
  } else {
    args <- list(y = ytrain, x = xtrain)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    filt_xtrain <- xtrain[, fset]
    filt_xtest <- xtest[, fset, drop = FALSE]
    filt_pen.factor <- penalty.factor[fset]
  }
  
  if (!is.null(balance)) {
    args <- list(y = ytrain, x = filt_xtrain)
    args <- append(args, balance_options)
    bal_dat <- do.call(balance, args)
    ytrain <- bal_dat$y
    filt_xtrain <- bal_dat$x
  }
  
  list(ytrain = ytrain, ytest = ytest,
       filt_xtrain = filt_xtrain, filt_xtest = filt_xtest,
       filt_pen.factor = filt_pen.factor)
}


#' Calculate weights for class imbalance
#' 
#' @param y Response vector
#' @return Vector of weights
#' @export
#' 
weight <- function(y) {
  if (is.numeric(y)) {
    message("y is numeric: this function is designed for classification")
  }
  tab <- c(table(y))
  props <- 1/tab
  weights <- props[as.numeric(y)]
  weights <- weights / sum(weights)
  weights
}

