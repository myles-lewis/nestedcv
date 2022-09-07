
# nested filtering and oversampling for balancing

nest_filt_bal <- function(test, y, x,
                          filterFUN, filter_options,
                          balance, balance_options,
                          penalty.factor = NULL) {
  ytrain <- y[-test]
  xtrain <- x[-test, ]
  ytest <- y[test]
  xtest <- x[test, ]
  
  if (!is.null(balance)) {
    args <- list(y = ytrain, x = xtrain)
    args <- append(args, balance_options)
    bal_dat <- do.call(balance, args)
    ytrain <- bal_dat$y
    xtrain <- bal_dat$x
  }
  
  if (is.null(filterFUN)) {
    filt_xtrain <- xtrain
    filt_xtest <- xtest
    filt_pen.factor <- penalty.factor
  } else {
    args <- list(y = ytrain, x = xtrain)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    filt_xtrain <- xtrain[, fset]
    filt_xtest <- xtest[, fset]
    filt_pen.factor <- penalty.factor[fset]
  }
  
  list(ytrain = ytrain, ytest = ytest,
       filt_xtrain = filt_xtrain, filt_xtest = filt_xtest,
       filt_pen.factor = filt_pen.factor)
}
