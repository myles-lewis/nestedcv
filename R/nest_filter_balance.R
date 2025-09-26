
# nested filtering and oversampling for balancing

nest_filt_bal <- function(test, y, x,
                          filterFUN, filter_options,
                          balance = NULL, balance_options,
                          modifyX = NULL, modifyX_useY, modifyX_options,
                          penalty.factor = NULL) {
  if (is.null(test)) {
    ytrain <- y
    xtrain <- x
    ytest <- NULL
    xtest <- NULL
  } else {
    xtrain <- x[-test, , drop = FALSE]
    xtest <- x[test, , drop = FALSE]
    if (is.matrix(y)) {
      ytrain <- y[-test, , drop = FALSE]
      ytest <- y[test, , drop = FALSE]
    } else {
      ytrain <- y[-test]
      ytest <- y[test]
    }
  }
  
  if (is.none(filterFUN)) {
    filt_xtrain <- xtrain
    filt_xtest <- xtest
    filt_pen.factor <- penalty.factor
  } else {
    args <- list(y = ytrain, x = xtrain)
    args <- append(args, filter_options)
    fset <- do.call(filterFUN, args)
    filt_xtrain <- xtrain[, fset, drop = FALSE]
    filt_xtest <- xtest[, fset, drop = FALSE]
    filt_pen.factor <- penalty.factor[fset]
  }
  
  if (!is.none(modifyX)) {
    if (!modifyX_useY) {
      # only modify X
      args <- list(x = filt_xtrain)
      args <- append(args, modifyX_options)
      filt_xtrain <- do.call(modifyX, args)
      if (!is.null(test)) {
        args <- list(x = filt_xtest)
        args <- append(args, modifyX_options)
        filt_xtest <- do.call(modifyX, args)
      }
    } else {
      # modify X using information from ytrain
      args <- list(y = ytrain, x = filt_xtrain)
      args <- append(args, modifyX_options)
      modfit <- do.call(modifyX, args)
      filt_xtrain <- predict(modfit, newdata = filt_xtrain)
      if (!is.null(test)) {
        filt_xtest <- predict(modfit, newdata = filt_xtest)
      }
    }
    if (!is.null(test) && !identical(colnames(filt_xtrain), colnames(filt_xtest)))
      message("Error in modifyX: different colnames in xtrain and xtest")
    if (!is.null(penalty.factor) && !identical(colnames(x), colnames(filt_xtrain))) {
      pf <- rep(1, ncol(filt_xtrain))
      ind <- match(colnames(x), colnames(filt_xtrain))
      ok <- !is.na(ind)
      pf[ind[ok]] <- penalty.factor[ok]
      filt_pen.factor <- pf
    }
  }
  
  if (!is.none(balance)) {
    args <- list(y = ytrain, x = filt_xtrain)
    args <- append(args, balance_options)
    bal_dat <- do.call(balance, args)
    ytrain <- bal_dat$y
    filt_xtrain <- bal_dat$x
  }
  
  out <- list(ytrain = ytrain, ytest = ytest,
              filt_xtrain = filt_xtrain, filt_xtest = filt_xtest,
              filt_pen.factor = filt_pen.factor)
  if (!is.null(modifyX) & modifyX_useY) out$modify_fit <- modfit
  out
}


is.none <- function(x) {
  is.null(x) || (is.character(x) && x == "none")
}


#' Calculate weights for class imbalance
#' 
#' @param y Factor or character response vector. If a character vector is
#'   supplied it is coerced into a factor. Unused levels are dropped.
#' @return Vector of weights
#' @export
#' 
weight <- function(y) {
  if (is.numeric(y)) {
    message("y is numeric: this function is designed for classification")
  }
  y <- factor(y)
  tab <- c(table(y))
  props <- 1/tab
  weights <- props[as.numeric(y)]
  weights <- weights / mean(weights, na.rm = TRUE)
  weights
}

