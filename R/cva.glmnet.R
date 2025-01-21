
#' Cross-validation of alpha for glmnet
#' 
#' Performs k-fold cross-validation for glmnet, including alpha mixing parameter.
#' 
#' @param x Matrix of predictors
#' @param y Response vector
#' @param nfolds Number of folds (default 10)
#' @param alphaSet Sequence of alpha values to cross-validate
#' @param foldid Optional vector of values between 1 and `nfolds` identifying
#'   what fold each observation is in.

#' @param parallel_method parallelization options "mclapply" (default) or "future"
#' @details With "future" option nested parallelisation is performed over `alphaSet` and then over
#'   the folds using the future framework. (Parallelization over lambda is not
#'   necessary due to the way glmnet works. See [glmnet::glmnet()].)

#' @param ... Other arguments passed to [glmnet::cv.glmnet]
#' @return Object of S3 class "cva.glmnet", which is a list of the cv.glmnet 
#' objects for each value of alpha and `alphaSet`.
#' \item{fits}{List of fitted [glmnet::cv.glmnet] objects}
#' \item{alphaSet}{Sequence of alpha values used}
#' \item{alpha_cvm}{The mean cross-validated error - a vector of length 
#' `length(alphaSet)`.}
#' \item{best_alpha}{Value of alpha giving lowest `alpha_cvm`.}
#' \item{which_alpha}{Index of `alphaSet` with lowest `alpha_cvm`}
#' @seealso [glmnet::cv.glmnet], [glmnet::glmnet]
#' @author Myles Lewis
#' @importFrom glmnet cv.glmnet
#' @importFrom utils tail
#' @export
#' 
cva.glmnet <- function(x, y, nfolds = 10, alphaSet = seq(0.1, 1, 0.1), parallel_method = "mclapply",
                       foldid = NULL, ...) {
  if (is.null(foldid)) {
    foldid <- sample(rep(seq_len(nfolds), length = NROW(y)))
  }
  if(parallel_method=="mclapply"){
    fit1 <- cv.glmnet(x = x, y = y, alpha = tail(alphaSet, 1), foldid = foldid, ...)
  }else{
 local_registerDoFuture()
  # Use run cv.glmnet inside a "useless" future_lapply here so that it is always
  # run at the same future parallel nesting level.
  fit1 <- future_lapply(tail(alphaSet, 1), function(alpha) {
    cv.glmnet(x = x, y = y,
              alpha = alpha, foldid = foldid, ..., parallel = TRUE)
  })[[1]]
  }
  if (length(alphaSet) > 1) {
    if(parallel_method=="mclapply"){
    fits <- lapply(alphaSet[1:(length(alphaSet)-1)], function(alpha) {
      cv.glmnet(x = x, y = y, alpha = alpha, foldid = foldid, lambda = fit1$lambda, ...)
    })
  }else{
      fits <- future_lapply(alphaSet[1:(length(alphaSet)-1)], function(alpha) {
      cv.glmnet(x = x, y = y, alpha = alpha, foldid = foldid, lambda = fit1$lambda, ..., parallel = TRUE)
  }    
    fits <- append(fits, list(fit1))
  } else fits <- list(fit1)
  if (fit1$name %in% c("AUC", "C-index")) {
    alpha_cvm <- unlist(lapply(fits, function(i) max(i$cvm)))
    which_alpha <- which.max(alpha_cvm)
  } else {
    alpha_cvm <- unlist(lapply(fits, function(i) min(i$cvm)))
    which_alpha <- which.min(alpha_cvm)
  }
  best_alpha <- alphaSet[which_alpha]
  cvafit <- list(fits = fits,
                 alphaSet = alphaSet,
                 alpha_cvm = alpha_cvm,
                 best_alpha = best_alpha,
                 which_alpha = which_alpha)
  class(cvafit) <- "cva.glmnet"
  cvafit
}


#' Extract coefficients from a cva.glmnet object
#' 
#' Extracts model coefficients from a fitted [cva.glmnet()] object.
#' 
#' @param object Fitted `cva.glmnet` object.
#' @param ... Other arguments passed to `coef.glmnet()` e.g. `s` the value of
#'   lambda at which coefficients are required.
#' @returns Sparse matrix containing coefficients from a `cv.glmnet` model
#' @export
coef.cva.glmnet <- function(object, ...) {
  coef(object$fits[[object$which_alpha]], ...)
}


#' Predict method for cva.glmnet models
#' 
#' Makes predictions from a cross-validated glmnet model with optimal value of 
#' lambda and alpha.
#' 
#' @param object Fitted `cva.glmnet` object.
#' @param newx Matrix of new values for `x` at which predictions are to be made.
#' @param s Value of penalty parameter lambda. Default value is `s="lambda.1se"`
#'   for consistency with glmnet. Alternatively `s="lambda.min"` can be used.
#' @param ... Other arguments passed to `predict.cv.glmnet()`.
#' @returns Object returned depends on arguments in `...` such as `type`.  
#' @export
predict.cva.glmnet <- function(object, newx,
                               s = "lambda.1se", ...) {
  w <- object$which_alpha
  fit <- object$fits[[w]]
  predict(fit, newx = as.matrix(newx), s = s, ...)
}


#' glmnet coefficients
#' 
#' Convenience function for retrieving coefficients from a [glmnet::cv.glmnet]
#' model at a specified lambda. Sparsity is removed and non-intercept
#' coefficients are ranked by absolute value.
#' 
#' @param fit A [glmnet::cv.glmnet] fitted model object.
#' @param s Value of lambda. See [glmnet::coef.glmnet] and
#'   [glmnet::predict.cv.glmnet]
#' @param ... Other arguments passed to [glmnet::coef.glmnet]
#' @return Vector or list of coefficients ordered with the intercept first, 
#' followed by highest absolute value to lowest.
#' @importFrom stats coef
#' @export
#' 
glmnet_coefs <- function(fit, s, ...) {
  if (length(fit) == 1 && is.na(fit)) return(NA)
  cf <- coef(fit, s = s, ...)
  if (is.list(cf)) {
    cf <- lapply(cf, function(i) {
      cf <- as.matrix(i)
      cf <- cf[cf != 0, ]
      cf2 <- cf[-1]
      cf2 <- cf2[order(abs(cf2), decreasing = TRUE)]
      c(cf[1], cf2)
    })
    return(cf)
  } 
  cf <- as.matrix(cf)
  cf <- cf[cf != 0, ]
  cf2 <- cf[-1]
  cf2 <- cf2[order(abs(cf2), decreasing = TRUE)]
  c(cf[1], cf2)  # keep intercept first
}


#' Extract coefficients from nestcv.glmnet object
#' 
#' Extracts coefficients from the final fit of a `"nestcv.glmnet"` object.
#' 
#' @param object Object of class `"nestcv.glmnet"`
#' @param s Value of penalty parameter lambda. Default is the mean of lambda 
#' values selected across each outer fold.
#' @param ... Other arguments passed to [glmnet::coef.glmnet()]
#' @return Vector or list of coefficients ordered with the intercept first, 
#' followed by highest absolute value to lowest.
#' @export
#'
coef.nestcv.glmnet <- function(object, s = object$final_param["lambda"], ...) {
  glmnet_coefs(object$final_fit, s = s, ...)
}


#' @export
print.nestcv.glmnet <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Nested cross-validation with glmnet\n")
  if (!is.null(x$call$filterFUN)) 
    cat("Filter: ", x$call$filterFUN, "\n") else cat("No filter\n")
  cat("\nFinal parameters:\n")
  if (length(x$final_param)==1 && is.na(x$final_param)) {
    cat("NA\n")
  } else print(x$final_param, digits = digits, print.gap = 2L)
  cat("\nFinal coefficients:\n")
  if (length(x$final_fit)==1 && is.na(x$final_fit)) {
    cat("NA\n")
  } else print(coef(x), digits = digits)
  cat("\nResult:\n")
  print(x$summary, digits = digits, print.gap = 2L)
}


#' @export
summary.nestcv.glmnet <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Nested cross-validation with glmnet\n")
  if (!is.null(object$call$filterFUN)) 
    cat("Filter: ", object$call$filterFUN, "\n") else cat("No filter\n")
  if (!is.null(object$call$modifyX))
    cat("Modifier: ", object$call$modifyX, "\n")
  cat("Outer loop: ", switch(object$outer_method,
                             cv = paste0(length(object$outer_folds), "-fold CV"),
                             LOOCV = "leave-one-out CV"))
  cat("\nInner loop: ", paste0(object$n_inner_folds, "-fold CV\n"))
  balance <- object$call$balance
  if (!is.null(balance)) {
    cat("Balancing: ", balance, "\n")
  }
  cat(object$dimx[1], "observations,", object$dimx[2], "predictors\n")
  if (!is.numeric(object$y)) print(c(table(object$y)))
  cat("\n")
  
  alpha <- unlist(lapply(object$outer_result, '[[', 'alpha'))
  lambda <- unlist(lapply(object$outer_result, '[[', 'lambda'))
  nfilter <- unlist(lapply(object$outer_result, '[[', 'nfilter'))
  foldres <- data.frame(alpha = alpha, lambda = lambda, n.filter = nfilter,
                        row.names = paste("Fold", seq_along(alpha)))
  
  print(foldres, digits = digits)
  cat("\nFinal parameters:\n")
  if (length(object$final_param)==1 && is.na(object$final_param)) {
    cat("NA\n")
  } else print(object$final_param, digits = digits, print.gap = 2L)
  cat("\nFinal coefficients:\n")
  if (length(object$final_fit)==1 && is.na(object$final_fit)) {
    cat("NA\n")
  } else print(coef(object), digits = digits)
  cat("\nResult:\n")
  print(object$summary, digits = digits, print.gap = 3L)
  out <- list(dimx = object$dimx, folds = foldres,
              final_param = object$final_param,
              coef = coef(object), result = object$summary)
  invisible(out)
}

#' Predict method for nestcv.glmnet fits
#'
#' Obtains predictions from the final fitted model from a [nestcv.glmnet]
#' object.
#' 
#' @param object Fitted `nestcv.glmnet` object
#' @param newdata New data to predict outcome on
#' @param s Value of lambda for glmnet prediction
#' @param modify Logical whether to modify `newdata` based on `modifyX`
#'   function. See `modifyX` and `modifyX_useY` arguments in [nestcv.glmnet()].
#' @param ... Other arguments passed to `predict.glmnet`.
#' @return Object returned depends on the `...` argument passed to predict
#'   method for `glmnet` objects.
#' @details Checks for missing predictors and if these are sparse (i.e. have
#'   zero coefficients) columns of 0 are automatically added to enable
#'   prediction to proceed.
#' @seealso [glmnet::glmnet]
#' @method predict nestcv.glmnet
#' @export
predict.nestcv.glmnet <- function(object, newdata,
                                  s = object$final_param["lambda"],
                                  modify = FALSE,
                                  ...) {
  if (modify) {
    if (is.null(object$modify_fit)) stop("`modify_fit` is missing")
    newdata <- predict(object$modify_fit, newdata)
  }
  newdata <- as.matrix(newdata)
  newx <- fix_cols(object$final_fit, newdata, s = s)
  predict(object$final_fit, newx = newx, s = unname(s), ...)
}


# fills in zero coefficent columns if missing
fix_cols <- function(x, newx, s) {
  cf <- coef(x, s = s)
  # check for multinomial
  if (is.list(cf)) {
    cf <- do.call(cbind, cf)
    cf <- as.matrix(cf)
  }
  final_vars <- rownames(cf)[-1]
  cf <- cf[-1,]
  # if full subset present
  if (all(final_vars %in% colnames(newx))) return(newx[, final_vars])
  # some cols missing
  nz <- if (is.vector(cf)) {
    cf != 0
  } else apply(cf, 1, function(i) any(i != 0))  # multinomial
  nonzeros <- final_vars[nz]
  if (!all(nonzeros %in% colnames(newx))) {
    stop("Some non-zero coefficient predictors are missing")}
  zeros <- final_vars[!nz]
  if (length(zeros) == 0) return(newx[, final_vars])
  miss <- zeros[!zeros %in% colnames(newx)] 
  if (length(miss) == 0) return(newx[, final_vars])
  mat <- matrix(0, nrow = nrow(newx), ncol = length(miss),
                dimnames = list(rownames(newx), miss))
  out <- cbind(newx, mat)
  out[, final_vars]  # col order seems to matter
}
