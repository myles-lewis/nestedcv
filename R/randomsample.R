
#' Oversampling and undersampling
#' 
#' Random oversampling of the minority group(s) or undersampling of the majority
#' group to compensate for class imbalance in datasets.
#' 
#' @param y Vector of response outcome as a factor
#' @param x Matrix of predictors
#' @param minor Amount of oversampling of the minority class. If set to `NULL`
#'   then all classes will be oversampled up to the number of samples in the
#'   majority class. To turn off oversampling set `minor = 1`.
#' @param major Amount of undersampling of the majority class
#' @param yminor Optional character value specifying the level in `y` which is
#'   to be oversampled. If `NULL`, this is set automatically to the class with
#'   the smallest sample size.
#' @details
#' `minor` < 1 and `major` > 1 are ignored.
#' @return List containing extended matrix `x` of synthesised data and extended
#'   response vector `y`
#' @examples
#' \donttest{
#' ## Imbalanced dataset
#' set.seed(1, "L'Ecuyer-CMRG")
#' x <- matrix(rnorm(150 * 2e+04), 150, 2e+04)  #' predictors
#' y <- factor(rbinom(150, 1, 0.2))  #' imbalanced binary response
#' table(y)
#' 
#' ## first 30 parameters are weak predictors
#' x[, 1:30] <- rnorm(150 * 30, 0, 1) + as.numeric(y)*0.5
#' 
#' ## Balance x & y outside of CV loop by random oversampling minority group
#' out <- randomsample(y, x)
#' y2 <- out$y
#' x2 <- out$x
#' table(y2)
#' 
#' ## Nested CV glmnet with unnested balancing by random oversampling on
#' ## whole dataset
#' fit1 <- nestcv.glmnet(y2, x2, family = "binomial", alphaSet = 1,
#'                       cv.cores=2,
#'                       filterFUN = ttest_filter)
#' fit1$summary
#' 
#' ## Balance x & y outside of CV loop by random oversampling minority group
#' out <- randomsample(y, x, minor=1, major=0.4)
#' y2 <- out$y
#' x2 <- out$x
#' table(y2)
#' 
#' ## Nested CV glmnet with unnested balancing by random undersampling on
#' ## whole dataset
#' fit1b <- nestcv.glmnet(y2, x2, family = "binomial", alphaSet = 1,
#'                        cv.cores=2,
#'                        filterFUN = ttest_filter)
#' fit1b$summary
#' 
#' ## Balance x & y outside of CV loop by SMOTE
#' out <- smote(y, x)
#' y2 <- out$y
#' x2 <- out$x
#' table(y2)
#' 
#' ## Nested CV glmnet with unnested balancing by SMOTE on whole dataset
#' fit2 <- nestcv.glmnet(y2, x2, family = "binomial", alphaSet = 1,
#'                       cv.cores=2,
#'                       filterFUN = ttest_filter)
#' fit2$summary
#' 
#' ## Nested CV glmnet with nested balancing by random oversampling
#' fit3 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1,
#'                       cv.cores=2,
#'                       balance = "randomsample",
#'                       filterFUN = ttest_filter)
#' fit3$summary
#' class_balance(fit3)
#' 
#' ## Plot ROC curves
#' plot(fit1$roc, col='green')
#' lines(fit1b$roc, col='red')
#' lines(fit2$roc, col='blue')
#' lines(fit3$roc)
#' legend('bottomright', legend = c("Unnested random oversampling", 
#'                                  "Unnested SMOTE",
#'                                  "Unnested random undersampling",
#'                                  "Nested balancing"), 
#'        col = c("green", "blue", "red", "black"), lty=1, lwd=2)
#' }
#' 
#' @export
#' 
randomsample <- function(y, x, minor = NULL, major = 1, yminor = NULL) {
  ytab <- table(y)
  ymajor <- names(ytab)[which.max(ytab)]
  nymajor <- round(max(ytab) * major)
  if (is.null(minor)) {
    # equalise
    yset <- names(ytab)[!names(ytab) %in% ymajor]
    add_samples <- unlist(lapply(yset, function(i) {
      size <- nymajor - ytab[i]
      ind <- which(y == i)
      sample(ind, size, replace = size > length(ind))
    }))
  } else if (minor > 1) {
    # manual
    if (is.null(yminor)) yminor <- names(ytab)[which.min(ytab)]
    size <- round(ytab[yminor] * (minor - 1))
    ind <- which(y == yminor)
    add_samples <- sample(ind, size, replace = (minor >= 2))
  } else add_samples <- NULL
  
  if (major < 1) {
    ind <- which(y == ymajor)
    size <- round(major * length(ind))
    major_samples <- sample(ind, size, replace = size > length(ind))
    rm_samples <- ind[!ind %in% major_samples]
    out <- c(seq_along(y)[-rm_samples], add_samples)
  } else {
    out <- c(seq_along(y), add_samples)
  }
  
  y <- y[out]
  x <- x[out, ]
  rownames(x) <- make.names(rownames(x), unique = TRUE)
  list(y = y, x = x)
}


sample_message <- function(minor = 1, major = 1) {
  if (is.null(minor)) {
    cat("Random oversampling to equalise classes\n")
  } else if (minor > 1) cat("Random oversampling minority x", minor, "\n")
  if (!is.null(major)) {
    if (major < 1) cat("Random undersampling majority x", major, "\n")
  }
}

