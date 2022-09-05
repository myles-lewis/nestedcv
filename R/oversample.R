
#' Oversampling and undersampling
#' 
#' Random oversampling to compensate for data imbalance
#' 
#' @param y Vector of response outcome as a factor
#' @param over Amount of oversampling of the minority class. If set to `NULL`
#'   then all classes will be oversampled up to the number of samples in the
#'   majority class. To turn off oversampling set `over = 1`.
#' @param under Amount of undersampling of the majority class
#' @param yminor Optional character value specifying the level in `y` which is
#'   to be oversampled. If `NULL`, this is set automatically to the class with
#'   the smallest sample size.
#' @return A vector of sample indices
#' @export
#' 
oversample <- function(y, over = NULL, under = 1, yminor = NULL) {
  ytab <- table(y)
  ymajor <- names(ytab)[which.max(ytab)]
  nymajor <- round(max(ytab) * under)
  if (is.null(over)) {
    # equalise
    yset <- names(ytab)[!names(ytab) %in% ymajor]
    add_samples <- unlist(lapply(yset, function(i) {
      size <- nymajor - ytab[i]
      ind <- which(y == i)
      sample(ind, size, replace = size > length(ind))
    }))
  } else if (over > 1) {
    # manual
    if (is.null(yminor)) yminor <- names(ytab)[which.min(ytab)]
    size <- round(ytab[yminor] * (over - 1))
    ind <- which(y == yminor)
    add_samples <- sample(ind, size, replace = (over >= 2))
  } else add_samples <- NULL
  
  if (under != 1) {
    ind <- which(y == ymajor)
    size <- round(under * length(ind))
    major_samples <- sample(ind, size, replace = size > length(ind))
    rm_samples <- ind[!ind %in% major_samples]
    out <- c(seq_along(y)[-rm_samples], add_samples)
  } else {
    out <- c(seq_along(y), add_samples)
  }
  
  out
}

