
#' Oversampling and undersampling
#' 
#' Random oversampling of the minority group(s) or undersampling of the majority
#' group to compensate for class imbalance in datasets.
#' 
#' @param y Vector of response outcome as a factor
#' @param minor Amount of oversampling of the minority class. If set to `NULL`
#'   then all classes will be oversampled up to the number of samples in the
#'   majority class. To turn off oversampling set `minor = 1`.
#' @param major Amount of undersampling of the majority class
#' @param yminor Optional character value specifying the level in `y` which is
#'   to be oversampled. If `NULL`, this is set automatically to the class with
#'   the smallest sample size.
#' @param verbose Logical whether to display messages
#' @details
#' `minor` < 1 and `major` > 1 are ignored.
#' @return A vector of sample indices
#' @export
#' 
randomsample <- function(y, minor = NULL, major = 1, yminor = NULL,
                         verbose = TRUE) {
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
  
  if (verbose) {
    sample_message(minor, major)
    print(c(table(y[out])))
  }
  
  out
}


sample_message <- function(minor = 1, major = 1) {
  if (is.null(minor)) {
    cat("Random oversampling to equalise classes\n")
  } else if (minor > 1) cat("Random oversampling minority x", minor, "\n")
  if (!is.null(major)) {
    if (major < 1) cat("Random undersampling majority x", major, "\n")
  }
}
