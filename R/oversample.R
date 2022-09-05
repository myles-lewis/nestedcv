


replicateData <- function(y, x, inc = 2, ymin = NULL) {
  
}



#' @export
#' 
oversample <- function(y, over = NULL, yminor = NULL, under = 1) {
  ytab <- table(y)
  ymajor <- names(ytab)[which.max(ytab)]
  if (is.null(over)) {
    # equalise
    yset <- names(ytab)[!names(ytab) %in% ymajor]
    add_samples <- unlist(lapply(yset, function(i) {
      size <- max(ytab) - ytab[i]
      ind <- which(y == i)
      sample(ind, size, replace = size > length(ind))
    }))
  } else if (over > 1) {
    # manual
    if (is.null(yminor)) yminor <- names(ytab)[which.min(ytab)]
    size <- round(ytab[yminor] * (over - 1))
    ind <- which(y == yminor)
    add_samples <- sample(ind, size, replace = (over >= 2))
  }
  
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

