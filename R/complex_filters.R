
#' Boruta filter
#'
#' Filter using Boruta algorithm.
#'
#' @param y Response vector
#' @param x Matrix of predictors
#' @param select Which type of features to retain. Options include "Confirmed"
#'   and/or "Tentative".
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names, "full" returns a named vector of variable 
#' importance.
#' @param ... Other arguments passed to [Boruta::Boruta()]
#' @details
#' Boruta works differently from other filters in that it does not rank
#' variables by variable importance, but tries to determine relevant features
#' and divides features into Rejected, Tentative or Confirmed.
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters. If
#'   `type` is `"full"` full output from `Boruta` is returned.
#' 
#' @export

boruta_filter <- function(y, x, select = c('Confirmed', 'Tentative'),
                           type = c("index", "names", "full"), ...) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    stop("Package 'Boruta' must be installed", call. = FALSE)
  }
  type <- match.arg(type)
  ref <- Boruta::Boruta(x, y, ...)$finalDecision
  out <- which(ref %in% select)
  if (length(out) == 0) stop("No predictors left after filtering")
  switch(type,
         index = out,
         names = colnames(x)[out],
         full = ref)
}


#' Multilayer filter
#' 
#' Experimental filter designed for use with imbalanced datasets. Each round a
#' simple t-test is used to rank predictors and keep a certain number. After
#' each round a set number of cases are culled determined as the most outlying
#' cases - those which if used as a cutoff for classification have the smallest
#' number of misclassified cases. The t-test is repeated on the culled dataset
#' so that after successive rounds the most influential outlying samples have
#' been removed and different samples drive the t-test filter.
#' 
#' @param y Response vector
#' @param x Matrix of predictors
#' @param nfilter Vector of number of target predictors to keep at each round.
#'   The length of this vector determines the number of rounds of culling.
#' @param imbalance Logical whether to assume the dataset is imbalanced, in
#'   which case samples are only culled from the majority class.
#' @param cull number of samples to cull at each round
#' @param force_vars not implemented yet
#' @param verbose whether to show sample IDs of culled individuals at each round
#' @param type Type of vector returned. Default "index" returns indices,
#' "names" returns predictor names.
#' 
#' @return Integer vector of indices of filtered parameters (type = "index") or
#'   character vector of names (type = "names") of filtered parameters.
#' @export

layer_filter <- function(y, x,
                         nfilter = NULL,
                         imbalance = TRUE,
                         cull = 5,
                         force_vars = NULL,
                         verbose = FALSE,
                         type = c("index", "names", "full")) {
  type <- match.arg(type)
  if (imbalance) {
    tab <- table(y)
    maj_class <- names(tab)[which.max(tab)]
  }
  x <- as.matrix(x)
  out <- NULL
  
  for (nf in nfilter) {
    tt <- ttest_filter(y = y, x = x, nfilter = nf, p_cutoff = NULL,
                       type = "full")
    tt <- tt[order(tt[, 'pvalue']), ]
    tt <- tt[1:nf, ]
    maj_index <- y == maj_class
    min_index <- !maj_index
    find_clean <- sapply(rownames(tt), function(i) {
      xset <- x[maj_index, i]
      if (tt[i, 'stat'] > 0) {
        out <- sapply(xset, function(xcut) sum(x[min_index, i] > xcut))
      } else {
        out <- sapply(xset, function(xcut) sum(x[min_index, i] < xcut))
      }
      out
    })
    cleansum <- rowSums(find_clean)
    cullset <- names(cleansum)[order(cleansum)[1:cull]]
    if (verbose) print(cullset)
    ok <- !rownames(x) %in% cullset
    y <- y[ok]
    x <- x[ok,]
    out <- c(out, rownames(tt))
  }
  if (type == "index") return(which(colnames(x) %in% unique(out)))
  if (type == "names") return(unique(out))
  out
}


pls_filter <- function(y, x,
                       force_vars = NULL,
                       nfilter,
                       ncomp = 5,
                       scale_x = TRUE,
                       type = c("index", "names", "full"), ...) {
  type <- match.arg(type)
  if (is.factor(y) && nlevels(y) > 2) stop("Classes > 2 not supported")
  y <- as.numeric(y)
  if (scale_x) x <- scale(x)
  fit <- pls::plsr(y ~ x, ncomp = ncomp, ...)
  cf <- fit$coefficients
  cf <- lapply(seq_len(ncomp), function(i) {
    cfi <- cf[,, i]
    cfi[order(abs(cfi), decreasing = TRUE)]
  })
  if (type == "full") return(cf)
  
  if (length(nfilter) == 1) {
    # find sufficient vars from each comp
    topvars <- ""
    n <- floor(nfilter / ncomp)
    while (length(topvars) < nfilter) {
      topvars <- unique(unlist(lapply(seq_len(ncomp), function(i) {
        names(cf[[i]][seq_len(n)])
      })))
      n <- n +1
    }
    topvars <- topvars[seq_len(nfilter)]
  } else {
    # nfilter as vector
    if (length(nfilter) != ncomp) stop("nfilter is not the same length as ncomp")
    topvars <- unique(unlist(lapply(seq_len(ncomp), function(i) {
      names(cf[[i]][1:nfilter[i]])
    })))
  }
  
  topvars <- unique(c(topvars, force_vars))
  if (type == "names") return(topvars)
  which(colnames(x) %in% topvars)
}

