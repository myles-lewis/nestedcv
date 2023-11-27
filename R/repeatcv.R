
#' Repeated nested CV
#'
#' Performs repeated calls to a `nestedcv` model to determine performance across
#' repeated runs of nested CV.
#' 
#' @param expr An expression containing a call to [nestcv.glmnet()],
#'   [nestcv.train()], [nestcv.SuperLearner()] or [outercv()].
#' @param n Number of repeats
#' @param repeat_folds Optional list containing fold indices to be applied to
#'   the outer CV folds.
#' @param progress Logical whether to show progress.
#' @details
#' When comparing models, it is recommended to fix the sets of outer CV folds
#' used across each repeat for comparing performance between models. The
#' function [repeatfolds()] can be used to create a fixed set of outer CV folds
#' for each repeat.
#' @returns Matrix of performance metrics of S3 class 'repeatcv'.
#' @importFrom magrittr pipe_nested
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' \donttest{
#' data("iris")
#' dat <- iris
#' y <- dat$Species
#' x <- dat[, 1:4]
#'
#' res <- repeatcv(n = 3, nestcv.glmnet(y, x,
#'                                      family = "multinomial", alphaSet = 1,
#'                                      n_outer_folds = 4, cv.cores = 2))
#' res
#' 
#' ## using magrittr nested pipe
#' `%|>%` <- magrittr::pipe_nested
#' res <- nestcv.glmnet(y, x, family = "multinomial", alphaSet = 1,
#'                      n_outer_folds = 4, cv.cores = 2) %|>%
#'        repeatcv(3)
#' res
#' 
#' ## set up fixed fold indices
#' set.seed(123, "L'Ecuyer-CMRG")
#' folds <- repeatfolds(y, repeats = 3, n_outer_folds = 4)
#' res <- nestcv.glmnet(y, x, family = "multinomial", alphaSet = 1,
#'                      n_outer_folds = 4, cv.cores = 2) %|>%
#'        repeatcv(3, repeat_folds = folds)
#' res
#' }
#' @export

repeatcv <- function(expr, n = 5, repeat_folds = NULL, progress = TRUE) {
  cl <- match.call()
  if (!is.null(repeat_folds) && length(repeat_folds) != n)
    stop("mismatch between n and repeat_folds")
  start <- Sys.time()
  ex0 <- ex <- substitute(expr)
  # modify args in expr call
  ex$verbose <- FALSE
  d <- deparse(ex[[1]])
  if (d == "nestcv.glmnet" | d == "nestcv.train") ex$finalCV <- NA
  if (d == "nestcv.SuperLearner") ex$final <- FALSE
  if (progress) pb <- txtProgressBar(style = 3)
  out <- lapply(seq_len(n), function(i) {
    if (!is.null(repeat_folds)) ex$outer_folds <- repeat_folds[[i]]
    fit <- try(eval.parent(ex), silent = TRUE)
    if (progress) setTxtProgressBar(pb, i / n)
    if (inherits(fit, "try-error")) {
      if (progress) warning(fit[1])
      return(NA)
    }
    s <- fit$summary
    if (is.list(s)) return(s[[2]])
    s
  })
  
  if (progress) {
    close(pb)
    end <- Sys.time()
    cat("Time", format(end - start, digits = 3), "\n")
  }
  out <- do.call(rbind, out)
  rownames(out) <- seq_len(n)
  attr(out, "call") <- ex0
  class(out) <- c("repeatcv", class(out))
  out
}


#' Create folds for repeated nested CV
#' 
#' @param y Outcome vector
#' @param repeats Number of repeats
#' @param n_outer_folds Number of outer CV folds
#' @returns List containing indices of outer CV folds
#' @examples
#' \donttest{
#' data("iris")
#' dat <- iris
#' y <- dat$Species
#' x <- dat[, 1:4]
#' 
#' ## using magrittr nested pipe
#' `%|>%` <- magrittr::pipe_nested
#' 
#' ## set up fixed fold indices
#' set.seed(123, "L'Ecuyer-CMRG")
#' folds <- repeatfolds(y, repeats = 3, n_outer_folds = 4)
#' 
#' res <- nestcv.glmnet(y, x, family = "multinomial", alphaSet = 1,
#'                      n_outer_folds = 4, cv.cores = 2) %|>%
#'        repeatcv(3, repeat_folds = folds)
#' res
#' }
#' @export
#' 
repeatfolds <- function(y, repeats = 5, n_outer_folds = 10) {
  rfolds <- lapply(seq_len(repeats), function(i) createFolds(y, k = n_outer_folds))
  names(rfolds) <- paste0("Rep", seq_len(repeats))
  rfolds
}


#' @export
print.repeatcv <- function(x, ...) {
  class(x) <- class(x)[class(x) != "repeatcv"]
  attr(x, "call") <- NULL
  print(x)
}


#' @export
summary.repeatcv <- function(object, ...) {
  m <- colMeans(object, na.rm = TRUE)
  sd <- apply(object, 2, sd, na.rm = TRUE)
  sem <- sd / sqrt(nrow(object))
  data.frame(mean = m, sd = sd, sem = sem)
}
