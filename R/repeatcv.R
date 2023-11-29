
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
#' @param keep Logical whether to save repeated outer CV predictions for ROC
#'   curves etc.
#' @param progress Logical whether to show progress.
#' @details
#' When comparing models, it is recommended to fix the sets of outer CV folds
#' used across each repeat for comparing performance between models. The
#' function [repeatfolds()] can be used to create a fixed set of outer CV folds
#' for each repeat.
#' @returns List of S3 class 'repeatcv' containing the model call, matrix of
#'   performance metrics, and if `keep = TRUE` a matrix or dataframe containing
#'   the outer CV predictions from each repeat.
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
#' summary(res)
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

repeatcv <- function(expr, n = 5, repeat_folds = NULL, keep = FALSE, progress = TRUE) {
  cl <- match.call()
  if (!is.null(repeat_folds) && length(repeat_folds) != n)
    stop("mismatch between n and repeat_folds")
  ex0 <- ex <- substitute(expr)
  # modify args in expr call
  ex$verbose <- FALSE
  d <- deparse(ex[[1]])
  if (d == "nestcv.glmnet" | d == "nestcv.train") ex$finalCV <- NA
  if (d == "nestcv.SuperLearner") ex$final <- FALSE
  if (d == "nestcv.train") d <- ex$method
  d <- gsub("nestcv.", "", d)
  if (progress) pb <- txtProgressBar2(title = d)
  res <- lapply(seq_len(n), function(i) {
    if (!is.null(repeat_folds)) ex$outer_folds <- repeat_folds[[i]]
    fit <- try(eval.parent(ex), silent = TRUE)
    if (progress) setTxtProgressBar(pb, i / n)
    if (inherits(fit, "try-error")) {
      if (progress) warning(fit[1])
      if (keep) return(list(NA, NA))
      return(NA)
    }
    s <- fit$summary
    if (is.list(s)) s <- s[[2]]
    if (keep) return(list(s, fit$output))
    s
  })
  if (progress) close(pb)
  
  if (keep) {
    res1 <- lapply(res, "[[", 1)
    result <- do.call(rbind, res1)
    rownames(result) <- seq_len(n)
    res2 <- lapply(res, "[[", 2)
    output <- do.call(rbind, res2)
    out <- list(call = ex0, result = result, output = output)
  } else {
    result <- do.call(rbind, res)
    rownames(result) <- seq_len(n)
    out <- list(call = ex0, result = result)
  }
  class(out) <- c("repeatcv")
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
summary.repeatcv <- function(object,
                             digits = max(3L, getOption("digits") - 3L),
                             ...) {
  cat("Call:\n")
  print(object$call)
  cat(nrow(object$result), "repeats\n")
  m <- colMeans(object$result, na.rm = TRUE)
  sd <- apply(object$result, 2, sd, na.rm = TRUE)
  sem <- sd / sqrt(nrow(object$result))
  df <- data.frame(mean = m, sd = sd, sem = sem)
  print(df, digits = digits)
  invisible(df)
}
