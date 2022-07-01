
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
#' @param ... Other arguments passed to [Boruta::Boruta]
#' @details
#' Boruta works differently from other filters in that it does not rank
#' variables by variable importance, but tries to determine relevant features
#' and divides features into Rejected, Tentative or Confirmed.
#' @importFrom Boruta Boruta
#' @export

boruta_filter <- function(y, x, select = c('Confirmed', 'Tentative'),
                           type = c("index", "names", "full"), ...) {
  type <- match.arg(type)
  ref <- Boruta::Boruta(x, y, ...)$finalDecision
  out <- which(ref %in% select)
  if (length(out) == 0) stop("No predictors left after filtering")
  switch(type,
         index = out,
         names = colnames(x)[out],
         full = ref)
}

