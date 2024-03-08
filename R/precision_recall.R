
#' Build precision-recall curve
#' 
#' Builds a precision-recall curve for a 'nestedcv' model using `prediction()`
#' and `performance()` functions from the ROCR package and returns an object of
#' class 'prc' for plotting.
#' 
#' @param response binary factor vector of response of default order controls,
#'   cases.
#' @param predictor numeric vector of probabilities
#' @param output data.frame with columns `testy` containing observed response
#'   from test folds, and `predyp` predicted probabilities for classification
#' @param object a 'nestcv.glmnet', 'nestcv.train', 'nestcv.SuperLearn',
#'   'outercv' or 'repeatcv' S3 class results object.
#' @param ... other arguments
#' @returns An object of S3 class 'prc' containing the following fields:
#' \item{recall}{vector of recall values}
#' \item{precision}{vector of precision values}
#' \item{auc}{area under precision-recall curve value using trapezoid method}
#' @export
prc <- function(...) {
  UseMethod("prc")
}

#' @rdname prc
#' @export
prc.default <- function(response, predictor, ...) {
  pred_obj <- ROCR::prediction(predictor, response)
  perf_obj <- ROCR::performance(pred_obj, measure = "prec", x.measure = "rec")
  x <- perf_obj@x.values[[1]]
  y <- perf_obj@y.values[[1]]
  auc <- auc_calc(x, y)
  out <- list(recall = x, precision = y, auc = auc)
  class(out) <- "prc"
  out
}

#' @rdname prc
#' @export
prc.data.frame <- function(output, ...) {
  if (!all(c("testy", "predyp") %in% colnames(output)))
    stop("not a classification output dataframe")
  prc.default(output$testy, output$predyp)
}

#' @rdname prc
#' @export
prc.nestcv.glmnet <- function(object, ...) {
  prc.data.frame(object$output)
}

#' @rdname prc
#' @export
prc.nestcv.train <- function(object, ...) {
  prc.data.frame(object$output)
}

#' @rdname prc
#' @export
prc.nestcv.SuperLearner <- function(object, ...) {
  prc.data.frame(object$output)
}

#' @rdname prc
#' @export
prc.outercv <- function(object, ...) {
  prc.data.frame(object$output)
}

#' @rdname prc
#' @export
prc.repeatcv <- function(object, ...) {
  prc.data.frame(object$output)
}

auc_calc <- function(x, y) {
  if (is.unsorted(x)) {
    y <- y[order(x)]
    x <- x[order(x)]
  }
  sum((y[-1] + y[-length(y)]) /2 * diff(x), na.rm = TRUE)
}



#' Plot precision-recall curve
#' 
#' Plots a precision-recall curve using base graphics. It accepts an S3 object
#' of class 'prc', see [prc()].
#' 
#' @param x An object of class 'prc'
#' @param ... Optional graphical arguments passed to [plot()]
#' @return No return value
#' @seealso [prc()]
#' @export
plot.prc <- function(x, ...) {
  new.args <- list(...)
  plot.args <- list(x = x$recall, y = x$precision, type = "l",
                    las = 1, lwd = 2,
                    xlab = "Recall", ylab = "Precision")
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
}

#' Add precision-recall curve to a plot
#' 
#' Adds a precision-recall curve to a base graphics plot. It accepts an S3
#' object of class 'prc', see [prc()].
#' 
#' @param x An object of class 'prc'
#' @param ... Optional graphical arguments passed to [lines()]
#' @return No return value
#' @seealso [prc()] [plot.prc()]
#' @export
lines.prc <- function(x, ...) {
  new.args <- list(...)
  plot.args <- list(x = x$recall, y = x$precision, lwd = 2)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("lines", plot.args)
}

