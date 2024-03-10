
#' Build precision-recall curve
#' 
#' Builds a precision-recall curve for a 'nestedcv' model using `prediction()`
#' and `performance()` functions from the ROCR package and returns an object of
#' class 'prc' for plotting.
#' 
#' @param response binary factor vector of response of default order controls,
#'   cases.
#' @param predictor numeric vector of probabilities
#' @param positive Either an integer 1 or 2 for the level of response factor
#'   considered to be 'positive' or 'relevant', or a character value for that
#'   factor.
#' @param output data.frame with columns `testy` containing observed response
#'   from test folds, and `predyp` predicted probabilities for classification
#' @param object a 'nestcv.glmnet', 'nestcv.train', 'nestcv.SuperLearn',
#'   'outercv' or 'repeatcv' S3 class results object.
#' @param ... other arguments
#' @returns An object of S3 class 'prc' containing the following fields:
#' \item{recall}{vector of recall values}
#' \item{precision}{vector of precision values}
#' \item{auc}{area under precision-recall curve value using trapezoid method}
#' \item{baseline}{baseline precision value}
#' @examples
#' \donttest{
#' library(mlbench)
#' data(Sonar)
#' y <- Sonar$Class
#' x <- Sonar[, -61]
#' 
#' fit1 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1, cv.cores = 2)
#' 
#' fit1$prc <- prc(fit1)  # calculate precision-recall curve
#' fit1$prc$auc  # precision-recall AUC value
#' 
#' fit2 <- nestcv.train(y, x, method = "gbm", cv.cores = 2)
#' fit2$prc <- prc(fit2)
#' fit2$prc$auc
#' 
#' plot(fit1$prc, ylim = c(0, 1))
#' lines(fit2$prc, col = "red")
#' 
#' ## use magritte pipe
#' `%|>%` <- magrittr::pipe_nested
#' res <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1) %|>%
#'   repeatcv(n = 4, rep.cores = 2)
#' 
#' res$prc <- prc(res)  # precision-recall curve on repeated predictions
#' plot(res$prc)
#' }
#' @export
prc <- function(...) {
  UseMethod("prc")
}

#' @rdname prc
#' @export
prc.default <- function(response, predictor, positive = 2, ...) {
  if (positive == 1 | positive == levels(response)[1]) {
    pred_obj <- ROCR::prediction(predictor, response,
                                 label.ordering = rev(levels(response)))
  } else pred_obj <- ROCR::prediction(predictor, response)
  perf_obj <- ROCR::performance(pred_obj, measure = "prec", x.measure = "rec")
  bl <- pred_obj@n.pos[[1]] / length(response)
  x <- perf_obj@x.values[[1]]
  y <- perf_obj@y.values[[1]]
  auc <- auc_calc(x, y)
  out <- list(recall = x, precision = y, auc = auc, baseline = bl)
  class(out) <- "prc"
  out
}

#' @rdname prc
#' @export
prc.data.frame <- function(output, ...) {
  if (!all(c("testy", "predyp") %in% colnames(output)))
    stop("not binary classification")
  prc.default(output$testy, output$predyp, ...)
}

#' @rdname prc
#' @export
prc.nestcv.glmnet <- function(object, ...) {
  prc.data.frame(object$output, ...)
}

#' @rdname prc
#' @export
prc.nestcv.train <- function(object, ...) {
  prc.data.frame(object$output, ...)
}

#' @rdname prc
#' @export
prc.nestcv.SuperLearner <- function(object, ...) {
  prc.data.frame(object$output, ...)
}

#' @rdname prc
#' @export
prc.outercv <- function(object, ...) {
  prc.data.frame(object$output, ...)
}

#' @rdname prc
#' @export
prc.repeatcv <- function(object, ...) {
  prc.data.frame(object$output, ...)
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
#' @examples
#' \donttest{
#' library(mlbench)
#' data(Sonar)
#' y <- Sonar$Class
#' x <- Sonar[, -61]
#' 
#' fit1 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1, cv.cores = 2)
#' fit1$prc <- prc(fit1)  # calculate precision-recall curve
#' 
#' fit2 <- nestcv.train(y, x, method = "gbm", cv.cores = 2)
#' fit2$prc <- prc(fit2)
#' 
#' plot(fit1$prc)
#' lines(fit2$prc, col = "red")
#' }
#' @export
plot.prc <- function(x, ...) {
  new.args <- list(...)
  plot.args <- list(x = x$recall, y = x$precision, type = "l",
                    las = 1, lwd = 2,
                    xlab = "Recall", ylab = "Precision")
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  abline(h = x$baseline, col = "grey", lty = 2)
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

