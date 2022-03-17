# nestedcv plots

#' Plot cross-validated glmnet alpha
#' 
#' Plot of cross-validated glmnet alpha parameter against deviance.
#' 
#' @param x Fitted "nestcv.glmnet" object
#' @param col Optional vector of line colours for each fold
#' @param ... other arguments passed to plot
#' @return No return value
#' @seealso [nestcv.glmnet]
#' @author Myles Lewis
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @export
#' 
plot_alphas <- function(x,
                        col = NULL,
                        ...) {
  cv_alpha <- lapply(x$outer_result, function(i) i$cvafit$alpha_cvm)
  alphaSet <- x$outer_result[[1]]$cvafit$alphaSet
  n <- length(cv_alpha)
  if (is.null(col)) {
    col <- rainbow(n)
  } else {
    col <- rep_len(col, n)
  }
  new.args <- list(...)
  plot.args <- list(y = cv_alpha[[1]], x = alphaSet,
                    type = 'l',
                    ylim = range(unlist(cv_alpha)),
                    xlab = 'Alpha',
                    ylab = x$outer_result[[1]]$cvafit$fits[[1]]$name,
                    col = col[1])
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  for (i in 2:n) {
    lines.args <- list(y = cv_alpha[[i]], x = alphaSet, col = col[i])
    if (length(new.args)) lines.args[names(new.args)] <- new.args
    do.call("lines", lines.args)
  }
}


#' Plot cross-validated glmnet lambdas across outer folds
#' 
#' Plot of cross-validated glmnet lambda parameter against deviance for each 
#' outer CV fold.
#' 
#' @param x Fitted "nestcv.glmnet" object
#' @param scheme colour scheme
#' @param palette palette name (one of `hcl.pals()`) which is passed to 
#' [hcl.colors]
#' @param showLegend Either a keyword to position the legend or `NULL` to hide 
#' the legend.
#' @param ... other arguments passed to plot. Use `type = 'p'` to plot a 
#' scatter plot instead of a line plot.
#' @return No return value
#' @seealso [nestcv.glmnet]
#' @author Myles Lewis
#' @importFrom graphics lines abline par legend
#' @importFrom grDevices hcl.colors
#' @export
#' 
plot_lambdas <- function(x,
                         scheme = NULL,
                         palette = "Dark 3",
                         showLegend = if(x$outer_method == "cv") "topright" else NULL,
                         ...) {
  cvms <- lapply(x$outer_result, function(fold) fold$cvafit$fits[[fold$cvafit$which_alpha]]$cvm)
  lambdas <- lapply(x$outer_result, function(fold) fold$cvafit$fits[[fold$cvafit$which_alpha]]$lambda)
  n <- length(cvms)
  if (is.null(scheme)) scheme <- hcl.colors(n, palette)
  new.args <- list(...)
  plot.args <- list(y = cvms[[1]], x = log(lambdas[[1]]),
                    type = 'l',
                    ylim = range(unlist(cvms)),
                    xlim = range(log(unlist(lambdas))),
                    xlab = expression(Log(lambda)),
                    ylab = x$outer_result[[1]]$cvafit$fits[[1]]$name,
                    col = scheme[1])
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  for (i in 2:n) {
    lines.args <- list(y = cvms[[i]], x = log(lambdas[[i]]), col = scheme[i])
    if (length(new.args)) lines.args[names(new.args)] <- new.args
    do.call("lines", lines.args)
  }
  abline(v = log(x$final_param["lambda"]), lty = 2, col = 'grey')
  if (!is.null(showLegend)) {
    if (plot.args$type == 'p') {
      legend.pch <- if (is.null(plot.args$pch)) par("pch") else plot.args$pch
      legend(showLegend, bty = 'n',
             legend = paste("Fold", 1:n),
             col = scheme, pch = legend.pch)
    } else {
      legend.lwd <- if (is.null(plot.args$lwd)) par("lwd") else plot.args$lwd
      legend(showLegend, bty = 'n',
             legend = paste("Fold", 1:n),
             col = scheme, lty = 1, lwd = legend.lwd)
    }
  }
}


#' Plot lambda across range of alphas
#'
#' Different types of plot showing cross-validated tuning of alpha and lambda
#' from elastic net regression via [glmnet]. If `xaxis` is set to `"lambda"`,
#' log lambda is on the x axis while the tuning metric (log loss, deviance,
#' accuracy, AUC etc) is on the y axis. Multiple alpha values are shown by
#' different colours. If `xaxis` is set to `"alpha"`, alpha is on the x axis
#' with the tuning metric on y, with error bars showing metric SD. if `xaxis` is
#' set to `"nvar"` the number of non-zero coefficients is shown on x and how
#' this relates to model deviance/ accuracy on y.
#'
#' @param x Object of class 'cva.glmnet'
#' @param xaxis String specifying what is plotted on the x axis, either log
#'   lambda, alpha or the number of non-zero coefficients.
#' @param errorBar Logical whether to show error bars for the standard deviation
#'   of model deviance. Error bars are interleaved to avoid overlap.
#' @param errorWidth Width of error bars.
#' @param min.pch Plotting 'character' for the minimum point of each curve. Not
#'   shown if set to `NULL`. See [points]
#' @param scheme Colour scheme
#' @param palette Palette name (one of `hcl.pals()`) which is passed to
#'   [hcl.colors]
#' @param showLegend Either a keyword to position the legend or `NULL` to hide
#'   the legend.
#' @param ... Other arguments passed to [plot]. Use `type = 'p'` to plot a
#'   scatter plot instead of a line plot.
#' @return No return value
#' @seealso [nestcv.glmnet]
#' @author Myles Lewis
#' @importFrom graphics lines legend par
#' @importFrom grDevices hcl.colors
#' @export
#' 
plot.cva.glmnet <- function(x,
                            xaxis = c('lambda', 'alpha', 'nvar'),
                            errorBar = (xaxis == "alpha"),
                            errorWidth = 0.01,
                            min.pch = NULL,
                            scheme = NULL,
                            palette = "zissou",
                            showLegend = "bottomright",
                            ...) {
  xaxis <- match.arg(xaxis)
  new.args <- list(...)
  if (xaxis == "alpha") {
    px <- x$alphaSet
    y <- unlist(lapply(x$fits, function(i) min(i$cvm)))
    ylo <- unlist(lapply(x$fits, function(i) {w <- which.min(i$cvm)
      i$cvlo[w]}))
    yup <- unlist(lapply(x$fits, function(i) {w <- which.min(i$cvm)
      i$cvup[w]}))
    ylim <- range(c(y, ylo, yup))
    plot.args <- list(y = y, x = px,
                      ylim = ylim,
                      xlab = "Alpha",
                      ylab = x$fits[[1]]$name,
                      pch = 21,
                      bg = "white")
    if (length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    if (errorBar) arrows(px, ylo, px, yup, length = errorWidth, angle = 90, code = 3)
    do.call("points", plot.args)
    return(invisible())
  }
  cvms <- lapply(x$fits, function(i) i$cvm)
  n <- length(cvms)
  if (is.null(scheme)) scheme <- hcl.colors(n, palette)
  px <- switch(xaxis,
               lambda = lapply(x$fits, function(i) log(i$lambda)),
               nvar = lapply(x$fits, function(i) i$nzero))
  ylim <- range(unlist(cvms))
  if (errorBar) {
    cvlo <- lapply(x$fits, function(i) i$cvlo)
    cvup <- lapply(x$fits, function(i) i$cvup)
    ylim <- range(c(unlist(cvlo), unlist(cvup)))
  }
  type <- switch(xaxis,
                lambda = 'l',
                nvar = 'p')
  cex <- switch(xaxis,
               lambda = 0.5,
               nvar = 0.8)
  plot.args <- list(y = cvms[[1]], x = px[[1]],
                    type = type,
                    ylim = ylim,
                    xlim = range(unlist(px)),
                    xlab = switch(xaxis,
                                  lambda = expression(Log(lambda)),
                                  nvar = "Number of non-zero coefficients"),
                    ylab = x$fits[[1]]$name,
                    cex = cex,
                    col = scheme[1])
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  if (errorBar) {
    ind <- seq(1, length(cvms[[1]]), 10)
    for (i in 1:n) {
      arrows(px[[i]][ind], cvlo[[i]][ind], px[[i]][ind], cvup[[i]][ind],
             length = errorWidth, angle = 90, code = 3, col = scheme[i])
    }
  }
  for (i in 2:n) {
    lines.args <- list(y = cvms[[i]], x = px[[i]], col = scheme[i],
                       type = type,
                       cex = cex)
    if (length(new.args)) lines.args[names(new.args)] <- new.args
    do.call("lines", lines.args)
  }
  if (errorBar) {
    for (i in 1:n) {
      points(px[[i]][ind], cvms[[i]][ind], pch = 19, col = scheme[i],
             cex = lines.args$cex)
    }
  }
  # show minima
  if (!is.null(min.pch)) {
    wy <- unlist(lapply(cvms, which.min))
    for (i in 1:n) {
      points(x = px[[i]][wy[i]], y = cvms[[i]][wy[i]], col = scheme[i], pch = min.pch)
    }
  }
  if (!is.null(showLegend)) {
    if (plot.args$type == 'p') {
      legend.pch <- if (is.null(plot.args$pch)) par("pch") else plot.args$pch
      legend(showLegend, bty = 'n',
             legend = parse(text = paste("alpha ==", x$alphaSet)),
             col = scheme, pch = legend.pch)
    } else {
      legend.lwd <- if (is.null(plot.args$lwd)) par("lwd") else plot.args$lwd
      legend(showLegend, bty = 'n',
             legend = parse(text = paste("alpha ==", x$alphaSet)),
             col = scheme, lty = 1, lwd = legend.lwd)
    }
  }
}


#' Boxplot model predictors
#' 
#' Boxplots to show range of model predictors to identify exceptional predictors
#' with excessively low or high values.
#' 
#' @param fit "nestedcv" object
#' @param x matrix of predictors
#' @param scheme colour scheme
#' @param palette palette name (one of `hcl.pals()`) which is passed to 
#' [hcl.colors]
#' @param ... other arguments passed to [boxplot].
#' @return No return value
#' @seealso [nestcv.glmnet]
#' @author Myles Lewis
#' @importFrom graphics boxplot par
#' @importFrom grDevices hcl.colors
#' @importFrom stats formula
#' @export
#' 
boxplot_model <- function(fit, x,
                          scheme = NULL, palette = "Dark 3", ...) {
  m <- coef(fit)
  m <- names(m[-1])
  df <- data.frame(vars = rep(m, each = nrow(x)), 
                   y = unlist(lapply(m, function(i) x[, i])))
  x_med <- Rfast::colMedians(x[, m])
  names(x_med) <- m
  x_med <- sort(x_med, decreasing = TRUE)
  df$vars <- factor(df$vars, levels = names(x_med))
  if (is.null(scheme)) scheme <- hcl.colors(length(m), palette)
  new.args <- list(...)
  plot.args <- list(formula = formula(df$y ~ df$vars),
                    xlab = "", ylab = "",
                    boxcol = scheme, medcol = scheme,
                    whiskcol = scheme, staplecol = scheme, outcol = scheme,
                    col = NA,
                    las = 2, outpch = 20, outcex = 0.5,
                    whisklty = 1, 
                    cex.lab = 0.75, cex.axis = 0.75)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("boxplot", plot.args)
}

#' Plot caret tuning
#' 
#' Plots the main tuning parameter in models built using [caret::train]
#' 
#' @param x Object of class 'train' generated by `caret` function [train]
#' @param error.col Colour of error bars
#' @param ... Other arguments passed to [plot()]
#' @return No return value
#' @importFrom graphics arrows points
#' @export
#' 
plot_caret <- function(x, error.col = "darkgrey", ...) {
  res <- x$results
  nrepeat <- x$control$number * x$control$repeats
  w <- switch(x$method, glmnet = "lambda", 1)
  x1 <- res[, w]
  y <- res[, x$metric]
  sem <- res[, paste0(x$metric, "SD")] / sqrt(nrepeat)
  xlab <- colnames(res)[w]
  if (x$method == "glmnet") {
    x1 <- log(x1)
    xlab <- "log lambda"
  }
  new.args = list(...)
  plot.args = list(x = x1, y = y,
                   xlab = xlab,
                   ylab = x$metric,
                   ylim = range(c(y - sem, y + sem)))
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  arrows(x1, y - sem, x1, y + sem,
         length = 0.025, angle = 90, code = 3, col = error.col)
  plot.args = list(x = x1, y = y,
                   pch = 21, bg = "white")
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("points", plot.args)
}


