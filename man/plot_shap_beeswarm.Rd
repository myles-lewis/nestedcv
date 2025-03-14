% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shap.R
\name{plot_shap_beeswarm}
\alias{plot_shap_beeswarm}
\title{SHAP importance beeswarm plot}
\usage{
plot_shap_beeswarm(
  shap,
  x,
  cex = 0.25,
  corral = "random",
  corral.width = 0.7,
  scheme = c("deepskyblue2", "purple3", "red"),
  sort = TRUE,
  top = NULL,
  ...
)
}
\arguments{
\item{shap}{a matrix of SHAP values}

\item{x}{a matrix or dataframe of feature values containing only features
values from the training data. The rows must match rows in \code{shap}. If a
dataframe is supplied it is converted to a numeric matrix using
\code{\link[=data.matrix]{data.matrix()}}.}

\item{cex}{Scaling for adjusting point spacing. See
\code{ggbeeswarm::geom_beeswarm()}.}

\item{corral}{String specifying method used to corral points. See
\code{ggbeeswarm::geom_beeswarm()}.}

\item{corral.width}{Numeric specifying width of corral, passed to
\code{geom_beeswarm}}

\item{scheme}{Colour scheme as a vector of 3 colours}

\item{sort}{Logical whether to sort predictors by mean absolute SHAP value.}

\item{top}{Sets a limit on the number of variables plotted or \code{NULL} to plot
all variables. If \code{top} is set then variables are sorted and \code{sort} is
overrode.}

\item{...}{Other arguments passed to \code{ggbeeswarm::geom_beeswarm()} e.g.
\code{size}.}
}
\value{
A ggplot2 plot
}
\description{
SHAP importance beeswarm plot
}
