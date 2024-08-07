% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cva.glmnet.R
\name{coef.nestcv.glmnet}
\alias{coef.nestcv.glmnet}
\title{Extract coefficients from nestcv.glmnet object}
\usage{
\method{coef}{nestcv.glmnet}(object, s = object$final_param["lambda"], ...)
}
\arguments{
\item{object}{Object of class \code{"nestcv.glmnet"}}

\item{s}{Value of penalty parameter lambda. Default is the mean of lambda
values selected across each outer fold.}

\item{...}{Other arguments passed to \code{\link[glmnet:predict.glmnet]{glmnet::coef.glmnet()}}}
}
\value{
Vector or list of coefficients ordered with the intercept first,
followed by highest absolute value to lowest.
}
\description{
Extracts coefficients from the final fit of a \code{"nestcv.glmnet"} object.
}
