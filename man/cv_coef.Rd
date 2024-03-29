% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varImp.R
\name{cv_coef}
\alias{cv_coef}
\title{Coefficients from outer CV glmnet models}
\usage{
cv_coef(x, level = 1)
}
\arguments{
\item{x}{a \code{nestcv.glmnet} fitted object}

\item{level}{For multinomial models only, either an integer specifying which
level of outcome is being examined, or the level can be specified as a
character value}
}
\value{
matrix of coefficients from outer CV glmnet models plus the final
glmnet model. Coefficients for variables which are not present in a
particular outer CV fold model are set to 0.
}
\description{
Extracts coefficients from outer CV glmnet models from a \code{nestcv.glmnet}
fitted object.
}
\seealso{
\code{\link[=cv_varImp]{cv_varImp()}}
}
