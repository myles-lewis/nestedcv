% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shap.R
\name{pred_nestcv_glmnet}
\alias{pred_nestcv_glmnet}
\alias{pred_nestcv_glmnet_class1}
\alias{pred_nestcv_glmnet_class2}
\alias{pred_nestcv_glmnet_class3}
\alias{pred_train}
\alias{pred_train_class1}
\alias{pred_train_class2}
\alias{pred_train_class3}
\alias{pred_SuperLearner}
\title{Prediction wrappers to use fastshap with nestedcv}
\usage{
pred_nestcv_glmnet(x, newdata)

pred_nestcv_glmnet_class1(x, newdata)

pred_nestcv_glmnet_class2(x, newdata)

pred_nestcv_glmnet_class3(x, newdata)

pred_train(x, newdata)

pred_train_class1(x, newdata)

pred_train_class2(x, newdata)

pred_train_class3(x, newdata)

pred_SuperLearner(x, newdata)
}
\arguments{
\item{x}{a \code{nestcv.glmnet} or \code{nestcv.train} object}

\item{newdata}{a matrix of new data}
}
\value{
prediction wrapper function designed for use with
\code{\link[fastshap:explain]{fastshap::explain()}}
}
\description{
Prediction wrapper functions to enable the use of the \code{fastshap} package for
generating SHAP values from \code{nestedcv} trained models.
}
\details{
These prediction wrapper functions are designed to be used with the
\code{fastshap} package. The functions \code{pred_nestcv_glmnet} and \code{pred_train} work
for \code{nestcv.glmnet} and \code{nestcv.train} models respectively for either binary
classification or regression.

For multiclass classification use \code{pred_nestcv_glmnet_class1}, \code{2} and \code{3}
for the first 3 classes. Similarly \code{pred_train_class1} etc for \link{nestcv.train}
objects. These functions can be inspected and easily modified to analyse
further classes.
}
\examples{
library(fastshap)

# Boston housing dataset
library(mlbench)
data(BostonHousing2)
dat <- BostonHousing2
y <- dat$cmedv
x <- subset(dat, select = -c(cmedv, medv, town, chas))

# Fit a glmnet model using nested CV
# Only 3 outer CV folds and 1 alpha value for speed
fit <- nestcv.glmnet(y, x, family = "gaussian", n_outer_folds = 3, alphaSet = 1)

# Generate SHAP values using fastshap::explain
# Only using 5 repeats here for speed, but recommend higher values of nsim
sh <- explain(fit, X=x, pred_wrapper = pred_nestcv_glmnet, nsim = 1)

# Plot overall variable importance
plot_shap_bar(sh, x)

# Plot beeswarm plot
plot_shap_beeswarm(sh, x, size = 1)

}
