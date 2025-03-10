% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nestedcv.R
\name{nestcv.glmnet}
\alias{nestcv.glmnet}
\title{Nested cross-validation with glmnet}
\usage{
nestcv.glmnet(
  y,
  x,
  family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
  filterFUN = NULL,
  filter_options = NULL,
  balance = NULL,
  balance_options = NULL,
  modifyX = NULL,
  modifyX_useY = FALSE,
  modifyX_options = NULL,
  outer_method = c("cv", "LOOCV"),
  n_outer_folds = 10,
  n_inner_folds = 10,
  outer_folds = NULL,
  pass_outer_folds = FALSE,
  alphaSet = seq(0.1, 1, 0.1),
  min_1se = 0,
  keep = TRUE,
  outer_train_predict = FALSE,
  weights = NULL,
  penalty.factor = rep(1, ncol(x)),
  parallel_mode = NULL,
  cv.cores = 1,
  finalCV = TRUE,
  na.option = "omit",
  verbose = FALSE,
  ...
)
}
\arguments{
\item{y}{Response vector or matrix. Matrix is only used for
\code{family = 'mgaussian'} or \code{'cox'}.}

\item{x}{Matrix of predictors. Dataframes will be coerced to a matrix as
is necessary for glmnet.}

\item{family}{Either a character string representing one of the built-in
families, or else a \code{glm()} family object. Passed to \link[glmnet:cv.glmnet]{glmnet::cv.glmnet}
and \link[glmnet:glmnet]{glmnet::glmnet}}

\item{filterFUN}{Filter function, e.g. \link{ttest_filter} or \link{relieff_filter}.
Any function can be provided and is passed \code{y} and \code{x}. Must return a
numeric vector with indices of filtered predictors.}

\item{filter_options}{List of additional arguments passed to the filter
function specified by \code{filterFUN}.}

\item{balance}{Specifies method for dealing with imbalanced class data.
Current options are \code{"randomsample"} or \code{"smote"}. See \code{\link[=randomsample]{randomsample()}} and
\code{\link[=smote]{smote()}}}

\item{balance_options}{List of additional arguments passed to the balancing
function}

\item{modifyX}{Character string specifying the name of a function to modify
\code{x}. This can be an imputation function for replacing missing values, or a
more complex function which alters or even adds columns to \code{x}. The
required return value of this function depends on the \code{modifyX_useY}
setting.}

\item{modifyX_useY}{Logical value whether the \code{x} modifying function makes
use of response training data from \code{y}. If \code{FALSE} then the \code{modifyX}
function simply needs to return a modified \code{x} object, which will be
coerced to a matrix as required by \code{glmnet}. If \code{TRUE} then the \code{modifyX}
function must return a model type object on which \code{predict()} can be
called, so that train and test partitions of \code{x} can be modified
independently.}

\item{modifyX_options}{List of additional arguments passed to the \code{x}
modifying function}

\item{outer_method}{String of either \code{"cv"} or \code{"LOOCV"} specifying whether
to do k-fold CV or leave one out CV (LOOCV) for the outer folds}

\item{n_outer_folds}{Number of outer CV folds}

\item{n_inner_folds}{Number of inner CV folds}

\item{outer_folds}{Optional list containing indices of test folds for outer
CV. If supplied, \code{n_outer_folds} is ignored.}

\item{pass_outer_folds}{Logical indicating whether the same outer folds are
used for fitting of the final model when final CV is applied. Note this can
only be applied when \code{n_outer_folds} and \code{n_inner_folds} are the same and
no balancing is applied.}

\item{alphaSet}{Vector of alphas to be tuned}

\item{min_1se}{Value from 0 to 1 specifying choice of optimal lambda from
0=lambda.min to 1=lambda.1se}

\item{keep}{Logical indicating whether inner CV predictions are retained for
calculating left-out inner CV fold accuracy etc. See argument \code{keep} in
\link[glmnet:cv.glmnet]{glmnet::cv.glmnet}.}

\item{outer_train_predict}{Logical whether to save predictions on outer
training folds to calculate performance on outer training folds.}

\item{weights}{Weights applied to each sample. Note \code{weights} and \code{balance}
cannot be used at the same time. Weights are only applied in glmnet and not
in filters.}

\item{penalty.factor}{Separate penalty factors can be applied to each
coefficient. Can be 0 for some variables, which implies no shrinkage, and
that variable is always included in the model. Default is 1 for all
variables. See \link[glmnet:glmnet]{glmnet::glmnet}. Note this works separately from filtering.
For some \code{nestedcv} filter functions you might need to set \code{force_vars} to
avoid filtering out features.}

\item{parallel_mode}{Either "mclapply", "parLapply" or "future". This
determines which parallel backend to use. The default is
\code{parallel::mclapply} on unix/mac and \code{parallel::parLapply} on windows.}

\item{cv.cores}{Number of cores for parallel processing of the outer loops.
Ignored if \code{parallel_mode = "future"}.}

\item{finalCV}{Logical whether to perform one last round of CV on the whole
dataset to determine the final model parameters. If set to \code{FALSE}, the
median of hyperparameters from outer CV folds are used for the final model.
Performance metrics are independent of this last step. If set to \code{NA},
final model fitting is skipped altogether, which gives a useful speed boost
if performance metrics are all that is needed.}

\item{na.option}{Character value specifying how \code{NA}s are dealt with.
\code{"omit"} (the default) is equivalent to \code{na.action = na.omit}. \code{"omitcol"}
removes cases if there are \code{NA} in 'y', but columns (predictors) containing
\code{NA} are removed from 'x' to preserve cases. Any other value means that
\code{NA} are ignored (a message is given).}

\item{verbose}{Logical whether to print messages and show progress}

\item{...}{Optional arguments passed to \link[glmnet:cv.glmnet]{glmnet::cv.glmnet}}
}
\value{
An object with S3 class "nestcv.glmnet"
\item{call}{the matched call}
\item{output}{Predictions on the left-out outer folds}
\item{outer_result}{List object of results from each outer fold containing
predictions on left-out outer folds, best lambda, best alpha, fitted glmnet
coefficients, list object of inner fitted cv.glmnet and number of filtered
predictors at each fold.}
\item{outer_method}{the \code{outer_method} argument}
\item{n_inner_folds}{number of inner folds}
\item{outer_folds}{List of indices of outer test folds}
\item{dimx}{dimensions of \code{x}}
\item{xsub}{subset of \code{x} containing all predictors used in both outer CV
folds and the final model}
\item{y}{original response vector}
\item{yfinal}{final response vector (post-balancing)}
\item{final_param}{Final mean best lambda
and alpha from each fold}
\item{final_fit}{Final fitted glmnet model}
\item{final_coef}{Final model coefficients and mean expression. Variables
with coefficients shrunk to 0 are removed.}
\item{final_vars}{Column names of filtered predictors entering final model.
This is useful for subsetting new data for predictions.}
\item{roc}{ROC AUC for binary classification where available.}
\item{summary}{Overall performance summary. Accuracy and balanced accuracy
for classification. ROC AUC for binary classification. RMSE for
regression.}
}
\description{
This function enables nested cross-validation (CV) with glmnet including
tuning of elastic net alpha parameter. The function also allows the option of
embedded filtering of predictors for feature selection nested within the
outer loop of CV. Predictions on the outer test folds are brought back
together and error estimation/ accuracy determined. The default is 10x10
nested CV.
}
\details{
glmnet does not tolerate missing values, so \code{na.option = "omit"} is the
default.
}
\examples{
\donttest{
## Example binary classification problem with P >> n
x <- matrix(rnorm(150 * 2e+04), 150, 2e+04)  # predictors
y <- factor(rbinom(150, 1, 0.5))  # binary response

## Partition data into 2/3 training set, 1/3 test set
trainSet <- caret::createDataPartition(y, p = 0.66, list = FALSE)

## t-test filter using whole dataset
filt <- ttest_filter(y, x, nfilter = 100)
filx <- x[, filt]

## Train glmnet on training set only using filtered predictor matrix
library(glmnet)
fit <- cv.glmnet(filx[trainSet, ], y[trainSet], family = "binomial")
plot(fit)

## Predict response on test partition
predy <- predict(fit, newx = filx[-trainSet, ], s = "lambda.min", type = "class")
predy <- as.vector(predy)
predyp <- predict(fit, newx = filx[-trainSet, ], s = "lambda.min", type = "response")
predyp <- as.vector(predyp)
output <- data.frame(testy = y[-trainSet], predy = predy, predyp = predyp)

## Results on test partition
## shows bias since univariate filtering was applied to whole dataset
predSummary(output)

## Nested CV
## n_outer_folds reduced to speed up example
fit2 <- nestcv.glmnet(y, x, family = "binomial", alphaSet = 1,
                      n_outer_folds = 3,
                      filterFUN = ttest_filter,
                      filter_options = list(nfilter = 100),
                      cv.cores = 2)
summary(fit2)
plot_lambdas(fit2, showLegend = "bottomright")

## ROC plots
library(pROC)
testroc <- roc(output$testy, output$predyp, direction = "<")
inroc <- innercv_roc(fit2)
plot(fit2$roc)
lines(inroc, col = 'blue')
lines(testroc, col = 'red')
legend('bottomright', legend = c("Nested CV", "Left-out inner CV folds", 
                                 "Test partition, non-nested filtering"), 
       col = c("black", "blue", "red"), lty = 1, lwd = 2, bty = "n")
}
}
\author{
Myles Lewis
}
