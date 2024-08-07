% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{metrics}
\alias{metrics}
\title{Model performance metrics}
\usage{
metrics(object, extra = FALSE, innerCV = FALSE, positive = 2)
}
\arguments{
\item{object}{A 'nestcv.glmnet', 'nestcv.train', 'nestcv.SuperLearner' or
'outercv' object.}

\item{extra}{Logical whether additional performance metrics are gathered for
classification models: area under precision recall curve (PR.AUC, binary
classification only), Cohen's kappa, F1 score, Matthews correlation
coefficient (MCC).}

\item{innerCV}{Whether to calculate metrics for inner CV folds. Only
available for 'nestcv.glmnet' and 'nestcv.train' objects.}

\item{positive}{For binary classification, either an integer 1 or 2 for the
level of response factor considered to be 'positive' or 'relevant', or a
character value for that factor. This affects the F1 score. See
\code{\link[caret:confusionMatrix]{caret::confusionMatrix()}}.}
}
\value{
A named numeric vector of performance metrics.
}
\description{
Returns model metrics from nestedcv models. Extended metrics including
}
\details{
Area under precision recall curve is estimated by trapezoidal estimation
using \code{MLmetrics::PRAUC()}.

For multi-class classification models, Matthews correlation coefficient is
calculated using Gorodkin's method. Multi-class F1 score (macro F1) is
calculated as the arithmetic mean of the class-wise F1 scores.
}
\references{
Gorodkin, J. (2004). \emph{Comparing two K-category assignments by a
K-category correlation coefficient}. Computational Biology and Chemistry. 28
(5): 367–374.
}
\seealso{
\code{\link[=mcc]{mcc()}}
}
