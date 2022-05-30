# Train and predict wrapper functions for calling hsstan with filters from the
# nestedcv package. Only outer.cv is needed

#' hsstan model to be fitted with cross-validation
#'
#' This function applies a cross-validation (CV) procedure for training Bayesian
#' models with hierarchical shrinkage priors using the `hsstan` package.
#' The function allows the option of embedded filtering of predictors for
#' feature selection within the CV loop. Within each training fold, an optional
#' filtering of predictors is performed, followed by fitting of an `hsstsan`
#' model. Predictions on the testing folds are brought back together and error
#' estimation/ accuracy determined. The default is 10-fold CV.
#' The function is implemented within the `nestedcv` package. The `hsstan`
#' models do not require tuning of meta-parameters and therefore only a single
#' CV procedure is needed to evaluate performance. This is implemented using the
#' `outer` CV procedure in the `nestedcv` package.
#'
#' @param y Response vector. For classification this should be a factor.
#' @param x Matrix of predictors
#' @param unpenalized Vector of column names `x` which are always retained into
#'   the model (i.e. not penalized). Default `NULL` means the parameters for all
#'   predictors will be drawn from a hierarchical prior distribution, i.e. will
#'   be penalized. Note: if filtering of predictors is specified, then the
#'   vector of `unpenalized` predictors should also be passed to the filter
#'   function using the `filter_options$force_vars` argument. Filters currently
#'   implementing this option are the `partial_ttest_filter` for binary outcomes
#'   and the `lm_filter` for continuous outcomes.
#'
#' @return An object of class `hsstan`
#'
#' @author Athina Spiliopoulou
#' @importFrom hsstan hsstan
#' @importFrom data.table as.data.table
#' @export
#'
model.hsstan <- function(y, x, unpenalized = NULL, ...) {

    ## reformat outcome and predictors to work with hsstan
    if(nlevels(y) == 2) {
        ## for binary outcomes convert factor to numeric (0, 1)
        ## need to keep names of factor levels for caret confusion matrix (the
        ## factor levels are applied back to the outcome when we use predict)
        y.levels <- levels(y)
        y <- as.integer(y) - 1
        family = "binomial"
    } else {
        family = "gaussian"
    }

    dt <- as.data.table(cbind(y, x))

    ## handle unpenalized predictors
    if(!is.null(unpenalized)) {
        covs.model <- paste(c("y ~ 1", unpenalized), collapse = " + ")
        penalized <- setdiff(colnames(x), unpenalized)
    } else {
        covs.model <- "y ~ 1"
        penalized <- colnames(x)
    }

    ## fit the model
    fit <- hsstan(x = dt, covs.model = covs.model, penalized = penalized,
                  family = family, ...)

    if(exists("y.levels")) {
        fit$y.levels <- y.levels
    }

    return(fit)
}

#' predict from hsstan model fitted within cross-validation
#'
#' Draws from the posterior predictive distribution of the outcome.
#'
#' @param object: An object of class `hsstan`.
#' @param newdata: Optional data frame containing the variables to use to
#'   predict. If `NULL` (default), the model matrix is used. If specified, its
#'   continuous variables should be standardized, since the model coefficients
#'   are learnt on standardized data.
#' @param type: Option for binary outcomes only. Default `NULL` will return a
#'   class with the highest probability for each sample. If set to `probs`, it
#'   will return the probabilities for outcome = 0 and for outcome = 1 for each
#'   sample.
#'
#' @return For a binary outcome and type = `NULL`, a character vector with the
#'   name of the class that has the highest probability for each sample.
#'   For a binary outcome and type = `prob`, a 2-dimensional matrix with the
#'   probability of class 0 and of class 1 for each sample.
#'   For a continuous outcome a numeric vector with the predicted value for
#'   each sample.
#'
#' @author Athina Spiliopoulou
#' @importFrom hsstan posterior_predict
#' @export
#'
predict.hsstan <- function(object, newdata = NULL, type = NULL, ...) {
    if (object$family[1] == "gaussian") {
        preds <- colMeans(posterior_predict(object, newdata = newdata, ...))
    } else if (object$family[1] == "binomial") {
        probs <- colMeans(posterior_predict(object, newdata = newdata,
                                            transform = TRUE, ...))
        if(is.null(type)) {
            # convert to binary
            preds <- ifelse(probs < 0.5, object$y.levels[1], object$y.levels[2])
        } else if (type == "prob") {
            # return probabilities for outcome = 0 and outcome = 1
            preds <- cbind(1 - probs, probs)
            colnames(preds) <- c("prob_class0", "prob_class1")
        }
    }
    return(preds)
}
