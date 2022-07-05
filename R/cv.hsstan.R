# Train and predict wrapper functions for calling hsstan with filters from the
# nestedcv package. Only outer.cv is needed

#' hsstan model for cross-validation
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
#' @param ... Optional arguments passed to `hsstan`
#'
#' @return An object of class `hsstan`
#'
#' @author Athina Spiliopoulou
#' @importFrom hsstan hsstan
#' @importFrom data.table as.data.table
#' @examples
#' \donttest{
#' # Cross-validation is used to apply univariate filtering of predictors.
#' # only one CV split is needed (outercv) as the Bayesian model does not 
#' # require learning of meta-parameters.
#' 
#' # load iris dataset and simulate a continuous outcome
#' data(iris)
#' dt <- iris[, 1:4]
#' colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
#' dt <- as.data.frame(apply(dt, 2, scale))
#' dt$outcome.cont <- -3 + 0.5 * dt$marker1 + 2 * dt$marker2 + rnorm(nrow(dt), 0, 2)
#' 
#' # unpenalised covariates: always retain in the prediction model
#' uvars <- "marker1"
#' # penalised covariates: coefficients are drawn from hierarchical shrinkage
#' # prior
#' pvars <- c("marker2", "marker3", "marker4") # penalised covariates
#' # run cross-validation with univariate filter and hsstan
#' # dummy sampling for fast execution of example
#' # recommend 4 chains, warmup 1000, iter 2000 in practice
#' oldopt <- options(mc.cores = 2)
#' res.cv.hsstan <- outercv(y = dt$outcome.cont, x = dt[, c(uvars, pvars)],
#'                          model = model.hsstan,
#'                          filterFUN = lm_filter,
#'                          filter_options = list(force_vars = uvars,
#'                                                nfilter = 2,
#'                                                p_cutoff = NULL,
#'                                                rsq_cutoff = 0.9),
#'                          n_outer_folds = 3, chains = 2,
#'                          unpenalized = uvars, warmup = 100, iter = 200)
#' # view prediction performance based on testing folds
#' res.cv.hsstan$summary
#' # view coefficients for the final model
#' res.cv.hsstan$final_fit
#' # view covariates selected by the univariate filter
#' res.cv.hsstan$final_vars
#' 
#' # load hsstan package to examine the Bayesian model
#' library(hsstan)
#' sampler.stats(res.cv.hsstan$final_fit)
#' print(projsel(res.cv.hsstan$final_fit), digits = 4) # adding marker2
#' options(oldopt)
#' 
#' # Here adding `marker2` improves the model fit: substantial decrease of
#' # KL-divergence from the full model to the submodel. Adding `marker3` does 
#' # not improve the model fit: no decrease of KL-divergence from the full model 
#' # to the submodel.
#' }
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

#' Predict from hsstan model fitted within cross-validation
#'
#' Draws from the posterior predictive distribution of the outcome.
#'
#' @param object An object of class `hsstan`.
#' @param newdata Optional data frame containing the variables to use to
#'   predict. If `NULL` (default), the model matrix is used. If specified, its
#'   continuous variables should be standardized, since the model coefficients
#'   are learnt on standardized data.
#' @param type Option for binary outcomes only. Default `NULL` will return a
#'   class with the highest probability for each sample. If set to `probs`, it
#'   will return the probabilities for outcome = 0 and for outcome = 1 for each
#'   sample.
#' @param ... Optional arguments passed to `hsstan::posterior_predict`
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
