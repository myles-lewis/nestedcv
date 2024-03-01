
#' Model performance metrics
#' 
#' Returns model metrics from nestedcv models. Extended metrics including 
#' 
#' @param object A 'nestcv.glmnet', 'nestcv.train', 'nestcv.SuperLearner' or
#'   'outercv' object.
#' @param extra Logical whether additional performance metrics are gathered for
#'   binary classification models: area under precision recall curve (AUC.PR),
#'   Cohen's kappa, F1 score, Matthew's correlation coefficient (MCC).
#' @param innerCV Whether to calculate metrics for inner CV folds. Only
#'   available for 'nestcv.glmnet' and 'nestcv.train' objects.
#' @param positive For binary classification, either an integer 1 or 2 for the
#'   level of response factor considered to be 'positive' or 'relevant', or a
#'   character value for that factor. This affects the F1 score. See
#'   [caret::confusionMatrix()].
#' @returns A named numeric vector of performance metrics.
#' @export
#'
metrics <- function(object, extra = FALSE, innerCV = FALSE, positive = 2) {
  if (!is.list(object$summary)) {
    # regression
    met <- object$summary
    if (innerCV && inherits(object, c("nestcv.glmnet", "nestcv.train"))) {
      inner_met <- innercv_summary(object)
      names(inner_met) <- paste0("innerCV.", names(inner_met))
      met <- c(met, inner_met)
    }
    return(c(met, nvar(object)))
  }
  met <- object$summary$metrics
  if (extra && nlevels(object$y) == 2) {
    # binary classification
    if (!requireNamespace("MLmetrics", quietly = TRUE)) {
      stop("Package 'MLmetrics' must be installed", call. = FALSE)
    }
    aucpr <- MLmetrics::PRAUC(object$output$predyp, object$output$testy)
    tab <- object$summary$table
    mcc <- mcc(tab)
    if (is.numeric(positive)) positive <- colnames(tab)[positive]
    ccm <- caret::confusionMatrix(tab, mode = "everything", positive = positive)
    extra <- setNames(c(aucpr, ccm$overall["Kappa"], ccm$byClass["F1"], mcc), 
                      c("AUC.PR", "Kappa", paste("F1", positive, sep = "."),
                        "MCC"))
    met <- c(met, extra)
  }
  if (innerCV && inherits(object, c("nestcv.glmnet", "nestcv.train"))) {
    inner_met <- innercv_summary(object)$metrics
    names(inner_met) <- paste0("innerCV.", names(inner_met))
    met <- c(met, inner_met)
  }
  c(met, nvar(object))
}


# obtain number of vars in models
nvar <- function(object) {
  nvar <- if (inherits(object, "nestcv.glmnet")) {
    if (identical(object$final_coef, NA)) {
      # from outer CV results
      ncf <- unlist(lapply(object$outer_result, function(i) length(i$coef)))
      round(median(ncf) -1)
    } else {
      if (!is.data.frame(object$final_coef)) {
        # multinomial glmnet
        cf <- coef(object$final_fit)
        cflist <- lapply(cf, as.matrix)
        cf <- do.call(cbind, cflist)
        ok <- rowSums(cf != 0)
        sum(ok != 0) -1
      } else nrow(object$final_coef) -1
    }
  } else if (inherits(object, c("nestcv.train", "nestcv.SuperLearner")) &&
             (is.null(object$final_vars) | identical(object$final_vars, NA))) {
    # from outer CV results
    nfilter <- unlist(lapply(object$outer_result, '[[', 'nfilter'))
    round(median(nfilter))
  } else length(object$final_vars)
  setNames(nvar, "nvar")
}


# Matthew's correlation coefficient
mcc <- function(cm) {
  tp <- cm[2, 2]
  tn <- cm[1, 1]
  fp <- cm[2, 1]
  fn <- cm[1, 2]
  (tp*tn - fp*fn) /
    sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn))
}
