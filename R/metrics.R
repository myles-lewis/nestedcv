
#' Model performance metrics
#' 
#' Returns model metrics from nestedcv models. Extended metrics including 
#' 
#' @param object A 'nestcv.glmnet', 'nestcv.train', 'nestcv.SuperLearner' or
#'   'outercv' object.
#' @param extra Logical whether additional performance metrics are gathered for
#'   classification models: area under precision recall curve (PR.AUC, binary
#'   classification only), Cohen's kappa, F1 score, Matthews correlation
#'   coefficient (MCC).
#' @param innerCV Whether to calculate metrics for inner CV folds. Only
#'   available for 'nestcv.glmnet' and 'nestcv.train' objects.
#' @param positive For binary classification, either an integer 1 or 2 for the
#'   level of response factor considered to be 'positive' or 'relevant', or a
#'   character value for that factor. This affects the F1 score. See
#'   [caret::confusionMatrix()].
#' @details
#' Area under precision recall curve is estimated by trapezoidal estimation 
#' using `MLmetrics::PRAUC()`.
#' 
#' For multi-class classification models, Matthews correlation coefficient is
#' calculated using Gorodkin's method. Multi-class F1 score (macro F1) is
#' calculated as the arithmetic mean of the class-wise F1 scores.
#' @returns A named numeric vector of performance metrics.
#' @references
#' Gorodkin, J. (2004). \emph{Comparing two K-category assignments by a
#' K-category correlation coefficient}. Computational Biology and Chemistry. 28
#' (5): 367–374.
#' @seealso [mcc()]
#' @export
#'
metrics <- function(object, extra = FALSE, innerCV = FALSE, positive = 2) {
  if (!inherits(object, c("nestcv.glmnet", "nestcv.train",
                          "nestcv.SuperLearner", "outercv")))
    stop("not a nestedcv model")
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
  if (inherits(object, "nestcv.glmnet") && object$call$family == "mgaussian") {
    # mgaussian
    return(object$summary)
  }
  met <- object$summary$metrics
  if (extra) {
    if (nlevels(object$y) == 2) {
      # binary classification
      aucpr <- prc(object, positive = positive)$auc
      tab <- object$summary$table
      mcc <- mcc(tab)
      if (is.numeric(positive)) positive <- colnames(tab)[positive]
      ccm <- caret::confusionMatrix(tab, mode = "everything", positive = positive)
      extra <- setNames(c(aucpr, ccm$overall["Kappa"], ccm$byClass["F1"], mcc), 
                        c("PR.AUC", "Kappa", "F1", "MCC"))
      met <- c(met, extra)
    } else if (nlevels(object$y) > 2) {
      # multiclass
      tab <- object$summary$table
      mcc <- mcc_multi(tab)
      if (is.numeric(positive)) positive <- colnames(tab)[positive]
      ccm <- caret::confusionMatrix(tab, mode = "everything", positive = positive)
      f1 <- ccm$byClass[, "F1"]
      f1.macro <- mean(f1)
      extra <- setNames(c(ccm$overall["Kappa"], f1.macro, mcc),
                        c("Kappa", "F1.macro", "MCC"))
      met <- c(met, extra)
    }
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
      round(median(ncf) -1L)
    } else {
      if (!is.data.frame(object$final_coef)) {
        # multinomial glmnet
        cf <- coef(object$final_fit)
        cflist <- lapply(cf, as.matrix)
        cf <- do.call(cbind, cflist)
        ok <- rowSums(cf != 0)
        sum(ok != 0) -1L
      } else nrow(object$final_coef) -1L
    }
  } else if (inherits(object, c("nestcv.train", "nestcv.SuperLearner")) &&
             (is.null(object$final_vars) | identical(object$final_vars, NA))) {
    # from outer CV results
    nfilter <- unlist(lapply(object$outer_result, '[[', 'nfilter'))
    round(median(nfilter))
  } else length(object$final_vars)
  setNames(nvar, "nvar")
}


#' Matthews correlation coefficient
#' 
#' Calculates Matthews correlation coefficient (MCC) which is in essence a
#' correlation coefficient between the observed and predicted binary
#' classifications. It has also been generalised to multi-class classification.
#' 
#' @param cm A contingency table or matrix of predicted vs observed classes with
#'   reference classes in columns and predicted classes in rows.
#' @details Use `mcc()` for 2x2 tables (binary classification). `mcc_multi()` is
#'   for multi-class classification with k x k tables and is calculated using
#'   Gorodkin's method.
#' @returns Returns a value between -1 and +1. A coefficient of +1
#' represents a perfect prediction, 0 no better than random prediction and -1
#' indicates total disagreement between prediction and observation.
#' @references
#' Gorodkin, J. (2004). \emph{Comparing two K-category assignments by a
#' K-category correlation coefficient}. Computational Biology and Chemistry. 28
#' (5): 367–374.
#' @export
mcc <- function(cm) {
  tp <- cm[2, 2]
  tn <- cm[1, 1]
  fp <- cm[2, 1]
  fn <- cm[1, 2]
  (tp*tn - fp*fn) /
    sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn))
}

#' @rdname mcc
#' @export
mcc_multi <- function(cm) {
  N <- sum(cm)
  tt <- colSums(cm)
  pp <- rowSums(cm)
  RK <- (N * sum(diag(cm)) - sum(tt * pp)) / (sqrt(N^2 - sum(tt^2)) * sqrt(N^2 - sum(pp^2)))
  RK
}
