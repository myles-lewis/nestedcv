# SHAP analysis using shapr


#' Prediction wrappers to use shapr with nestedcv
#'
#' Prediction wrapper functions to enable the use of the `shapr` package for
#' generating SHAP values from `nestedcv` trained models.
#'
#' @param x a `nestcv.glmnet` or `nestcv.train` object
#' @param newdata a matrix or data frame of new data
#' @param cl integer representing which class to predict
#' @return prediction wrapper function designed for use with
#'   [nestcv.explain()] or [shapr::explain()]
#' @details
#' These prediction wrapper functions are designed to be used with the
#' `shapr` package via [nestcv.explain()]. The functions
#' `pred_nestcv_glmnet` and `pred_train` work for `nestcv.glmnet` and
#' `nestcv.train` models respectively for either binary classification or
#' regression.
#'
#' For multiclass classification use `pred_nestcv_glmnet_class(1)`,
#' `pred_nestcv_glmnet_class(2)` etc for each class. Similarly
#' `pred_train_class(1)`, `pred_train_class(2)` etc for [nestcv.train]
#' objects.
#'
#' @examples
#' if (requireNamespace("shapr") & requireNamespace("mlbench")) {
#'   library(shapr)
#'
#'   # Boston housing dataset
#'   library(mlbench)
#'   data(BostonHousing2)
#'   dat <- BostonHousing2
#'   y <- dat$cmedv
#'   x <- subset(dat, select = -c(cmedv, medv, town, chas))
#'
#'   # Fit a glmnet model using nested CV
#'   # Only 3 outer CV folds and 1 alpha value for speed
#'   fit <- nestcv.glmnet(y, x, family = "gaussian", n_outer_folds = 3, alphaSet = 1)
#'
#'   # Generate SHAP values using shapr
#'   # Reduced `max_n_coalitions` for demo; best left as default (NULL)
#'   sh <- nestcv.explain(fit, pred_nestcv_glmnet,
#'                        max_n_coalitions = 20)
#'
#'   # Plot overall variable importance
#'   plot_shap_bar(sh, x)
#'
#'   # Plot beeswarm plot
#'   plot_shap_beeswarm(sh, x, size = 1)
#' }
#' @export
#'
pred_nestcv_glmnet <- function(x, newdata) {
  predict(x, newdata)[,1]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_nestcv_glmnet_class <- function(cl) {
  function(x, newdata) {
    predict(x, newdata)[, cl, 1]
  }
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_train <- function(x, newdata) {
  newdata <- as.data.frame(newdata)
  if (x$final_fit$modelType == "Classification") {
    predict(x, newdata, type="prob")[,2]
  } else {
    predict(x, newdata)
  }
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_train_class <- function(cl) {
  function(x, newdata) {
    newdata <- as.data.frame(newdata)
    predict(x, newdata, type="prob")[,cl]
  }
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_SuperLearner <- function(x, newdata) {
  newdata <- as.data.frame(newdata)
  predict(x, newdata)$pred[,1]
}


select_approach <- function(x_train, min_obs_per_feature) {
  is_categorical <- vapply(
    x_train, function(col) is.factor(col) || is.character(col), logical(1)
  )
  if (all(is_categorical)) {
    return("categorical")
  }
  if (any(is_categorical)) {
    return("ctree")
  }
  if (nrow(x_train) / ncol(x_train) >= min_obs_per_feature) {
    return("empirical")
  }
  return("independence")
}

#' Generate SHAP values from nestedcv models using shapr
#'
#' Convenience wrapper around [shapr::explain()] that works with `nestedcv`
#' fitted models. Returns the `shapr` object directly, compatible with
#' [plot_shap_bar()], [plot_shap_beeswarm()], and `shapr`'s own `print()`/
#' `plot()` methods.
#'
#' @param model A `nestcv.glmnet`, `nestcv.train`, or other `nestedcv` model
#'   object.
#' @param predict_model Prediction wrapper function with signature
#'   `function(model, newdata)` returning a numeric vector of predictions.
#'   Use [pred_nestcv_glmnet()], [pred_train()], [pred_nestcv_glmnet_class()],
#'   [pred_train_class()], or [pred_SuperLearner()] as appropriate.
#' @param x_explain A matrix or data frame of feature values to compute SHAP
#'   values for.
#' @param x_train A matrix or data frame of feature values used as the
#'   background training data. Defaults to `model$xsub[, model$final_vars]`
#'   when `model` stores these (e.g. a fitted `nestcv.glmnet`/`nestcv.train`);
#'   otherwise falls back to `x_explain` (matching `fastshap` behaviour of
#'   using the same data for both).
#' @param approach Character string specifying the shapr estimation approach.
#'   Defaults to `NULL`, in which case it is auto-selected from `x_train`'s
#'   column types: `"categorical"` if all columns are factor/character,
#'   `"ctree"` if mixed numeric and factor/character, or for all-numeric
#'   `x_train`, `"empirical"` if there are at least `min_obs_per_feature`
#'   observations per feature, otherwise `"independence"`.
#'   Pass any other value accepted by [shapr::explain()] (e.g. `"vaeac"`,
#'   `"gaussian"`, `"copula"`) to override.
#' @param phi0 Numeric scalar; the baseline (null) prediction (i.e. the
#'   expected model output when no features are known). Defaults to `NULL`,
#'   in which case it is automatically computed as
#'   `mean(predict_model(model, x_explain))` on the supplied data. For
#'   regression this equals `mean(y_train)`; for classification it equals the
#'   mean predicted probability. Override this argument if you want to use a
#'   different reference value, e.g. computed on a held-out set.
#' @param min_obs_per_feature Numeric scalar specifying the minimum ratio of
#'   `nrow(x_train)` to `ncol(x_train)` required to auto-select `"empirical"`
#'   for all-numeric `x_train` (see `approach` above). Defaults to `20`.
#'   Ignored if `approach` is supplied explicitly.
#' @param ... Additional arguments passed to [shapr::explain()], e.g.
#'   `verbose`.
#' @return the `shapr` object returned by [shapr::explain()]. Pass it
#'   directly to [plot_shap_bar()] or [plot_shap_beeswarm()], or use
#'   `shapr`'s own `print()`/`plot()` methods on it.
#' @importFrom stats predict
#' @export
nestcv.explain <- function(model, predict_model,
                           x_explain = NULL,
                           x_train  = x_explain,
                           approach = NULL,
                           phi0     = NULL,
                           min_obs_per_feature = 20,
                           ...) {
  if (!requireNamespace("shapr", quietly = TRUE)) {
    stop("Package 'shapr' must be installed to use nestcv.explain()", call. = FALSE)
  }
  
  if (is.null(x_explain)) {
    if (!is.null(model$xsub) && is.character(model$final_vars)) {
      x_explain <- model$xsub[, model$final_vars, drop = FALSE]
    }
  }
  x_explain <- as.data.frame(x_explain)
  x_train <- as.data.frame(x_train)

  if (is.null(approach)) {
    approach <- select_approach(x_train, min_obs_per_feature)
  }

  # Compute phi0 if not supplied
  if (is.null(phi0)) {
    phi0 <- mean(model$y, na.rm = TRUE)
    # for classification
    if (is.na(phi0)) phi0 <- mean(predict_model(model, x_explain), na.rm = TRUE)
  }
  
  shapr::explain(
    model         = model,
    x_explain     = x_explain,
    x_train       = x_train,
    approach      = approach,
    phi0          = phi0,
    predict_model = predict_model,
    ...
  )
}


#' SHAP importance beeswarm plot
#'
#' @param shap a matrix of SHAP values, or the `shapr` object returned by
#'   [nestcv.explain()]/[shapr::explain()]
#' @param x a matrix or dataframe of feature values containing only features
#'   values from the training data. The rows must match rows in `shap`. If a
#'   dataframe is supplied it is converted to a numeric matrix using
#'   [data.matrix()].
#' @param cex Scaling for adjusting point spacing. See
#'   `ggbeeswarm::geom_beeswarm()`.
#' @param corral String specifying method used to corral points. See
#'   `ggbeeswarm::geom_beeswarm()`.
#' @param corral.width Numeric specifying width of corral, passed to
#'   `geom_beeswarm`
#' @param scheme Colour scheme as a vector of 3 colours
#' @param sort Logical whether to sort predictors by mean absolute SHAP value.
#' @param top Sets a limit on the number of variables plotted or `NULL` to plot
#'   all variables. If `top` is set then variables are sorted and `sort` is
#'   overrode.
#' @param ... Other arguments passed to `ggbeeswarm::geom_beeswarm()` e.g.
#'   `size`.
#' @return A ggplot2 plot
#' @importFrom ggplot2 scale_color_gradient2 guide_colorbar rel element_line element_text
#' @importFrom utils stack
#' @export
#'
plot_shap_beeswarm <- function(shap, x,
                               cex = 0.25,
                               corral = "random",
                               corral.width = 0.7,
                               scheme = c("deepskyblue2", "purple3", "red"),
                               sort = TRUE,
                               top = NULL, ...) {
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Package 'ggbeeswarm' must be installed", call. = FALSE)
  }
  if (inherits(shap, "shapr")) {
    sv_dt <- as.data.frame(shap$shapley_values_est)
    keep_cols <- setdiff(names(sv_dt), c("explain_id", "none"))
    shap <- as.matrix(sv_dt[, keep_cols, drop = FALSE])
    x <- data.matrix(x[, colnames(shap), drop = FALSE])
  } else {
    shap <- as.matrix(shap)
  }  
  if (!identical(dim(shap), dim(x))) stop("`shap` and `x` are misaligned")
  meanshap <- colMeans(abs(as.matrix(shap)))
  zeros <- if (sort) meanshap == 0 else FALSE
  if (any(zeros)) {
    message("Variables with mean(|SHAP|)=0: ",
            paste(colnames(shap)[zeros], collapse = ", "))
  }
  keep <- !zeros
  if (!is.null(top)) {
    ord <- order(meanshap, decreasing = TRUE)
    if (top < length(meanshap)) keep <- ord[1:top]
  }
  shap_stack <- stack(as.data.frame(shap[, keep]))
  x_stack <- stack(as.data.frame(clip_scale(x[, keep])))
  df <- data.frame(variable = shap_stack$ind, SHAP = shap_stack$values,
                   val = x_stack$values)

  if (sort) {
    ord <- sort(meanshap[keep], decreasing = TRUE)
    df$variable <- factor(df$variable, levels = names(ord))
  }

  ggplot(df, aes(y=.data$variable, x=.data$SHAP, col=.data$val)) +
    geom_vline(xintercept = 0) +
    ggbeeswarm::geom_beeswarm(cex = cex, corral = corral,
                              corral.width = corral.width,
                              orientation = "y", ...) +
    scale_color_gradient2(low = scheme[1], mid = scheme[2], high = scheme[3],
                          breaks = c(-1.5, 1.5),
                          labels = c("Low", "High"), name="Feature value",
                          guide = guide_colorbar(
                            barwidth = 0.5,
                            barheight = 8,
                            title.position = "left")
                          ) +
    scale_y_discrete(limits = rev) +
    ylab("") +
    xlab("SHAP value") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"),
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0,
                                                       size = rel(0.9)))
}


clip_scale <- function(x) {
  scx <- scale(x)
  scx[which(scx < -1.5, arr.ind = TRUE)] <- -1.5
  scx[which(scx > 1.5, arr.ind = TRUE)] <- 1.5
  scx
}


#' SHAP importance bar plot
#'
#' @param shap a matrix of SHAP values, or the `shapr` object returned by
#'   [nestcv.explain()]/[shapr::explain()]
#' @param x a matrix or dataframe of feature values containing only features
#'   values from the training data. The rows must match rows in `shap`. If a
#'   dataframe is supplied it is converted to a numeric matrix using
#'   [data.matrix()].
#' @param sort Logical whether to sort predictors by mean absolute SHAP value
#' @param labels Character vector of labels for directionality
#' @param top Sets a limit on the number of variables plotted or `NULL` to plot
#'   all variables. If `top` is set then variables are sorted and `sort` is
#'   overrode.
#' @return A ggplot2 plot
#' @importFrom ggplot2 geom_bar element_line element_text
#' @export
#'
plot_shap_bar <- function(shap, x,
                          sort = TRUE,
                          labels = c("Negative", "Positive"),
                          top = NULL) {
  if (inherits(shap, "shapr")) {
    sv_dt <- as.data.frame(shap$shapley_values_est)
    keep_cols <- setdiff(names(sv_dt), c("explain_id", "none"))
    shap <- sv_dt[, keep_cols, drop = FALSE]
    x <- x[, colnames(shap), drop = FALSE]
  } else {
    shap <- as.matrix(shap)
  }
  if (!identical(dim(shap), dim(x))) stop("`shap` and `x` are misaligned")
  meanshap <- colMeans(abs(as.matrix(shap)))
  x <- data.matrix(x)
  cor1 <- diag(suppressWarnings(cor(shap, x)))
  sign1 <- sign(cor1)
  sign1[is.na(sign1)] <- 1
  sign1 <- factor(sign1, levels = c(-1, 1), labels = labels)
  df <- data.frame(var = names(meanshap), meanshap = meanshap,
                   Direction = sign1)
  if (!is.null(top) && top < ncol(x)) {
    ord <- order(meanshap, decreasing = TRUE)[1:top]
    df <- df[ord, ]
  } else if (sort) {
    zeros <- meanshap == 0
    if (any(zeros)) {
      message("Variables with mean(|SHAP|)=0: ",
              paste(colnames(shap)[zeros], collapse = ", "))
    }
    df <- df[!zeros, ]
    ord <- order(df$meanshap, decreasing = TRUE)
    df <- df[ord, ]
  }
  df$var <- factor(df$var, levels = rev(df$var))

  ggplot(df, aes(y = .data$var, x = .data$meanshap, fill = .data$Direction)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = c("royalblue", "red")) +
    xlab("mean(|SHAP|)") +
    ylab("") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"))
}
