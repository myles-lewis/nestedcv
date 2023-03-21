# Help for shap analysis using fastshap


#' Prediction wrappers to use fastshap with nestedcv
#' 
#' Prediction wrapper functions to enable the use of the `fastshap` package for
#' generating SHAP values from `nestedcv` trained models.
#' 
#' @param x a `nestcv.glmnet` or `nestcv.train` object
#' @param newdata a matrix of new data
#' @return prediction wrapper function designed for use with
#'   [fastshap::explain()]
#' @details
#' These prediction wrapper functions are designed to be used with the
#' `fastshap` package. The functions `pred_nestcv_glmnet` and `pred_train` work
#' for `nestcv.glmnet` and `nestcv.train` models respectively for either binary
#' classification or regression.
#' 
#' For multiclass classification use `pred_nestcv_glmnet_class1`, `2` and `3`
#' for the first 3 classes. Similarly `pred_train_class1` etc for [nestcv.train]
#' objects. These functions can be inspected and easily modified to analyse
#' further classes.
#' 
#' @export
#' 
pred_nestcv_glmnet <- function(x, newdata) {
  predict(x, newdata)[,1]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_nestcv_glmnet_class1 <- function(x, newdata) {
  predict(x, newdata)[, 1, 1]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_nestcv_glmnet_class2 <- function(x, newdata) {
  predict(x, newdata)[, 2, 1]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_nestcv_glmnet_class3 <- function(x, newdata) {
  predict(x, newdata)[, 3, 1]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_train <- function(x, newdata) {
  if (x$final_fit$modelType == "Classification") {
    predict(x, newdata, type="prob")[,2]
  } else {
    predict(x, newdata)
  }
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_train_class1 <- function(x, newdata) {
  predict(x, newdata, type="prob")[,1]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_train_class2 <- function(x, newdata) {
  predict(x, newdata, type="prob")[,2]
}

#' @rdname pred_nestcv_glmnet
#' @export
pred_train_class3 <- function(x, newdata) {
  predict(x, newdata, type="prob")[,3]
}



#' SHAP importance beeswarm plot
#'
#' @param shap a matrix of SHAP values
#' @param x a matrix of feature values containing only features values from the
#'   training data. The rows must match rows in `shap`.
#' @param bee.cex Scaling for adjusting point spacing. See `cex` in
#'   `ggbeeswarm::geom_beeswarm()`.
#' @param sort Logical whether to sort predictors by mean absolute SHAP value.
#' @importFrom ggplot2 scale_color_gradient2 guide_colorbar
#' @importFrom reshape2 melt
#' @export
#' 
plot_shap_importance <- function(shap, x,
                                 bee.cex = 0.5, sort = TRUE) {
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Package 'ggbeeswarm' must be installed", call. = FALSE)
  }
  if (!identical(dim(shap), dim(x))) stop("`shap` and `x` are misaligned")
  meanshap <- colMeans(abs(as.matrix(shap)))
  zeros <- if (sort) meanshap == 0 else FALSE
  if (any(zeros)) {
    message(paste(colnames(shap[zeros]), collapse = ", "),
            " have mean abs SHAP value 0")
  }
  df <- suppressMessages(melt(shap[, !zeros], value.name = "SHAP"))
  df$val <- melt(scale2(x[, !zeros]))[, "value"]
  if (sort) {
    ord <- sort(meanshap[!zeros], decreasing = TRUE)
    df$variable <- factor(df$variable, levels = names(ord))
  }
  
  ggplot(df, aes(y=.data$variable, x=.data$SHAP, col=.data$val)) +
    geom_vline(xintercept = 0) +
    ggbeeswarm::geom_beeswarm(cex = bee.cex) +
    scale_color_gradient2(low="deepskyblue2", mid="purple3", high="red",
                          breaks = c(-1.5, 1.5),
                          labels = c("Low", "High"), name="Feature\nvalue\n",
                          guide = guide_colorbar(
                            barwidth = 0.5,
                            barheight = 8)
                          #   title.theme = element_text(angle = 90, hjust = 0.5, vjust = 0),
                          #   title.position = "left")
                          ) +
    scale_y_discrete(limits = rev) +
    ylab("") +
    xlab("SHAP value") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))
}


scale2 <- function(x) {
  scx <- scale(x)
  scx[which(scx < -1.5, arr.ind = TRUE)] <- -1.5
  scx[which(scx > 1.5, arr.ind = TRUE)] <- 1.5
  scx
}
