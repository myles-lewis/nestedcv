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
#' `fastshap` package. They currently work for both `nestcv.glmnet` and
#' `nestcv.train` models for binary classification or regression. For multiclass
#' classification they need to be recoded to specify which class is being
#' investigated.
#' @export
#' 
pred_nestcv_glmnet <- function(x, newdata) {
  predict(x, newdata)[,1]
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


#' SHAP importance beeswarm plot
#'
#' @param shap a matrix of SHAP values
#' @param x a matrix of feature values containing only features values from the
#'   training data. The rows must match rows in `shap`.
#' @param bee.cex Scaling for adjusting point spacing. See `cex` in
#'   `ggbeeswarm::geom_beeswarm()`.
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom reshape2 melt
#' @export
#' 
plot_shap_importance <- function(shap, x, bee.cex = 0.5) {
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("ggbeeswarm package is not installed", call. = FALSE)
  }
  if (!identical(dim(shap), dim(x))) stop("`shap` and `x` are misaligned")
  meanshap <- colMeans(abs(as.matrix(shap)))
  zeros <- meanshap == 0
  if (any(zeros)) {
    message(paste(colnames(shap[zeros]), collapse = ", "),
            " have mean abs SHAP value 0")
  }
  ord <- sort(meanshap[!zeros], decreasing = TRUE)
  df <- suppressMessages(melt(shap[, !zeros], value.name = "SHAP"))
  df$val <- melt(scale2(x[, !zeros]))[, "value"]
  df$variable <- factor(df$variable, levels = names(ord))
  
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
