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
#' @examples
#' library(fastshap)
#' 
#' # Boston housing dataset
#' library(mlbench)
#' data(BostonHousing2)
#' dat <- BostonHousing2
#' y <- dat$cmedv
#' x <- subset(dat, select = -c(cmedv, medv, town, chas))
#' 
#' # Fit a glmnet model using nested CV
#' # Only 3 outer CV folds and 1 alpha value for speed
#' fit <- nestcv.glmnet(y, x, family = "gaussian", n_outer_folds = 3, alphaSet = 1)
#' 
#' # Generate SHAP values using fastshap::explain
#' # Only using 5 repeats here for speed, but recommend higher values of nsim
#' sh <- explain(fit, X=x, pred_wrapper = pred_nestcv_glmnet, nsim = 1)
#' 
#' # Plot overall variable importance
#' plot_shap_bar(sh, x)
#' 
#' # Plot beeswarm plot
#' plot_shap_beeswarm(sh, x, size = 1)
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
#' @param ... Other arguments passed to `ggbeeswarm::geom_beeswarm()`
#' @return A ggplot2 plot
#' @importFrom ggplot2 scale_color_gradient2 guide_colorbar rel
#' @importFrom reshape2 melt
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
  if (!identical(dim(shap), dim(x))) stop("`shap` and `x` are misaligned")
  meanshap <- colMeans(abs(as.matrix(shap)))
  zeros <- if (sort) meanshap == 0 else FALSE
  if (any(zeros)) {
    message("Variables with mean(|SHAP|)=0: ",
            paste(colnames(shap[zeros]), collapse = ", "))
  }
  keep <- !zeros
  if (!is.null(top)) {
    ord <- order(meanshap, decreasing = TRUE)
    if (top < length(meanshap)) keep <- ord[1:top]
  }
  df <- suppressMessages(melt(shap[, keep], value.name = "SHAP"))
  df$val <- melt(scale2(x[, keep]))[, "value"]
  if (sort) {
    ord <- sort(meanshap[keep], decreasing = TRUE)
    df$variable <- factor(df$variable, levels = names(ord))
  }
  
  ggplot(df, aes(y=.data$variable, x=.data$SHAP, col=.data$val)) +
    geom_vline(xintercept = 0) +
    ggbeeswarm::geom_beeswarm(cex = cex, corral = corral,
                              corral.width = corral.width,
                              ...) +
    scale_color_gradient2(low=scheme[1], mid=scheme[2], high=scheme[3],
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
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0,
                                                       size = rel(0.9)))
}


scale2 <- function(x) {
  scx <- scale(x)
  scx[which(scx < -1.5, arr.ind = TRUE)] <- -1.5
  scx[which(scx > 1.5, arr.ind = TRUE)] <- 1.5
  scx
}


#' SHAP importance bar plot
#'
#' @param shap a matrix of SHAP values
#' @param x a matrix of feature values containing only features values from the
#'   training data. The rows must match rows in `shap`.
#' @param sort Logical whether to sort predictors by mean absolute SHAP value
#' @param labels Character vector of labels for directionality
#' @param top Sets a limit on the number of variables plotted or `NULL` to plot
#'   all variables. If `top` is set then variables are sorted and `sort` is
#'   overrode.
#' @return A ggplot2 plot
#' @importFrom ggplot2 geom_bar
#' @export
#' 
plot_shap_bar <- function(shap, x,
                          sort = TRUE,
                          labels = c("Negative", "Positive"),
                          top = NULL) {
  if (!identical(dim(shap), dim(x))) stop("`shap` and `x` are misaligned")
  meanshap <- colMeans(abs(as.matrix(shap)))
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
              paste(colnames(shap[zeros]), collapse = ", "))
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
    theme(axis.text = element_text(colour = "black"))
}

