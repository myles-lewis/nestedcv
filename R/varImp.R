
#' Extract coefficients from outer CV glmnet models
#' 
#' Extracts coefficients from outer CV glmnet models from a `nestcv.glmnet`
#' fitted object.
#' 
#' @param x a `nestcv.glmnet` fitted object
#' @return matrix of coefficients from outer CV glmnet models. Coefficients for
#'   variables which are not present in a particular outer CV fold model are set
#'   to 0.
#' @seealso [cv_varImp()]
#' @export
cv_coef <- function(x) {
  if (!inherits(x, "nestcv.glmnet")) stop("Not a `nestcv.glmnet` object")
  cfset <- lapply(x$outer_result, function(i) {
    i$coef[-1]
  })
  m <- list2matrix(cfset)
  mr <- rowMeans(m)
  m[order(abs(mr), decreasing = TRUE), ]
}


#' Extract coefficients from outer CV caret models
#' 
#' Extracts coefficients from outer CV glmnet models from a `nestcv.train`
#' fitted object.
#' 
#' @param x a `nestcv.train` fitted object
#' @return matrix of variable importance from outer CV caret models. Variable
#'   importance for variables which are not present in a particular outer CV
#'   fold model is set to 0.
#' @details
#' Note that [caret::varImp()] may require the model package to be fully loaded
#' in order to function. During the fitting process `caret` often only loads the
#' package by namespace.
#' @seealso [cv_coefs()]
#' @importFrom caret varImp
#' @export
cv_varImp <- function(x) {
  if (!inherits(x, "nestcv.train")) stop("Not a `nestcv.train` object")
  vset <- lapply(x$outer_result, function(i) {
    v <- varImp(i$fit)$importance
    if (ncol(v) > 2) {
      v <- rowMeans(v, na.rm = TRUE)
      return(v)
    } else if (ncol(v) == 2) {
      v <- v[,1, drop = FALSE]
    }
    setNames(unlist(v), rownames(v))
  })
  m <- list2matrix(vset)
  mr <- rowMeans(m, na.rm = TRUE)
  m[order(abs(mr), decreasing = TRUE), ]
}


list2matrix <- function(x, na.val = 0) {
  varset <- unique(unlist(lapply(x, names)))
  m <- matrix(na.val, nrow=length(varset), ncol=length(x),
              dimnames=list(varset, names(x)))
  for (i in seq_along(x)) {
    m[match(names(x[[i]]), varset), i] <- x[[i]]
  }
  m
}


#' Variable stability
#' 
#' Uses variable importance across models trained and tested across outer CV
#' folds to assess stability of variable importance.
#' 
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @param percent Logical for `nestcv.glmnet` objects only, whether to scale
#'   coefficients to percentage of the largest coefficient in each model
#' @param ... Optional arguments for compatibility
#' @return Dataframe containing mean, sd, sem of variable importance and
#'   frequency by which each variable is selected in outer folds.
#' @details
#' Note that for caret models [caret::varImp()] may require the model package to
#' be fully loaded in order to function. During the fitting process `caret`
#' often only loads the package by namespace.
#' @seealso [cv_coefs()] [cv_varImp()]
#' @export
var_stability <- function(x, ...) {
  UseMethod("var_stability")
}


#' @rdname var_stability
#' @importFrom stats sd
#' @export
var_stability.nestcv.glmnet <- function(x, percent = FALSE, ...) {
  m <- cv_coef(x)
  if (percent) {
    cm <- Rfast::colMaxs(m, value = TRUE)
    m <- t(t(m) / cm * 100)
  }
  mm <- rowMeans(m)
  msd <- apply(m, 1, sd)
  msem <- msd/sqrt(ncol(m))
  freq <- apply(m, 1, function(i) sum(i!=0))
  df <- data.frame(mean = mm, sd = msd, sem = msem, freq = freq, check.names = FALSE)
  df[order(abs(df$mean), decreasing = TRUE), ]
}


#' @rdname var_stability
#' @export
var_stability.nestcv.train <- function(x, ...) {
  m <- cv_varImp(x)
  mm <- rowMeans(m)
  msd <- apply(m, 1, sd)
  msem <- msd/sqrt(ncol(m))
  freq <- apply(m, 1, function(i) sum(i!=0))
  df <- data.frame(mean = mm, sd = msd, sem = msem, freq = freq, check.names = FALSE)
  df[df$freq > 0, ]
}


#' Plot variable stability
#'
#' Produces a ggplot2 plot of stability of variable importance across models trained and
#' tested across outer CV folds.
#'
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @param abs Logical for `nestcv.glmnet` objects only, whether to show absolute
#'   coefficients
#' @param final Logical whether to restrict variables to only those which ended
#'   up in the final fitted model or to include all variables selected across
#'   all outer folds.
#' @param top Limits number of variables plotted. Ignored if `final = TRUE`.
#' @param breaks Vector of continuous breaks for legend colour/size
#' @param percent Logical for `nestcv.glmnet` objects only, whether to scale
#'   coefficients to percentage of the largest coefficient in each model
#' @return A ggplot2 plot
#' @seealso [var_stability()]
#' @importFrom ggplot2 geom_vline geom_errorbarh scale_fill_distiller scale_size
#' @export
plot_var_stability <- function(x, abs = TRUE,
                               final = TRUE,
                               top = 25,
                               breaks = c(1,2,5,10),
                               percent = TRUE) {
  df <- var_stability(x, percent = percent)
  df$name <- factor(rownames(df), levels = rownames(df))
  if (final) {
    fv <- if (inherits(x, "nestcv.glmnet")) {
      names(coef(x))[-1]
    } else if (inherits(x, "nestcv.train")) {
      x$final_vars
    }
    if (!all(fv %in% rownames(df))) {
      message(paste(fv[!fv %in% rownames(df)], collapse = ", "),
              " not in final model")
      fv <- fv[fv %in% rownames(df)]
    }
    df <- df[rownames(df) %in% fv, ]
  } else {
    top <- min(top, nrow(df))
    df <- df[1:top, ]
  }
  if (abs) df$mean <- abs(df$mean)
  xtitle <- if (!percent & inherits(x, "nestcv.glmnet")) {
    if (abs) "Mean absolute coefficient" else "Mean coefficient"
  } else "Mean variable importance"
  ggplot(df, aes(x = .data$mean, y = .data$name)) +
    (if (!abs) geom_vline(xintercept = 0)) +
    geom_errorbarh(aes(xmin = .data$mean - .data$sem,
                       xmax = .data$mean + .data$sem), height = 0.2) +
    geom_point(aes(size = .data$freq,
                   fill = .data$freq), shape = 21) +
    scale_fill_distiller(guide = "legend", palette = "Blues", direction = 1,
                         breaks = breaks) +
    scale_size(guide = "legend", breaks = breaks) +
    scale_y_discrete(limits=rev) + ylab("") +
    xlab(xtitle) +
    theme_minimal() +
    theme(axis.text = element_text(colour = "black"))
}
