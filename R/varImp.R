
#' Coefficients from outer CV glmnet models
#' 
#' Extracts coefficients from outer CV glmnet models from a `nestcv.glmnet`
#' fitted object.
#' 
#' @param x a `nestcv.glmnet` fitted object
#' @param level For multinomial models only, either an integer specifying which
#'   level of outcome is being examined, or the level can be specified as a
#'   character value
#' @return matrix of coefficients from outer CV glmnet models plus the final
#'   glmnet model. Coefficients for variables which are not present in a
#'   particular outer CV fold model are set to 0.
#' @seealso [cv_varImp()]
#' @export
cv_coef <- function(x, level = 1) {
  if (!inherits(x, "nestcv.glmnet")) stop("Not a `nestcv.glmnet` object")
  cfset <- lapply(x$outer_result, function(i) {
    if (is.list(i$coef)) {
      coef(i$cvafit)[[level]][,1][-1]  # multinomial
    } else i$coef[-1]
  })
  if (inherits(x$final_fit, "cv.glmnet")) {
    if (is.list(coef(x))) {
      cfset$Final <- coef(x)[[level]][-1]
    } else cfset$Final <- coef(x)[-1]
  }
  list2matrix(cfset)
}


#' Extract variable importance from outer CV caret models
#' 
#' Extracts variable importance or coefficients from outer CV glmnet models from
#' a `nestcv.train` fitted object.
#' 
#' @param x a `nestcv.train` fitted object
#' @return matrix of variable importance from outer CV fold caret models as well
#'   as the final model. Variable importance for variables which are not present
#'   in a particular outer CV fold model is set to 0.
#' @details
#' Note that [caret::varImp()] may require the model package to be fully loaded
#' in order to function. During the fitting process `caret` often only loads the
#' package by namespace.
#' @seealso [cv_coef()]
#' @importFrom caret varImp
#' @export
cv_varImp <- function(x) {
  if (!inherits(x, "nestcv.train")) stop("Not a `nestcv.train` object")
  vset <- lapply(x$outer_result, function(i) {
    extractImp(i$fit)
  })
  if (inherits(x$final_fit, "train")) {
    vset$Final <- extractImp(x$final_fit)
  }
  list2matrix(vset)
}


# extract importance from a caret model
extractImp <- function(x) {
  v <- caret::varImp(x)$importance
  if (ncol(v) > 2) {
    v <- rowMeans(v, na.rm = TRUE)
    return(v)
  } else if (ncol(v) == 2) {
    v <- v[,1, drop = FALSE]
  }
  setNames(unlist(v), rownames(v))
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
#' folds to assess stability of variable importance. For glmnet, variable
#' importance is measured as the absolute model coefficients, optionally scaled
#' as a percentage. The frequency with which each variable is selected in outer
#' folds as well as the final model is also returned which is helpful for sparse
#' models or with filters to determine how often variables end up in the model
#' in each fold. For glmnet, the direction of effect is taken directly from the
#' sign of model coefficients. For `caret` models, direction of effect is not
#' readily available, so as a substitute, the directionality of each predictor
#' is determined by the function [var_direction()] using the sign of a t-test
#' for binary classification or the sign of regression coefficient for
#' continuous outcomes (not available for multiclass caret models). To better
#' understand direction of effect of each predictor within the final model, we
#' recommend using SHAP values - see the vignette "Explaining nestedcv models
#' with Shapley values". See [pred_train()] for an example.
#' 
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @param percent Logical for `nestcv.glmnet` objects only, whether to scale
#'   coefficients to percentage of the largest coefficient in each model
#' @param level For multinomial `nestcv.glmnet` models only, either an integer
#'   specifying which level of outcome is being examined, or the level can be
#'   specified as a character value
#' @param sort Logical whether to sort variables by mean importance
#' @param ... Optional arguments for compatibility
#' @return Dataframe containing mean, sd, sem of variable importance and
#'   frequency by which each variable is selected in outer folds.
#' @details
#' Note that for caret models [caret::varImp()] may require the model package to
#' be fully loaded in order to function. During the fitting process `caret`
#' often only loads the package by namespace.
#' @seealso [cv_coef()] [cv_varImp()] [pred_train()]
#' @export
var_stability <- function(x, ...) {
  UseMethod("var_stability")
}


#' @rdname var_stability
#' @importFrom stats sd
#' @export
var_stability.nestcv.glmnet <- function(x,
                                        percent = TRUE,
                                        level = 1,
                                        sort = TRUE, ...) {
  m <- cv_coef(x, level)
  vdir <- sign(rowMeans(m))
  if (percent) {
    # cm <- Rfast::colMaxs(m, value = TRUE)
    cm <- colSums(abs(m))
    if (any(cm == 0)) {
      m <- m[, cm != 0]
      cm <- cm[cm != 0]
    }
    m <- t(t(abs(m)) / cm * 100)
  }
  mm <- rowMeans(m)
  msd <- apply(m, 1, sd)
  msem <- msd/sqrt(ncol(m))
  freq <- apply(m, 1, function(i) sum(i!=0))
  df <- data.frame(mean = mm, sd = msd, sem = msem, frequency = freq,
                   check.names = FALSE)
  df$sign <- vdir[rownames(df)]
  df$direction <- factor(df$sign, levels = c(-1, 1),
                         labels = c("negative", "positive"))
  df$final <- "no"
  df$final[m[, "Final"] != 0] <- "yes"
  df$final <- factor(df$final)
  if (!sort) return(df)
  df[order(abs(df$mean), decreasing = TRUE), ]
}


#' @rdname var_stability
#' @export
var_stability.nestcv.train <- function(x,
                                       sort = TRUE, ...) {
  m <- cv_varImp(x)
  mm <- rowMeans(m)
  msd <- apply(m, 1, sd)
  msem <- msd/sqrt(ncol(m))
  freq <- apply(m, 1, function(i) sum(i!=0))
  df <- data.frame(mean = mm, sd = msd, sem = msem, frequency = freq, check.names = FALSE)
  vdir <- var_direction(x)
  if (!is.null(vdir)) {
    df$sign <- vdir[rownames(df)]
    df$direction <- if (nlevels(x$y) == 2) {
      factor(df$sign, levels = c(-1, 1), labels = paste("Up in", levels(x$y)))
    } else {
      factor(df$sign, levels = c(-1, 1), labels = c("negative", "positive"))
    }
  }
  df$final <- "no"
  df$final[rownames(df) %in% x$final_vars] <- "yes"
  df$final <- factor(df$final)
  df <- df[df$freq > 0, ]
  if (!sort) return(df)
  df[order(abs(df$mean), decreasing = TRUE), ]
}


#' Plot variable stability
#'
#' Produces a ggplot2 plot of stability (as SEM) of variable importance across
#' models trained and tested across outer CV folds. Overlays frequency with
#' which variables are selected across the outer folds and optionally overlays
#' directionality for binary response outcome.
#'
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @param final Logical whether to restrict variables to only those which ended
#'   up in the final fitted model or to include all variables selected across
#'   all outer folds.
#' @param top Limits number of variables plotted. Set to `NULL` to plot all
#'   variables.
#' @param direction Integer controlling plotting of directionality for binary or
#'   regression models. `0` means no directionality is shown, `1` means
#'   directionality is overlaid as a colour, `2` means directionality is
#'   reflected in the sign of variable importance. Not available for multiclass
#'   caret models.
#' @param dir_labels Character vector for controlling the legend when
#'   `direction = 1`
#' @param scheme Vector of 2 colours for directionality when `direction = 1`
#' @param breaks Vector of continuous breaks for legend colour/size
#' @param percent Logical for `nestcv.glmnet` objects only, whether to scale
#'   coefficients to percentage of the largest coefficient in each model. If set
#'   to `FALSE`, model coefficients are shown and `direction` is ignored.
#' @param level For multinomial `nestcv.glmnet` models only, either an integer
#'   specifying which level of outcome is being examined, or the level can be
#'   specified as a character value.
#' @param sort Logical whether to sort by mean variable importance. Passed to
#'   [var_stability()].
#' @return A ggplot2 plot
#' @seealso [var_stability()]
#' @importFrom ggplot2 geom_vline geom_errorbarh scale_fill_distiller
#'   scale_fill_manual scale_radius xlim
#' @export
plot_var_stability <- function(x,
                               final = TRUE,
                               top = NULL,
                               direction = 0,
                               dir_labels = NULL,
                               scheme = c("royalblue", "red"),
                               breaks = NULL,
                               percent = TRUE,
                               level = 1,
                               sort = TRUE) {
  df <- var_stability(x, percent = percent, level = level, sort = sort)
  df$name <- factor(rownames(df), levels = rownames(df))
  if (!sort) final <- FALSE
  if (final) {
    fv <- if (inherits(x, "nestcv.glmnet")) {
      if (is.list(coef(x))) {
        # multinomial
        names(coef(x)[[level]])[-1]
      } else names(coef(x))[-1]
    } else if (inherits(x, "nestcv.train")) {
      x$final_vars
    }
    if (!all(fv %in% rownames(df))) {
      message(paste(fv[!fv %in% rownames(df)], collapse = ", "),
              " not in final model")
      fv <- fv[fv %in% rownames(df)]
    }
    df <- df[rownames(df) %in% fv, ]
  }
  if (!is.null(top) && top < nrow(df)) df <- df[1:top, ]
  if (!percent & inherits(x, "nestcv.glmnet")) {
    if (direction == 1) {
      df$mean <- abs(df$mean)
      xtitle <- "|coefficient|"
    } else xtitle <- "Coefficient"
  } else xtitle <- "Variable importance"
  if (is.null(breaks)) {
    nof <- length(x$outer_folds)
    pr <- unique(round(pretty(c(1, nof), n = 4)))
    breaks <- setNames(c(pr, nof+1), c(as.character(pr), "all"))
  }
  
  if (direction != 0 && !"direction" %in% colnames(df)) {
    message("Missing directionality information")
    direction <- 0
  }
  if (direction == 2 && "sign" %in% colnames(df) && percent) {
    df$mean <- df$mean * df$sign
  }
  if (direction > 0 && !is.null(dir_labels)) {
    df$direction <- factor(df$direction, labels = dir_labels)
  }
  
  if (direction == 1) {
    ggplot(df, aes(x = .data$mean, y = .data$name)) +
      (if (min(df$mean) < 0) geom_vline(xintercept = 0)) +
      geom_errorbarh(aes(xmin = .data$mean - .data$sem,
                         xmax = .data$mean + .data$sem), height = 0.2) +
      geom_point(aes(size = .data$frequency,
                     fill = .data$direction), shape = 21) +
      scale_fill_manual(values = scheme, na.translate = FALSE) +
      scale_radius(breaks = breaks,
                   limits = c(1, NA)) +
      (if (min(df$mean) > 0) xlim(0, NA)) +
      scale_y_discrete(limits=rev) + ylab("") +
      xlab(xtitle) +
      theme_minimal() +
      theme(axis.text = element_text(colour = "black"))
  } else {
    ggplot(df, aes(x = .data$mean, y = .data$name)) +
      (if (min(df$mean) < 0) geom_vline(xintercept = 0)) +
      geom_errorbarh(aes(xmin = .data$mean - .data$sem,
                         xmax = .data$mean + .data$sem), height = 0.2) +
      geom_point(aes(size = .data$frequency,
                     fill = .data$frequency), shape = 21) +
      scale_fill_distiller(guide = "legend", palette = "Spectral", direction = 1,
                           breaks = breaks, limits = c(1, NA)) +
      scale_radius(guide = "legend", breaks = breaks,
                   limits = c(1, NA)) +
      (if (min(df$mean) > 0) xlim(0, NA)) +
      scale_y_discrete(limits=rev) + ylab("") +
      xlab(xtitle) +
      theme_minimal() +
      theme(axis.text = element_text(colour = "black"))
  }
}


#' Barplot variable stability
#'
#' Produces a ggplot2 plot of stability (as SEM) of variable importance across
#' models trained and tested across outer CV folds. Optionally overlays
#' directionality for binary response or regression outcomes.
#'
#' @param x a `nestcv.glmnet` or `nestcv.train` fitted object
#' @param final Logical whether to restrict variables to only those which ended
#'   up in the final fitted model or to include all variables selected across
#'   all outer folds.
#' @param top Limits number of variables plotted. Set to `NULL` to plot all
#'   variables.
#' @param direction Integer controlling plotting of directionality for binary or
#'   regression models. `0` means no directionality is shown, `1` means
#'   directionality is overlaid as a colour, `2` means directionality is
#'   reflected in the sign of variable importance. Not available for multiclass
#'   caret models.
#' @param dir_labels Character vector for controlling the legend when
#'   `direction = 1`
#' @param scheme Vector of 2 colours for directionality when `direction = 1`
#' @param breaks Vector of continuous breaks for legend colour/size
#' @param percent Logical for `nestcv.glmnet` objects only, whether to scale
#'   coefficients to percentage of the largest coefficient in each model. If set
#'   to `FALSE`, model coefficients are shown and `direction` is ignored.
#' @param level For multinomial `nestcv.glmnet` models only, either an integer
#'   specifying which level of outcome is being examined, or the level can be
#'   specified as a character value.
#' @param sort Logical whether to sort by mean variable importance. Passed to
#'   [var_stability()].
#' @return A ggplot2 plot
#' @seealso [var_stability()]
#' @importFrom ggplot2 geom_vline geom_errorbarh scale_fill_distiller
#'   scale_fill_manual scale_radius xlim
#' @export
barplot_var_stability <- function(x,
                                  final = TRUE,
                                  top = NULL,
                                  direction = 0,
                                  dir_labels = NULL,
                                  scheme = c("royalblue", "red"),
                                  breaks = NULL,
                                  percent = TRUE,
                                  level = 1,
                                  sort = TRUE) {
  df <- var_stability(x, percent = percent, level = level, sort = sort)
  df$name <- factor(rownames(df), levels = rownames(df))
  if (!sort) final <- FALSE
  if (final) {
    fv <- if (inherits(x, "nestcv.glmnet")) {
      if (is.list(coef(x))) {
        # multinomial
        names(coef(x)[[level]])[-1]
      } else names(coef(x))[-1]
    } else if (inherits(x, "nestcv.train")) {
      x$final_vars
    }
    if (!all(fv %in% rownames(df))) {
      message(paste(fv[!fv %in% rownames(df)], collapse = ", "),
              " not in final model")
      fv <- fv[fv %in% rownames(df)]
    }
    df <- df[rownames(df) %in% fv, ]
  }
  if (!is.null(top) && top < nrow(df)) df <- df[1:top, ]
  if (!percent & inherits(x, "nestcv.glmnet")) {
    if (direction == 1) {
      df$mean <- abs(df$mean)
      xtitle <- "|coefficient|"
    } else xtitle <- "Coefficient"
  } else xtitle <- "Variable importance"
  
  if (direction != 0 && !"direction" %in% colnames(df)) {
    message("Missing directionality information")
    direction <- 0
  }
  if (direction == 2 && "sign" %in% colnames(df) && percent) {
    df$mean <- df$mean * df$sign
  }
  if (direction > 0 && !is.null(dir_labels)) {
    df$direction <- factor(df$direction, labels = dir_labels)
  }
  
  if (direction == 0) {
    ggplot(df, aes(x = .data$mean, y = .data$name)) +
      (if (min(df$mean) < 0) geom_vline(xintercept = 0)) +
      geom_col(width = 0.75) +
      geom_errorbarh(aes(xmin = .data$mean - .data$sem,
                         xmax = .data$mean + .data$sem), height = 0.3) +
      scale_y_discrete(limits=rev) + ylab("") +
      xlab(xtitle) +
      theme_minimal() +
      theme(axis.text = element_text(colour = "black"))
  } else {
    ggplot(df, aes(x = .data$mean, y = .data$name)) +
      (if (min(df$mean) < 0) geom_vline(xintercept = 0)) +
      geom_col(aes(fill = .data$direction), width = 0.75) +
      geom_errorbarh(aes(xmin = .data$mean - .data$sem,
                         xmax = .data$mean + .data$sem), height = 0.3) +
      scale_fill_manual(values = scheme, na.translate = FALSE) +
      scale_y_discrete(limits = rev) + ylab("") +
      xlab(xtitle) +
      theme_minimal() +
      theme(axis.text = element_text(colour = "black"))
  }
}
