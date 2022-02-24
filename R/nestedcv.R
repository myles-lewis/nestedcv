# nestedcv
# Myles Lewis
# 24-02-2022

# univariate t-test filter
unifilt <- function(y, data, nfilter = NULL, p_uni_cutoff = 0.05, 
                    return = "names") {
  indx1 <- as.numeric(y) == 1
  indx2 <- as.numeric(y) == 2
  res <- Rfast::ttests(data[indx1, ], data[indx2, ])
  rownames(res) <- colnames(data)
  if (return == "full") return(res[, 'stat'])
  out <- res[, 'pvalue']
  out <- sort(out[out < p_uni_cutoff])
  if (!is.null(nfilter)) out <- out[1:nfilter]
  names(out)
}

# obtains coefficients from a cv.glmnet model
# needs lambda to be specified
glmnet_coefs <- function(fit, lambda) {
  cf <- coef(fit, s = lambda)
  cf <- as.matrix(cf)
  cf <- cf[cf != 0, ]
  cf2 <- cf[-1]
  cf2 <- cf2[order(abs(cf2), decreasing = TRUE)]
  c(cf[1], cf2)  # keep intercept first
}

rf_filter <- function(y, x, nfilter = NULL, return = "names",
                      ntree = 1000,
                      mtry = ncol(x) * 0.2,
                      ...) {
  fit <- randomForest::randomForest(x, y, importance = TRUE,
                                    ntree = ntree, mtry = mtry, ...)
  vi <- as.vector(importance(fit, type = 1))
  names(vi) <- colnames(x)
  if (return == "full") return(vi)
  vi <- sort(vi, decreasing = TRUE)
  vi <- vi[vi != 0]
  if (!is.null(nfilter)) vi <- vi[1:min(nfilter, length(vi))]
  names(vi)
}

relieff_filter <- function(y, x, nfilter = NULL, 
                           estimator = "ReliefFequalK",
                           return = "names", ...) {
  df <- as.data.frame(x)
  df$y <- y
  ref <- CORElearn::attrEval('y', df, estimator = estimator, ...)
  names(ref) <- colnames(x)
  if (return == "full") return(ref)
  ref <- sort(ref, decreasing = TRUE)
  if (!is.null(nfilter)) ref <- ref[1:min(nfilter, length(ref))]
  names(ref)
}

combo_filter <- function(y, x, nfilter, return = "names", ...) {
  uni_set <- unifilt(y, x, nfilter, return = return)
  relf_set <- relieff_filter(y, x, nfilter, return = return, ...)
  if (return == "full") {
    return(list(unifilt = uni_set, relieff_filter = relf_set))
  }
  n <- round(nfilter / 2)
  unique(c(uni_set[1:n], relf_set[1:n]))
}

cv2.glmnet <- function(y, x,
                       filterFUN = NULL,
                       n_outer_folds = 10,
                       n_inner_folds = 10,
                       alphaSet = seq(0.8, 1, 0.05),
                       min_1se = 0,
                       keep_innerCV_pred = TRUE,
                       cores = 8, ...) {
  outer_folds <- createFolds(y, k = n_outer_folds, returnTrain = TRUE)
  outer_res <- mclapply(1:n_outer_folds, function(i) {
    trainIndex <- outer_folds[[i]]
    foldid <- sample(rep(seq_len(n_inner_folds), length = length(trainIndex)))
    # expand data with interactions
    filtx <- if (is.null(filterFUN)) x else {
      fset <- filterFUN(y[trainIndex], x[trainIndex, ], ...)
      x[, fset]
    }
    fit <- lapply(alphaSet, function(alpha) {
      cv.glmnet(x = filtx[trainIndex, ], y = y[trainIndex], family = "binomial", 
                alpha = alpha, nfolds = n_inner_folds, foldid = foldid, 
                keep = keep_innerCV_pred)
    })
    alphas <- unlist(lapply(fit, function(fitx) {
      w <- which.min(fitx$cvm)
      fitx$cvm[w]
    }))
    alpha.x <- which.min(alphas)
    s <- exp((log(fit[[alpha.x]]$lambda.min) * (1-min_1se) + log(fit[[alpha.x]]$lambda.1se) * min_1se))
    cf <- as.matrix(coef(fit[[alpha.x]], s = s))
    cf <- cf[cf != 0, ]
    # test on outer CV
    predy <- as.vector(predict(fit[[alpha.x]], newx = filtx[-trainIndex, ], s = s, type = "class"))
    predyp <- as.vector(predict(fit[[alpha.x]], newx = filtx[-trainIndex, ], s = s))
    preds <- data.frame(predy=predy, predyp=predyp, testy=y[-trainIndex])
    rownames(preds) <- rownames(x[-trainIndex, ])
    ret <- list(preds = preds,
                alpha = alphaSet[[alpha.x]], lambda = s, coef = cf, cvfit = fit[[alpha.x]],
                cv_alpha = alphas,
                nfilter = ncol(filtx))
    # inner CV predictions
    if (keep_innerCV_pred) {
      ind <- fit[[alpha.x]]$index["min", ]
      innerCV_preds <- fit[[alpha.x]]$fit.preval[, ind]
      ytrain <- y[trainIndex]
      ret <- append(ret, list(ytrain = ytrain, innerCV_preds = innerCV_preds))
    }
    ret
  }, mc.cores = cores)
  predslist <- lapply(outer_res, '[[', 'preds')
  output <- data.table::rbindlist(predslist)
  output <- as.data.frame(output)
  cm <- table(output$predy, output$testy)
  acc <- setNames(sum(diag(cm))/ sum(cm), "Accuracy")
  b_acc <- confusionMatrix(cm)$byClass[11]
  glmnet.roc <- pROC::roc(output$testy, output[, 2], direction = "<")
  auc <- glmnet.roc$auc
  # fit final glmnet
  lam <- mean(unlist(lapply(outer_res, '[[', 'lambda')))
  alph <- mean(unlist(lapply(outer_res, '[[', 'alpha')))
  filtx <- if (is.null(filterFUN)) x else {
    fset <- filterFUN(y, x, ...)
    x[, fset]
  }
  
  fit <- glmnet(filtx, y, family = "binomial", alpha = alph)
  list(output = output,
       outer_result = outer_res,
       mean_lambda = lam,
       mean_alpha = alph,
       final_fit = fit,
       roc = glmnet.roc,
       alphaSet = alphaSet,
       summary = setNames(c(auc, acc, b_acc), c("auc", "accuracy", "bal_accuracy")))
}

# plot cv alphas against deviance
plot_alphas <- function(cva, ...) {
  cv_alpha <- lapply(cva$outer_result, '[[', 'cv_alpha')
  coln <- length(cv_alpha)
  cols <- rainbow(coln)
  plot(cv_alpha[[1]], type = 'l', x = cva$alphaSet,
       ylim = range(unlist(cv_alpha)),
       xlab = 'Elastic net alpha',
       ylab = 'Binomial deviance',
       col = cols[1],
       las = 1, bty = 'l', ...)
  for (i in 2:10) {
    lines(cv_alpha[[i]], x = cva$alphaSet, col = cols[i])
  }
}

# extract roc from LO folds from inner CV
innercv_roc <- function(cva, direction = "<") {
  innerpreds <- unlist(lapply(cva$outer_result, '[[', 'innerCV_preds'))
  ytrain <- unlist(lapply(cva$outer_result, '[[', 'ytrain'))
  pROC::roc(ytrain, innerpreds, direction = direction)
}

