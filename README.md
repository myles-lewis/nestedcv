# nestedcv

Nested cross-validation (CV) for the 'glmnet' package, including cross-validation 
of elastic net alpha parameter and filter functions for feature selection.

# Installation

Install from Github (requires API token).

```
devtools::install_github("myles-lewis/nestedcv", auth_token = "API token...")
library(nestedcv)
```

# Example

Using R4RA data to predict CDAI 50% response to rituximab.

```
# set up data
load("/Users/myles/R/R4RA_shiny_pw7_update/1_raw_data/R4RA_270821.RData")

index <- r4ra.meta$Outliers_Detected_On_PCA != "outlier" & r4ra.meta$Visit == 3 
          & !is.na(r4ra.meta$Visit)
metadata <- r4ra.meta[index, ]
dim(metadata)  # 133 individuals

medians <- Rfast::rowMedians(as.matrix(r4ra.vst))
data <- t(as.matrix(r4ra.vst))
# remove low expressed genes
data <- data[index, medians > 6]  
dim(data)  # 16254 genes

# Rituximab cohort only
yrtx <- metadata$CDAI.response.status.V7[metadata$Randomised.medication == "Rituximab"]
yrtx <- factor(yrtx)
data.rtx <- data[metadata$Randomised.medication == "Rituximab", ]

# no filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0.5,
                         family = "binomial", cores = 8,
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
coef(res.rtx)

# t-test filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = ttest_filter,
                         filter_options = list(nfilter = 300, p_cutoff = NULL),
                         family = "binomial", cores = 8,
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
coef(res.rtx)
```

Output from the nested CV with glmnet can be plotted to show how deviance is 
affected by alpha and lambda.

```
plot_alphas(res.rtx)
plot_lambdas(res.rtx)
```

The tuning of alpha for each outer fold can be plotted.

```
plot(res.rtx$outer_result[[1]]$cvafit, showLegend = "bottomright")

# scatter plot
plot(res.rtx$outer_result[[1]]$cvafit, type = 'p', showLegend = "bottomright")
```

ROC curves from left-out folds from both outer and inner CV can be plotted.

```
# Outer CV ROC
plot(res.rtx$roc)

# Inner CV ROC
rtx.inroc <- innercv_roc(res.rtx)
plot(rtx.inroc)
pROC::auc(rtx.inroc)
```

The overall expression level of each gene selected in the final model can be 
compared with a boxplot.

```
boxplot_model(res.rtx, data.rtx, ylab = "VST")
```

Other filters include wilcoxon (Mann-Whitney) test, Pearson or Spearman 
correlation for regression modelling, random forest and ReliefF filters.

```
# random forest filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0.5, filterFUN = rf_filter,
                         filter_options = list(nfilter = 300),
                         family = "binomial", cores = 8, 
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
coef(res.rtx)

# ReliefF algorithm filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = relieff_filter,
                         filter_options = list(nfilter = 300),
                         family = "binomial", cores = 8, 
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
coef(res.rtx)
```

Leave-one-out cross-validation (LOOCV) can be performed on the outer folds.

```
# outer LOOCV
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = ttest_filter,
                         filter_options = list(nfilter = 300, p_cutoff = NULL),
                         outer_method = "loocv",
                         family = "binomial", cores = 8,
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
coef(res.rtx)
```
