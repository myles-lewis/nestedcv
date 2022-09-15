# nestedcv

Nested cross-validation (CV) for the glmnet and caret packages. With glmnet this
includes cross-validation of elastic net alpha parameter. A number of filter
functions (t-test, Wilcoxon test, ANOVA, Pearson/Spearman correlation, random
forest, ReliefF) for feature selection are provided and can be embedded within
the outer loop of the nested CV. Nested CV can be also be performed with the
caret package giving access to the large number of prediction methods available
in caret.

# Installation

Install from CRAN
```
install.packages("nestedcv")
library(nestedcv)
```

Install from Github
```
devtools::install_github("myles-lewis/nestedcv")
```

# Example

Using R4RA data to predict CDAI 50% response to rituximab.

```
# set up data
load("/../R4RA_270821.RData")

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
                         family = "binomial", cv.cores = 4,
                         alphaSet = seq(0.7, 1, 0.05))
res.rtx
```

Use `summary()` to see full information from the nested model fitting. `coef()`
can be used to show the coefficients of the final fitted model.
Filters can be used by setting the `filterFUN` argument. Options for the filter 
function are passed as a list through `filter_options`.

```
# t-test filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = ttest_filter,
                         filter_options = list(nfilter = 300, p_cutoff = NULL),
                         family = "binomial", cv.cores = 4,
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
```

Output from the nested CV with glmnet can be plotted to show how deviance is 
affected by alpha and lambda.

```
plot_alphas(res.rtx)
plot_lambdas(res.rtx)
```

The tuning of alpha for each outer fold can be plotted.

```
plot(res.rtx$outer_result[[1]]$cvafit)

# scatter plot
plot(res.rtx$outer_result[[1]]$cvafit, type = 'p')

# number of non-zero coefficients
plot(res.rtx$outer_result[[1]]$cvafit, xaxis = 'nvar')
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

Other filters include Wilcoxon (Mann-Whitney) test, Pearson or Spearman 
correlation for regression modelling, random forest and ReliefF filters.

```
# random forest filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0.5, filterFUN = rf_filter,
                         filter_options = list(nfilter = 300),
                         family = "binomial", cv.cores = 4, 
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)

# ReliefF algorithm filter
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = relieff_filter,
                         filter_options = list(nfilter = 300),
                         family = "binomial", cv.cores = 4, 
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
```

Leave-one-out cross-validation (LOOCV) can be performed on the outer folds.

```
# outer LOOCV
res.rtx <- nestcv.glmnet(y = yrtx, x = data.rtx, min_1se = 0, filterFUN = ttest_filter,
                         filter_options = list(nfilter = 300, p_cutoff = NULL),
                         outer_method = "loocv",
                         family = "binomial", cv.cores = 4,
                         alphaSet = seq(0.7, 1, 0.05))
summary(res.rtx)
```

Nested CV can also be performed using the caret package framework. Here we use
caret for tuning glmnet.

```
# nested CV using caret
tg <- expand.grid(lambda = exp(seq(log(2e-3), log(1e0), length.out = 100)),
                  alpha = seq(0.8, 1, 0.1))
ncv <- nestcv.train(y = yrtx, x = data.rtx,
               method = "glmnet",
               savePredictions = "final",
               filterFUN = ttest_filter, filter_options = list(nfilter = 300),
               tuneGrid = tg, cv.cores = 4)
ncv$summary

# Plot ROC on outer folds
plot(ncv$roc)

# Plot ROC on inner LO folds
inroc <- innercv_roc(ncv)
plot(inroc)
auc(inroc)

# Show example tuning plot for outer fold 1
plot(ncv$outer_result[[1]]$fit, xTrans = log)

# Extract coefficients of final fitted model
glmnet_coefs(ncv$final_fit$finalModel, s = ncv$finalTune$lambda)
```


## Linear regression with hsstan (continuous outcome)

Cross-validation is used to apply univariate filtering of predictors.
Only one CV split is needed (outercv) as the Bayesian model does not require
learning of meta-parameters.

```
# specify options for running the code
nfolds <- 3

# specify number of cores for parallelising computation
# the product of cv.cores and mc.cores (12 cores) will be used in total
# number of cores for parallelising over CV folds
cv.cores <- 3
# number of cores for parallelising stan sampling (default over 4 chains)
# in hsstan the number of cores is specified by setting options(mc.cores)
options(mc.cores = 4)


# load iris dataset and simulate a continuous outcome
data(iris)
dt <- iris[, 1:4]
colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
dt <- as.data.frame(apply(dt, 2, scale))
dt$outcome.cont <- -3 + 0.5 * dt$marker1 + 2 * dt$marker2 + rnorm(nrow(dt), 0, 2)

# unpenalised covariates: always retain in the prediction model
uvars <- "marker1"
# penalised covariates: coefficients are drawn from hierarchical shrinkage prior
pvars <- c("marker2", "marker3", "marker4") # penalised covariates
# run cross-validation with univariate filter and hsstan
res.cv.hsstan <- outercv(y = dt$outcome.cont, x = dt[, c(uvars, pvars)],
                         model = model.hsstan,
                         filterFUN = lm_filter,
                         filter_options = list(force_vars = uvars,
                                               nfilter = 2,
                                               p_cutoff = NULL,
                                               rsq_cutoff = 0.9),
                         n_outer_folds = nfolds, cv.cores = cv.cores,
                         unpenalized = uvars, warmup = 1000, iter = 2000)
# view prediction performance based on testing folds
summary(res.cv.hsstan)

# load hsstan package to examine the Bayesian model
library(hsstan)
sampler.stats(res.cv.hsstan$final_fit)
print(projsel(res.cv.hsstan$final_fit), digits = 4) # adding marker2
```

Here adding `marker2` improves the model fit: substantial decrease of
KL-divergence from the full model to the submodel. Adding `marker3` does not
improve the model fit: no decrease of KL-divergence from the full model to the
submodel.

## Logistic regression with hsstan (binary outcome)

Cross-validation is used to apply univariate filtering of predictors.
Only one CV split is needed (outercv) as the Bayesian model does not require
learning of meta-parameters.

```
# sigmoid function
sigmoid <- function(x) {1 / (1 + exp(-x))}

# specify options for running the code
nfolds <- 3

# specify number of cores for parallelising computation
# the product of cv.cores and mc.cores will be used in total
# number of cores for parallelising over CV folds
cv.cores <- 3
# number of cores for parallelising stan sampling (default over 4 chains)
# in hsstan the number of cores is specified by setting options(mc.cores)
options(mc.cores = 4)

# load iris dataset and create a binary outcome
set.seed(267)
data(iris)
dt <- iris[, 1:4]
colnames(dt) <- c("marker1", "marker2", "marker3", "marker4")
dt <- as.data.frame(apply(dt, 2, scale))
rownames(dt) <- paste0("sample", c(1:nrow(dt)))
dt$outcome.bin <- sigmoid(0.5 * dt$marker1 + 2 * dt$marker2) > runif(nrow(dt))
dt$outcome.bin <- factor(dt$outcome.bin)


# unpenalised covariates: always retain in the prediction model
uvars <- "marker1"
# penalised covariates: coefficients are drawn from hierarchical shrinkage prior
pvars <- c("marker2", "marker3", "marker4") # penalised covariates
# run cross-validation with univariate filter and hsstan
res.cv.hsstan <- outercv(y = dt$outcome.bin,
                         x = as.matrix(dt[, c(uvars, pvars)]),
                         model = model.hsstan,
                         filterFUN = ttest_filter,
                         filter_options = list(force_vars = uvars,
                                               nfilter = 2,
                                               p_cutoff = NULL,
                                               rsq_cutoff = 0.9),
                         n_outer_folds = nfolds, cv.cores = cv.cores,
                         unpenalized = uvars, warmup = 1000, iter = 2000)


# view prediction performance based on testing folds
summary(res.cv.hsstan)

# load hsstan package to examine the Bayesian model
library(hsstan)
sampler.stats(res.cv.hsstan$final_fit)
print(projsel(res.cv.hsstan$final_fit), digits = 4) # adding marker2
```

Here adding `marker2` improves the model fit: substantial decrease of
KL-divergence from the full model to the submodel. Adding `marker3` does not
improve the model fit: no decrease of KL-divergence from the full model to the
submodel.
