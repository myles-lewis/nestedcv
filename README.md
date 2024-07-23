# nestedcv

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/nestedcv)](https://cran.r-project.org/package=nestedcv)
[![Downloads](https://cranlogs.r-pkg.org/badges/nestedcv)](https://CRAN.R-project.org/package=nestedcv)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/nestedcv)](https://CRAN.R-project.org/package=nestedcv)

Nested cross-validation (CV) for the glmnet and caret packages. With glmnet this
includes cross-validation of elastic net alpha parameter. A number of feature
selection filter functions (t-test, Wilcoxon test, ANOVA, Pearson/Spearman
correlation, random forest, ReliefF) for feature selection are provided and can
be embedded within the outer loop of the nested CV. Nested CV can be also be
performed with the caret package giving access to the large number of prediction
methods available in caret.

# Installation

Install from CRAN
```
install.packages("nestedcv")
```

Install from Github
```
devtools::install_github("myles-lewis/nestedcv")
```

# Example

In this example using iris dataset (multinomial, 3 classes), we fit a glmnet
model, tuning both lambda and alpha with 10 x 10-fold nested CV.

```
library(nestedcv)
data(iris)
y <- iris$Species
x <- as.matrix(iris[, -5])

cores <- parallel::detectCores(logical = FALSE)  # detect physical cores

res <- nestcv.glmnet(y, x, family = "multinomial", cv.cores = cores)
summary(res)
```

Use `summary()` to see full information from the nested model fitting. `coef()`
can be used to show the coefficients of the final fitted model.
Filters can be used by setting the `filterFUN` argument. Options for the filter 
function are passed as a list through `filter_options`.

Output from the nested CV with glmnet can be plotted to show how deviance is 
affected by alpha and lambda.

```
plot_alphas(res)
plot_lambdas(res)
```

The tuning of lambda and alpha for each outer CV fold can be plotted. Here we
inspect outer CV fold 1.

```
plot(res$outer_result[[1]]$cvafit)
```

ROC curves from left-out folds from both outer and inner CV can be plotted for
binary comparisons (see vignette).

Nested CV can also be performed using the caret package framework. Here we use
caret for tuning random forest using the ranger package.

```
res <- nestcv.train(y, x, method = "ranger", cv.cores = cores)
summary(res)
```
