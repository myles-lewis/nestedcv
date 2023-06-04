# nestedcv

Nested cross-validation (CV) for the glmnet and caret packages. With glmnet this
includes cross-validation of elastic net alpha parameter. A number of feature selection filter
functions (t-test, Wilcoxon test, ANOVA, Pearson/Spearman correlation, random
forest, ReliefF) for feature selection are provided and can be embedded within
the outer loop of the nested CV. Nested CV can be also be performed with the
caret package giving access to the large number of prediction methods available
in caret.

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

Using iris dataset (multinomial, 3 classes).

```
library(nestedcv)
data(iris)
y <- iris$Species
x <- as.matrix(iris[, -5])
res <- nestcv.glmnet(y, x, family = "multinomial",
                     alphaSet = (8:10)/10, cv.cores = 4)
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

The tuning of alpha for each outer fold can be plotted.

```
plot(res$outer_result[[1]]$cvafit)
```

ROC curves from left-out folds from both outer and inner CV can be plotted for binary comparisons (see vignette).

Nested CV can also be performed using the caret package framework. Here we use
caret for tuning glmnet.

```
# nested CV using caret
res <- nestcv.train(y, x, method="rf", cv.cores = 8)
summary(res)
```
