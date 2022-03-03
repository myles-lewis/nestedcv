# nestedcv

Nested cross-validation for the 'glmnet' package, including cross-validation 
of elastic net alpha parameter and filter functions for feature selection.

# Installation

Install from Github (requires API token).

```
devtools::install_github("myles-lewis/nestedcv", auth_token = "API token...")
library(nestedcv)
```

# Known issues

* The filter functions return a character vector of column names, rather than a 
true index. This causes problems if the column names of the matrix of predictors 
are pure numbers or non-unique.
