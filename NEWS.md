News
=====

# nestedcv 0.6.9
###### 15/08/2023

## New features

* Added argument `verbose` in `nestcv.train()` and `nestcv.glmnet()` to show 
progress.
* Added argument `multicore_fork` in `nestcv.train()` to allow choice of 
parallelisation between forked multicore processing using `mclapply` or 
non-forked using `parLapply`. This can help prevent errors with some 
multithreaded caret models e.g. `model = "xgbTree"` when `tuneLength` is set to 
a high number.
* In `one_hot()` changed `all_levels` argument default to `FALSE` to be 
compatible with regression models by default.

# nestedcv 0.6.8
###### 16/07/2023

## New features

* Add coefficient column to `lm_filter()` full results table

## Bug fixes

* Fixed significant bug in `lm_filter()` where variables with zero variance were
incorrectly reporting very low p-values in linear models instead of returning
`NA`. This is due to how rank deficient models are handled by
`RcppEigen::fastLmPure`. Default method for `fastLmPure` has been changed to `0`
to allow detection of rank deficient models.
* Fixed bug in `weight()` caused by `NA`. Allow `weight()` to tolerate character 
vectors.

# nestedcv 0.6.7
###### 01/07/2023

## New features

* Better handling of dataframes in filters. `keep_factors` option has been added 
to filters to control filtering of factors with 3 or more levels.
* Added `one_hot()` for fast one-hot encoding of factors and character columns 
by creating dummy variables.
* Added `stat_filter()` which applies univariate filtering to dataframes with 
mixed datatype (continuous & categorical combined).
* Changed one-way ANOVA test in `anova_filter()` from `Rfast::ftests()` to 
`matrixTests::col_oneway_welch()` for much better accuracy

## Bug fixes

* Fixed bug caused by use of weights with `nestcv.train()` (Matt Siggins 
suggestion)

# nestedcv 0.6.6
###### 07/06/2023

## New features

* Added `n_inner_folds` argument to `nestcv.train()` to make it easier to set
the number of inner CV folds, and `inner_folds` argument which enables setting
the inner CV fold indices directly (suggestion Aline Wildberger)

## Bug fixes

* Fixed error in `plot_shap_beeswarm()` caused by change in fastshap 0.1.0 output 
from tibble to matrix
* Fixed bug with categorical features and `nestcv.train()`

# nestedcv 0.6.4
###### 29/05/2023

## New features

* Add argument `pass_outer_folds` to both `nestcv.glmnet` and `nestcv.train`: 
this enables passing of passing of outer CV fold indices stored in `outer_folds` 
to the final round of CV. Note this can only work if `n_outer_folds` = number of 
inner CV folds and balancing is not applied so that `y` is a consistent length.

## Bug fixes

* Fix: ensure `nfolds` for final CV equals `n_inner_folds` in `nestcv.glmnet()`

# nestedcv 0.6.3
###### 17/05/2023
* Improve `plot_var_stability()` to be more user friendly
* Add `top` argument to shap plots

# nestedcv 0.6.2
###### 15/05/2023
* Modified examples and vignette in anticipation of new version of fastshap 0.1.0

# nestedcv 0.6.1
###### 15/04/2023
* Add vignette for variable stability and SHAP value analysis
* Refine variable stability and shap plots

# nestedcv 0.6.0
###### 19/03/2023
* Switch some packages from Imports to Suggests to make basic installation 
simpler.
* Provide helper prediction wrapper functions to make it easier to use package 
`fastshap` for calculating SHAP values.
* Add `force_vars` argument to `glmnet_filter()`
* Add `ranger_filter()`

# nestedcv 0.5.2
###### 17/02/2023
* Disable printing in `nestcv.train()` from models such as `gbm`. This fixes 
multicore bug when using standard R gui on mac/linux.
* Bugfix if `nestcv.glmnet()` model has 0 or 1 coefficients.
* Add multiclass AUC for multinomial classification.

# nestedcv 0.5.0
###### 23/01/2023
* `nestedcv` models now return `xsub` containing a subset of the predictor
matrix `x` with filtered variables across outer folds and the final fit
* `boxplot_model()` no longer needs the predictor matrix to be specified as it 
is contained in `xsub` in `nestedcv` models
* `boxplot_model()` now works for all `nestedcv` model types
* Add function `var_stability()` to assess variance and stability of variable 
importance across outer folds, and directionality for binary outcome
* Add function `plot_var_stability()` to plot variable stability across outer 
folds
* Add `finalCV = NA` option which skips fitting the final model completely. This
gives a useful speed boost if performance metrics are all that is needed.
* `model` argument in `outercv` now prefers a character value instead of a 
function for the model to be fitted
* Bugfixes

# nestedcv 0.4.6
###### 07/12/2022
* Add check model exists in `outercv`
* Perform final model fit first in `nestcv.train` which improves error detection
in caret. So `nestcv.train` can be run in multicore mode straightaway.
* Removes predictors with variance = 0
* Fix bug caused by filter p-values = NA

# nestedcv 0.4.4
###### 05/12/2022
* Add confusion matrix to results summaries for classification
* Fix bugs in extraction of inner CV predictions for `nestcv.glmnet`
* Fix multinomial `nestcv.glmnet`
* Add `outer_train_predict` argument to enable saving of predictions on outer
training folds
* Add function `train_preds` to obtain outer training fold predictions
* Add function `train_summary` to show performance metrics on outer training
folds

# nestedcv 0.4.1
###### 12/11/2022
* Add examples of imbalance datasets
* Fix rowname bug in `smote()`

# nestedcv 0.4.0
###### 28/09/2022
* Add support for nested CV on ensemble models from `SuperLearner` package
* Final CV on whole data is now the default in `nestcv.train` and
`nestcv.glmnet`

# nestedcv 0.3.6
###### 18/09/2022
* Fix windows parallelisation bugs

# nestedcv 0.3.5
###### 16/09/2022
* Fix bug in `nestcv.train` for caret models with tuning parameters which are
factors
* Fix bug in `nestcv.train` for caret models using regression
* Add option in `nestcv.train` and `nestcv.glmnet` to tune final model
parameters using a final round of CV on the whole dataset
* Fix bugs in LOOCV
* Add balancing to final model fitting
* Add case weights to `nestcv.train` and `outercv`

# nestedcv 0.3.0
###### 07/09/2022
* Add `randomsample()` to handle class imbalance using random over/undersampling
* Add `smote()` for SMOTE algorithm for increasing minority class data
* Add bootstrap wrapper to filters, e.g. `boot_ttest()`

# nestedcv 0.2.6
###### 02/09/2022
* Final lambda in `nestcv.glmnet()` is mean of best lambdas on log scale
* Added `plot_varImp` for plotting variable importance for `nestcv.glmnet` final
models

# nestedcv 0.2.4
###### 19/07/2022
* Corrected handling of multinomial models in `nestcv.glmnet()`
* Align lambda in `cva.glmnet()`
* Improve plotting of error bars in `plot.cva.glmnet`
* Bugfix: plot of single `alphaSet` in `plot.cva.glmnet`
* Updated documentation and vignette

# nestedcv 0.2.1
###### 15/06/2022

* Parallelisation on windows added
* hsstan model has been added (Athina Spiliopoulou)
* outer_folds can be specified for consistent model comparisons
* Checks on x, y added
* NA handling
* summary and print methods
* Implemented LOOCV
* Collinearity filter
* Implement lm and glm as models in outercv()
* Runnable examples have been added throughout

# nestedcv 0.0.9100
###### 02/03/2022

* Major update to include nestedcv.train function which adds nested CV to the 
`train` function of `caret`
* Note passing of extra arguments to filter functions specified by `filterFUN`
is no longer done through `...` but with a list of arguments passed through a
new argument `filter_options`.

# nestedcv 0.0.9003
###### 02/03/2022

* Initial build of nestedcv
* Added outercv.rf function for measuring performance of rf
* Added cv.rf for tuning mtry parameter
* Added plot_caret for plotting caret objects with error bars on the tuning 
metric
