News
=====

# nestedcv 0.3.1
###### 15/09/2022
* Fix bug in `nestcv.train` for caret models with tuning parameters which are
factors
* Fix bug in `nestcv.train` for caret models using regression

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
