News
=====

# nestedcv 0.4.4
###### 05/12/2022
* Add contingency table to summary functions for classification
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
