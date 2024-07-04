## R CMD check results

0 errors | 0 warnings | 0 note

* Important fix to correct the calculation of R-squared as a performance metric.
The previous version relied on caret which actually calculates Pearson 
correlation coefficient. This has been corrected.
* Bugfixes for edge cases such as a `x` being a single predictor instead of a 
dataframe.
* Fixed Rd package anchor links as requested.
