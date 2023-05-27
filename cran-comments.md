## Update ODRF version to 0.0.4

We have fixed some known bugs and added some new functions for the ODRF package. Therefore we request to update the ODRF version. The main changes as follows.

* Fixed function VarImp(), adding the method of measuring the importance of variables with node purity, and now VarImp() can be used for both class ODT and ODRF. 
* When the argument “Xcat ! = 0”, i.e., the category variable in predictor X is transformed to one-of-K encode. however for the argument “NodeRotateFun=‘RotMatRF’ (‘RotMatRand’)“ run error, we have now fixed it. 
* Added predicted values of training data for class ODT and ODRF.
* Fixed issue related to function predict.ODRF() when argument "weight.tree = TRUE".
* Optimized some other known issues.


## R CMD check results

0 errors | 0 warnings | 0 note

* This is the second version of the package. It currently has no ERRORs, WARNINGs, or NOTEs from devtools::check().


## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


---
Thanks!

Yu Liu and Yingcun Xia
