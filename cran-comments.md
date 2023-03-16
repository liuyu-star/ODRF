## Update ODRF version to 0.0.3

We have optimized the ODRF package and made significant changes. Therefore we request to update the ODRF version. The main changes as follows.

* The function predicate.ODT() runs error when ODT is not split (depth=1), and we have fixed this bug.
* We have fixed the function predict.ODRF with arguments numOOB and weight.tree related issues.
* We have fixed the functions plot.ODT(), VarImp() and plot.VarImp().
* We have fixed the argument 'lambda' of the functions ODT() and ODRF().


## R CMD check results

0 errors | 0 warnings | 0 note

* This is the second version of the package. It currently has no ERRORs, WARNINGs, or NOTEs from devtools::check().


## Reverse dependencies

* We have run R CMD check on the downstream dependency and no problems were found related to this package.


---
Thanks!

Yu Liu and Yingcun Xia
