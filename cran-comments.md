## Update ODRF version to 0.0.5

Dear CRAN maintainers,

Thank you for your feedback, we cannot resolve the "checking CRAN incoming feasibility ... NOTE". We have fixed some known bugs and added some new functions for the ODRF package. Therefore we request to update the ODRF version. The main changes as follows.

- Added linear model tree. Specifically, use parameter Xsplit as the splitting variable for ODT and fit a linear model for each split using the function "gmlnet". The corresponding parameter is split="linear".
- Added the ensemble of ODT-based boosting trees，denoted by ODBT.
- Changed the parameter "leafnode" to "type" in the function predict.ODT(), and for classification tasks, added category probability output.
- Optimized some other known issues.

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
