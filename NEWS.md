# ODRF 0.0.2.9000 (development version)

* The function predicate.ODT() runs error when ODT is not split (depth=1), and we have fixed this bug.
* We have fixed the function predict.ODRF with arguments numOOB and weight.tree related issues.
* We have fixed the functions VarImp() and plot.VarImp().
* We have modified the argument 'lambda' of the functions ODT() and ODRF().


# ODRF 0.0.2

* We have now explained CART and Random Forest in the description text.
* We have changed the Date field to a more recent date.
* We have now exported the functions RandRot() and defaults(), and no longer need ODRF:::
* We have removed par from plot.VarImp() and added on.exit to plot.prune.ODT(), and checked the code to make sure that it does not change the user's options, including par or working directory.
* We have removed the random seed number in functions ODRF(), poune.ODRF(), online.ODRF() and plot_ODT_depth().


# ODRF 0.0.1

* Added a `NEWS.md` file to track changes to the package.
* This is the first fully-functioning version of the package. It currently has no ERRORs, WARNINGs, or NOTEs from devtools::check().

