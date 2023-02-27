# ODRF 0.0.1

* Added a `NEWS.md` file to track changes to the package.
* This is the first fully-functioning version of the package. It currently has no ERRORs, WARNINGs, or NOTEs from devtools::check().

# ODRF 0.0.2

* We explained CART and Random Forest in the description text.
* We changed the Date field to a more recent date.
* We have exported the functions RandRot() and defaults(), and already replaced ODRF:::f with f.
* We removed par from plot.VarImp() and added on.exit to plot.prune.ODT(). so we made sure not to change the user's options, par or working directory.
* We removed the random seed number in functions ODRF(), poune.ODRF(), online.ODRF() and plot_ODT_depth().
