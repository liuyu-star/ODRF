## Responses to Comments

We are glad to receive your comments, and we have fixed every one of them according to your request, as follows.

1. Please always explain all acronyms in the description text. -> CART

If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) 
authors (year) 
authors (year, ISBN:...)
or if those are not available: 
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking. (If you want to add a title as well please put it in 
quotes: "Title")

* **Response.** We explained CART and Random Forest in the description text.


2. We see: The Date field is over a month old. Please change this to a more recent date.

* **Response.** We changed the Date field to a more recent date.


3. "Using foo:::f instead of foo::f allows access to unexported objects. 
This is generally not recommended, as the semantics of unexported 
objects may be changed by the package author in routine maintenance."
Used ::: in documentation:
      man/defaults.Rd:
         (paramList <- ODRF:::defaults(paramList, split = "entropy"))
      man/RandRot.Rd:
         (ODRF:::RandRot(100))
Please omit one colon.

* **Response.** We have exported the functions RandRot() and defaults(), and already replaced ODRF:::f with f.


4. Please make sure that you do not change the user's options, par or 
working directory. If you really have to do so within functions, please 
ensure with an *immediate* call of on.exit() that the settings are reset 
when the function is exited.
e.g.:
...
oldpar <- par(no.readonly = TRUE) # code line i
on.exit(par(oldpar)) # code line i + 1
...
par(mfrow=c(2,2)) # somewhere after
...

e.g.: in most of your plot functions.
If you're not familiar with the function, please check ?on.exit. This 
function makes it possible to restore options before exiting a function 
even if the function breaks. Therefore it needs to be called immediately 
after the option change within a function.

* **Response.** We removed par from plot.VarImp() and added on.exit to plot.prune.ODT(). so we made sure not to change the user's options, par or working directory.


5. Please do not set a seed to a specific number within a function. -> 
R/plot_ODT_depth.R; R/prune_ODRF.R

* **Response.** We removed the random seed number in functions ODRF(), poune.ODRF(), online.ODRF() and plot_ODT_depth().


## R CMD check results

0 errors | 0 warnings | 0 note

* This is the second version of the package. It currently has no ERRORs, WARNINGs, or NOTEs from devtools::check().

## Reverse dependencies

* We have run R CMD check on the downstream dependency and no problems were found related to this package.


---
Thanks!
Yu Liu and Yingcun Xia
