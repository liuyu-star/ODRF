#' print ODRF
#'
#' Print contents of ODRF object.
#'
#' @param x An object of class \code{\link{ODRF}}.
#' @param ... Arguments to be passed to methods.
#'
#' @return OOB error, misclassification rate (MR) for classification or mean square error (MSE) for regression.
#'
#' @seealso \code{\link{ODRF}}
#'
#' @examples
#' data(iris)
#' forest <- ODRF(Species ~ ., data = iris, parallel = FALSE, ntrees = 50)
#' forest
#'
#' @keywords forest print
#' @rdname print.ODRF
#' @aliases print.ODRF
#' @method print ODRF
#' @export
print.ODRF <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("               Type of oblique decision random forest: ",
    ifelse(x$split %in% c("gini","entropy"), "classification", "regression"), "\n",
    sep = ""
  )
  cat("                                      Number of trees: ", x$forest$ntrees, "\n", sep = "")
  # cat("No. of variables tried at each split: ", x$mtrforest, "\n\n", sep="")
  if (x$forest$ratOOB == 0) {
    cat("The number of OOBs is 0")
  }
  if (!is.null(x$oobErr)) {
    cat("                           OOB estimate of error rate: ", round(x$oobErr * 100, 2), "%\n", sep = "")
  }
  if (!is.null(x$oobConfusionMat)) {
    cat("Confusion matrix:\n")
    print(x$oobConfusionMat)
  }
}
