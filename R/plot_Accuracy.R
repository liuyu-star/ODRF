#' plot method for \code{Accuracy} objects
#'
#' Draw the error graph of class \code{ODRF} at different tree sizes.
#'
#' @param x Object of class \code{\link{Accuracy}}.
#' @param lty A vector of line types, see \code{\link{par}}.
#' @param digits Integer indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @param main main title of the plot.
#' @param ... Arguments to be passed to methods.
#'
#' @return OOB error and test error, misclassification rate (MR) for classification or mean square error (MSE) for regression.
#'
#' @keywords plot forest
#'
#' @seealso \code{\link{ODRF}} \code{\link{Accuracy}}
#'
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train <- sample(1:569, 80)
#' train_data <- data.frame(breast_cancer[train, -1])
#' test_data <- data.frame(breast_cancer[-train, -1])
#'
#' forest <- ODRF(diagnosis ~ ., train_data, split = "gini",
#' parallel = FALSE, ntrees = 30)
#' (error <- Accuracy(forest, train_data, test_data))
#' plot(error)
#'
#' @rdname plot.Accuracy
#' @aliases plot.Accuracy
#' @method plot Accuracy
#' @export
plot.Accuracy <- function(x, lty = 1, digits = NULL, main = NULL, ...) {
  if (is.null(main)) {
    main <- paste0("Oblique ", ifelse(x$split == "mse", "Regression", "Classification"), " Forest")
  }
  Err <- x
  err <- cbind(Err$err.oob, Err$err.test)

  minErr <- strsplit(as.character(min(err)), "")[[1]]
  id <- which(minErr == "e")
  if (x$split != "mse") {
    digits <- 0
  } else if (is.null(digits)) {
    if (length(id) > 0) {
      digits <- sum(as.numeric(paste0(minErr[c(id + 2, length(minErr))])) * c(10, 1))
    } else {
      digits <- which(minErr[-seq(which(minErr == "."))] != 0)[2]
    }
  }

  if(digits==0){
    ylab = paste0("Error")
  }else if(digits==2){
    ylab = paste0("Error (%)")
  }else{
    ylab = substitute(paste("Error ("*10^{-dig},")"),list(dig = digits))
  }

  err <- round(err * 10^digits, 2)

  ntrees <- length(Err$err.oob)
  if (!is.null(Err$err.test)) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:ntrees, err, type = "l", lty = lty, xlab = "trees", ylab = ylab, col = c("black", "red"), main = main)
    legend("topright", legend = c("OOB", "Test"), lty = rep(lty, 2), col = c("black", "red"), bty = "n") # , bty = "n"
  } else {
    matplot(1:ntrees, err, type = "l", lty = lty, xlab = "trees", ylab = ylab, col = c("black"), main = main)
    legend("topright", legend = c("OOB"), lty = lty, col = c("black"), bty = "n")
  }
  # axis(1, seq(min(err[,2]),max(err[,2]),length.out = min(6,ntrees)),
  #     round(min(err[,2]),max(err[,2]),length.out = min(6,ntrees),2)*10^digits,
  #     cex.lab = 1.5,cex.axis = 1.25)

  return(invisible(err))
}
