#' Plot method for ODRF objects
#' 
#' Draw the error graph of class \code{ODRF} at different number of trees.
#' 
#' @param Err Object of class \code{\link{ODRF.error}}.
#' @param lty a vector of line types, see \code{\link{par}}.
#' @param main main title of the plot.
#' @param ... arguments to be passed to methods.
#' 
#' @return OOB error and test error, classification error rate for classification or RPE(MSE/mean((ytest-mean(y))^2)) for regression.
#' 
#' @keywords forest
#' 
#' @seealso \code{ODRF} \code{ODRF.error}
#' 
#' @examples
#' library(ODRF)
#' 
#' data(seeds)
#' set.seed(221212)
#' train = sample(1:209,100)
#' train_data = data.frame(seeds[train,])
#' test_data = data.frame(seeds[-train,])
#' 
#' forest = ODRF(varieties_of_wheat~.,train_data,type='i-classification')
#' error=ODRF.error(forest,train_data,test_data)
#' plot(error)
#' 
#' @aliases plot.ODRF.error
#' @rdname plot.ODRF.error
#' @method plot ODRF.error
#' @export
plot.ODRF.error <- function(Err, lty=1, main=paste0("Oblique ",
                      ifelse(Err$type=="regression","Regression","Classification")," Forest"), ...) {
  
  err <- cbind(Err$err.oob,Err$err.test)
  
  ntrees=length(Err$err.oob)
  if(!is.null(Err$err.test)) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:ntrees, err, type = "l",lty=lty, xlab="trees", ylab="Error",col = c("black","red"),main=main)
    legend("topright", legend = c("OOB","Test"),lty = rep(lty,2), col = c("black","red"))#, bty = "n"
  } else {
    matplot(1:ntrees, err, type = "l",lty=lty, xlab="trees", ylab="Error",col = c("black"),main=main)
    legend("topright", legend = c("OOB"),lty = lty, col = c("black"))
  }
  
  return(invisible(err))
}


