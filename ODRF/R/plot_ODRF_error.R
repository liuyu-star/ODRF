#' projection pursuit classification tree plot
#' 
#' Draw projection pursuit classification tree with tree structure. It is 
#' modified from a function in party library.
#' @title PPtree plot
#' @param x PPtreeclass object
#' @param font.size font size of plot
#' @param width.size size of eclipse in each node.
#' @param main main title
#' @param sub sub title
#' @param ... arguments to be passed to methods
#' @references Lee, EK(2017) 
#' PPtreeViz: An R Package for Visualizing Projection Pursuit Classification 
#' Trees, Journal of Statistical Software <doi:10.18637/jss.v083.i08>
#' @keywords tree
#' 
#' @aliases plot.ODRF.error
#' @rdname plot.ODRF.error
#' @method plot ODRF.error
#' @export
#' 
#' @examples
#' data(iris)
#' Tree.result <- PPTreeclass(Species~., data = iris,"LDA")
#' Tree.result
#' plot(Tree.result,xjust=3)
plot.ODRF.error <- function(Err, lty=1, main=paste0("Oblique ",
                      ifelse(Err$method=="regression","Regression","Classification")," Forest"), ...) {
  
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


