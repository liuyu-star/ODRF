#' Plot method for ODRF.accuracy objects
#' 
#' Draw the error graph of class \code{ODRF} at different number of trees.
#' 
#' @param Err Object of class \code{\link{ODRF_accuracy}}.
#' @param lty a vector of line types, see \code{\link{par}}.
#' @param digits integer indicating the number of decimal places (round) or significant digits (signif) to be used. 
#' @param main main title of the plot.
#' @param ... arguments to be passed to methods.
#' 
#' @return OOB error and test error, misclassification rate (MR) for classification or mean square error (MSE) for regression.
#' 
#' @keywords forest
#' 
#' @seealso \code{\link{ODRF}} \code{\link{ODRF_accuracy}}
#' 
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train = sample(1:569,200)
#' train_data = data.frame(breast_cancer[train,-1])
#' test_data = data.frame(breast_cancer[-train,-1])
#'
#' forest = ODRF(diagnosis~.,train_data,type='i-classification',parallel=FALSE)
#' (error=ODRF_accuracy(forest,train_data,test_data))
#' plot(error)
#' 
#' @rdname plot.ODRF_accuracy
#' @aliases plot.ODRF_accuracy
#' @method plot ODRF_accuracy
#' @export
plot.ODRF_accuracy <- function(Err, lty=1, digits=NULL, main=paste0("Oblique ",
                      ifelse(Err$type=="regression","Regression","Classification")," Forest"), ...) {
  
  err <- cbind(Err$err.oob,Err$err.test)
  
  minErr=strsplit(as.character(min(err)),"")[[1]]
  id=which(minErr=="e")
  if(Err$type!="regression"){
    digits=0
  }else if(is.null(digits)){
    if(length(id)>0){
      digits=sum(as.numeric(paste0(minErr[c(id+2,length(minErr))]))*c(10,1))
    }else{
      digits=which(minErr[-seq(which(minErr=="."))]!=0)[2]
    }
  }
  err=round(err*10^digits, 2)
  
  ntrees=length(Err$err.oob)
  if(!is.null(Err$err.test)) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:ntrees, err, type = "l",lty=lty, xlab="trees", ylab=paste0("Error (*",10^-digits,")"),col = c("black","red"),main=main)
    legend("topright", legend = c("OOB","Test"),lty = rep(lty,2), col = c("black","red"),bty="n")#, bty = "n"
  } else {
    matplot(1:ntrees, err, type = "l",lty=lty, xlab="trees", ylab=paste0("Error (*",10^-digits,")"),col = c("black"),main=main)
    legend("topright", legend = c("OOB"),lty = lty, col = c("black"),bty="n")
  }
  #axis(1, seq(min(err[,2]),max(err[,2]),length.out = min(6,ntrees)),
  #     round(min(err[,2]),max(err[,2]),length.out = min(6,ntrees),2)*10^digits,
  #     cex.lab = 1.5,cex.axis = 1.25)
  
  return(invisible(err))
}


