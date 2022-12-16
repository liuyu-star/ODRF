#' Print the projection pursuit classification tree result
#' @title Print PP.Tree.class result
#' @param x PPtreeclass object
#' @keywords tree
#' @aliases print.ODRF
#' @rdname print.ODRF
#' @method print ODRF
#' @export
#' 
### @method print ppRF
print.ODRF <-function(ppForest, ...) {
  cat("\nCall:\n", deparse(ppForest$call), "\n")
  cat("               Type of oblique decision random forest: ", 
      ifelse(ppForest$type=="regression","regression","classification"), "\n", sep="")
  cat("                                      Number of trees: ", ppForest$forest$ntrees, "\n",sep="")
  #cat("No. of variables tried at each split: ", x$mtrforest, "\n\n", sep="")
  if(ppForest$forest$numOOB==0){
  cat("The number of OOBs is 0")
  }
  if(!is.null(ppForest$oobErr)){
  cat("                           OOB estimate of error rate: ",round(ppForest$oobErr*100,2), "%\n", sep="")
  }
  if(!is.null(ppForest$oobConfusionMat)) {    
  cat("Confusion matrix:\n")
    print(ppForest$oobConfusionMat)
  }
}
