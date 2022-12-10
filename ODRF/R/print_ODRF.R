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
print.ODRF <-function(forest, ...) {
  cat("\nCall:\n", deparse(forest$call), "\n")
  cat("               Type of oblique decision random forest: ", 
      ifelse(forest$method=="regression","regression","classification"), "\n", sep="")
  cat("                                      Number of trees: ", forest$forest$ntrees, "\n",sep="")
  #cat("No. of variables tried at each split: ", x$mtrforest, "\n\n", sep="")
  if(forest$forest$numOOB==0){
  cat("The number of OOBs is 0")
  }
  if(!is.null(forest$oobErr)){
  cat("                           OOB estimate of error rate: ",round(forest$oobErr*100,2), "%\n", sep="")
  }
  if(!is.null(forest$oobConfusionMat)) {    
  cat("Confusion matrix:\n")
    print(forest$oobConfusionMat)
  }
}
