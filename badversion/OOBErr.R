#' carPPtree Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param X an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction
#'
#' @useDynLib ppRF
#' @import Rcpp
#' 
#' @export
#' 
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#'
#' ### Train RerF on one-of-K encoded categorical data ###
#' df1 <- as.data.frame(Titanic)
#' nc <- ncol(df1)
#' df2 <- df1[NULL, -nc]
#' for (i in which(df1$Freq != 0L)) {
#'   df2 <- rbind(df2, df1[rep(i, df1$Freq[i]), -nc])
#' }
## @rdname predict#'  #' @aliases predict_ppCART #' @method predict ppCART
OOBErr=function(ppForest,X,y){
  #OOBindx=unique(unlist(lapply(ppForest$ppTrees,function(trees)trees$OOBindx)))
  #OOBindx= Reduce(intersect, lapply(ppForest$ppTrees,function(trees)trees$OOBindx))
  #OOBindx=seq(n)
  #for (it in 1:ntrees) {
  #  OOBindx=intersect(OOBindx,ppForest$ppTrees[[it]]$OOBindx)
  # }
  #unique(unlist(sapply(ppForest$ppTrees,runOOBErr)))
  #ppForest$forest$numOOB=length(OOBindx)
  #nOOB=cumsum(sapply(OOBindx, length))
  
  X=as.matrix(X)
  y=as.vector(y)
  
  oobIndex=lapply(ppForest$ppTrees,function(trees)trees$oobIndex)
  indx=sort(Reduce(union,oobIndex))
  treeOOB=lapply(indx,function(i)which(sapply(oobIndex,function(treeOOB)i%in%treeOOB)))

  runOOBErr=function(idx,...){
    forest=ppForest
    forest$ppTrees=forest$ppTrees[treeOOB[[idx]]]
    pred <- predict_ppRF(forest,X[indx[idx],,drop = FALSE],weight = FALSE)$prediction
    return(pred)
  }
  
  pred=sapply(seq(length(indx)),runOOBErr)
  if(ppForest$method!="regression"){
    oobErr=mean(pred!=y[indx]);#ppForest$Levels[
  }else{
    oobErr=mean((pred-y[indx])^2);
  }
  
  XConfusionMat=table(pred,y)
  class_error=(rowSums(XConfusionMat)-diag(XConfusionMat))/rowSums(XConfusionMat)
  XConfusionMat=cbind(XConfusionMat,class_error)
  
  return(list(oobErr=oobErr,ConfusionMat=XConfusionMat))
}
