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
OOBErr=function(X,y,ppForest){
  #OOBindx=unique(unlist(lapply(ppForest$ppTrees,function(trees)trees$OOBindx)))
  #OOBindx= Reduce(intersect, lapply(ppForest$ppTrees,function(trees)trees$OOBindx))
  #OOBindx=seq(n)
  #for (it in 1:ntrees) {
  #  OOBindx=intersect(OOBindx,ppForest$ppTrees[[it]]$OOBindx)
  # }
  #unique(unlist(sapply(ppForest$ppTrees,runOOBErr)))
  #ppForest$forest$numOOB=length(OOBindx)
  
  X <- as.matrix(X)
  y=c(y)
  p=ncol(X)
  n=nrow(X)
  nC=length(ppForest$Levels)
  ntrees=length(ppForest$ppTrees);
  Levels=ppForest$Levels
  method=ppForest$method
  
  pred <- predict_ppRF(X,ppForest,weight = FALSE)$prediction
  XConfusionMat=table(pred,Levels[y])
  class_error=(rowSums(XConfusionMat)-diag(XConfusionMat))/rowSums(XConfusionMat)
  XConfusionMat=cbind(XConfusionMat,class_error)
  
  Xcat=ppForest$data$Xcat
  catLabel=ppForest$data$catLabel
  numCat=0
  if(sum(Xcat)>0){
    xj=1
    X1 <- matrix(0, nrow = n, ncol = length(unlist(catLabel))) # initialize training data matrix X
    # one-of-K encode each categorical feature and store in X
    for (j in 1:length(Xcat)) {
      catMap= which(catLabel[[j]]%in%unique(X[,Xcat[j]]))
      indC=catLabel[[j]][catMap]
      Xj <- (matrix(X[,Xcat[j]],n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
      
      if(length(indC)>length(catLabel[[j]])){
        Xj=Xj[,1:length(catLabel[[j]])]
      }
      
      xj1=xj+length(catLabel[[j]])
      X1[,(xj:(xj1-1))[catMap]]=Xj
      xj=xj1
    }
    
    X=cbind(X1,X[,-Xcat]) 
    p=ncol(X)
    numCat=length(unlist(catLabel))
    rm(X1)
    rm(Xj)
  }
  
  #Variable scaling.
  if(ppForest$data$Xscale!="No"){
    indp=(sum(numCat)+1):p
    X[,indp]=(X[,indp]-matrix(ppForest$data$minCol,n,length(indp),byrow = T))/
      matrix(ppForest$data$maxminCol,n,length(indp),byrow = T)
  }
  
  runOOBErr=function(tree,...){
    OOBindx=tree$OOBindx
    if(length(OOBindx)>0){
      pred <- predict_ppCART(X[OOBindx,],tree)
      if(method!="regression"){
        oobErr=mean(pred!=Levels[y[OOBindx]]);
      }else{
        oobErr=mean((pred-y[OOBindx])^2);
      }
    }else{
      oobErr=1
    }
  }
  
  oobErr=sapply(ppForest$ppTrees,runOOBErr)

  return(list(oobErr=oobErr,XConfusionMat=XConfusionMat))
}
