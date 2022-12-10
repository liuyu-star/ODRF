#' carPPtree Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param X an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction
#'
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
VarImp=function(forest,X,y){
  if(forest$method!="regression"){
    y <- factor(y,levels = forest$Levels)
  }
  #X=forest$data$X
  #y=forest$data$y
  ntrees=forest$tree$ntrees
  n=length(y);p=ncol(X)
  
  Xcat=forest$data$Xcat
  catLabel=forest$data$catLabel
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
  if(forest$data$Xscale!="No"){
    indp=(sum(numCat)+1):p
    X[,indp]=(X[,indp]-matrix(forest$data$minCol,n,length(indp),byrow = T))/
      matrix(forest$data$maxminCol,n,length(indp),byrow = T)
  }
  
  runOOBErr=function(tree,...){
    class(tree)="ODT"
    oobErrs=rep(0,p+1)
    oobIndex=tree$oobIndex
    
    Xi=X[oobIndex,]
    yi=y[oobIndex]
    yn=length(yi)
    if(forest$method=="regression"){
      e.0=mean((yi-mean(y[-oobIndex]))^2)
    }
    for (j in 1:(p+1)) {
        if(j!=1){
          Xi[,j-1]=Xi[sample(yn),j-1]#+rnorm(length(oobIndex))
        }
        pred <- predict(tree,Xi)
        if(forest$method!="regression"){
          oobErr=mean(pred!=yi);
        }else{
          oobErr=mean((pred-yi)^2)/e.0
        }
        
      if(j==1){
        oobErrs[1]=oobErr
      }else{
        oobErrs[j]=abs(oobErr-oobErrs[1])
      }
    }
    
    return(oobErrs)
  }
  
  oobErrVar=sapply(forest$ppTrees,runOOBErr)[-1,]
  oobErrVar=rowMeans(oobErrVar)
  varimport=cbind(varible=seq(p),increased_error=oobErrVar)
  rownames(varimport)=colnames(X)

  varimport=list(varImp=varimport[order(oobErrVar,decreasing = T),],method=forest$method)
  class(varimport)="VarImp"
  
  return(varimport)
}   
