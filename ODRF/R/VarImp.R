#' variable importance of oblique decision random forest
#' 
#' variable importance is computed from permuting OOB data.
#' 
#' @param ppForest an object of class \code{\link{ODRF}}.
#' @param X An n by d numerical matrix (preferably) or data frame is used in the \code{ODRF}.
#' @param y A response vector of length n is used in the \code{ODRF}.
#' 
#' @return A matrix of importance measure, first column is the predictors and second column is Increased error. Misclassification rate (MR) for classification or mean square error (MSE) for regression.
#'
#' @details A note from \code{randomForest}, here are the definitions of the variable importance measures. The first measure is computed from permuting OOB data: For each tree, the prediction error on the out-of-bag portion of the data is recorded.
#' Then the same is done after permuting each predictor variable. The difference between the two are then averaged over all trees, and normalized by the standard deviation of the differences. 
#' If the standard deviation of the differences is equal to 0 for a variable, the division is not done (but the average is almost always equal to 0 in that case).
#' 
#' @seealso \code{\link{ODRF}} \code{\link{plot.VarImp}}
#'
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train = sample(1:569,200)
#' train_data = data.frame(breast_cancer[train,-1])
#' test_data = data.frame(breast_cancer[-train,-1])
#'
#' forest = ODRF(diagnosis~.,train_data,type='i-classification',parallel=FALSE)
#' (varimp=VarImp(forest,train_data[,-1],train_data[,1]))
#' 
#' @export
VarImp=function(ppForest,X,y){
  #vars=all.vars(ppForest$terms)
  # address na values.
  #if (any(is.na(data))) {
  #  data=ppForest$data$na.action(data.frame(data))
  #  warning("NA values exist in data matrix")
  #}
  #y= data[,setdiff(colnames(data),vars[-1])]
  #X= data[,vars[-1]]
  X=as.matrix(X)
  
  if(!is.null(ppForest$data$subset))
    X=X[ppForest$data$subset,]
  #if(!is.null(weights))
  #  X <- X * matrix(weights0,length(y),ncol(X))
  #weights=weights0
  
  
  if(ppForest$type!="regression"){
    y <- factor(y,levels = ppForest$Levels)
  }
  #X=ppForest$data$X
  #y=ppForest$data$y
  ntrees=ppForest$tree$ntrees
  n=length(y);p=ncol(X)
  
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
    class(tree)="ODT"
    oobErrs=rep(0,p+1)
    oobIndex=tree$oobIndex
    
    Xi=X[oobIndex,]
    yi=y[oobIndex]
    yn=length(yi)
    #if(ppForest$type=="regression"){
    #  e.0=mean((yi-mean(y[-oobIndex]))^2)
    #}
    for (j in 1:(p+1)) {
        if(j!=1){
          Xi[,j-1]=Xi[sample(yn),j-1]#+rnorm(length(oobIndex))
        }
        pred <- predict(tree,Xi)
        if(ppForest$type!="regression"){
          oobErr=mean(pred!=yi);
        }else{
          oobErr=mean((pred-yi)^2)#/e.0
        }
        
      if(j==1){
        oobErrs[1]=oobErr
      }else{
        oobErrs[j]=abs(oobErr-oobErrs[1])
      }
    }
    
    return(oobErrs)
  }
  
  oobErrVar=sapply(ppForest$ppTrees,runOOBErr)[-1,]
  oobErrVar=rowMeans(oobErrVar)
  varimport=cbind(varible=seq(p),increased_error=oobErrVar)
  rownames(varimport)=colnames(X)

  varimport=list(varImp=varimport[order(oobErrVar,decreasing = T),],type=ppForest$type)
  class(varimport)="VarImp"
  
  return(varimport)
}   
