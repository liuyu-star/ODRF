#' RerF Tree Generator
#' 
#' Grows a CARTree using X(samplesXfeatures) and one of the following criteria which can be set via the parameter 'method'
#'
#' @param       'g-classification' : gini impurity index (classification)
#' @param       'i-classification' : information gain (classification, default)
#' @param       'regression' : squared error (regression)
#' 
#' Other parameters that can be set are:
#' 
#' @param paramList parameters in a named list to be used by FUN. If left unchanged,
#' @param catMap a list specifying which columns in X correspond to the same one-of-K encoded feature. Each element of catMap is a numeric vector specifying the K column indices of X corresponding to the same categorical feature after one-of-K encoding. All one-of-K encoded features in X must come after the numeric features. The K encoded columns corresponding to the same categorical feature must be placed contiguously within X. The reason for specifying catMap is to adjust for the fact that one-of-K encoding cateogorical features results in a dilution of numeric features, since a single categorical feature is expanded to K binary features. If catMap = NULL, then RerF assumes all features are numeric (i.e. none of the features have been one-of-K encoded).
#' @param FUN :FUN=RandMatPPR.
#' @param minparent    : the minimum amount of samples in an impure node for it to be considered for splitting (default 2)
#' @param MinLeaf      : the minimum amount of samples in a leaf (default 1)
#' @param Xweights      : a vector of values which weigh the samples when considering a split (default [])
#' @param nvartosample : the number of (randomly selected) variables to consider at each node (default all)
#' @return ppTree
#'
#' @useDynLib ppRF
#' @import Rcpp
#' 
#' @export
#'
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
#X=X;y=y;method='g-classification';NodeRotateFun="RandMatPPR";paramList=NULL;catLabel=NULL;MaxDepth=Inf;
#MinLeaf=ifelse(method=='regression',5,1);Xweights=1;numNode=Inf;
#Levels=if(method=='regression')NULL else levels(as.factor(y));Xcat=NULL;
#Xscale=c("Min-max","Quantile","No")[1];FunDir=getwd();TreeRandRotate=FALSE
#Xcat=c(NULL,0)[1]
ppCART=function(X,y,method='g-classification',NodeRotateFun="RandMatPPR",FunDir=getwd(),paramList=NULL,catLabel=NULL,
                 Xcat=0,MaxDepth=Inf,numNode=Inf,MinLeaf=ifelse(method=='regression',5,1),Xweights=1,
                 Levels=NULL,Xscale=c("Min-max","Quantile","No")[1],TreeRandRotate=FALSE)
{
  #ppTree=ppCART(X,y,method,NodeRotateFun,FunDir,paramList,catLabel,Xcat,MaxDepth,
  #       numNode,MinLeaf,weights,Levels,Xscale,TreeRandRotate)
  
  #ppRFOtherFun=c("RandMatBinary","RandMatContinuous","RandMatCustom","RandMatFRC","RandMatFRCN",
  #          "RandMatIden", "RandMatImageControl", "RandMatImagePatch" , "RandMatPoisson")
  #FUNName=FUN
  if(!NodeRotateFun%in%ls("package:ppRF")){
    source(paste0(FunDir,"/",NodeRotateFun,".R"))
  }
  
  FUN <- match.fun(NodeRotateFun, descend = TRUE)
  method0=strsplit(method,split = "")[[1]][1]
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  n = length(y);
  p = ncol(X);
  
  if(method!="regression"){
    if(is.null(Levels)){
      Levels=levels(as.factor(y))
    }
    if (!is.integer(y)) {
      y <- as.integer(as.factor(y));
    }
    maxLabel = length(Levels); 
  }else{
    y=c(y)
    maxLabel= 0;
  }
  
  
  if(is.null(Xcat)){
    Xcat=which(apply(X, 2, function(x){(length(unique(x))<10)&(n>20)}))
  }
  
  numCat=0
  if((sum(Xcat)>0)&(is.null(catLabel))){
    numCat <- apply(X[,Xcat,drop = FALSE], 2, function(x) length(unique(x)))
    X1 <- matrix(0, nrow = n, ncol = sum(numCat)) # initialize training data matrix X
    catLabel <- vector("list", length(Xcat))
    names(catLabel)<- colnames(X)[Xcat]
    col.idx <- 0L
    # one-of-K encode each categorical feature and store in X
    for (j in 1:length(Xcat)) {
      catMap <- (col.idx + 1L):(col.idx + numCat[j])
      # convert categorical feature to K dummy variables
      catLabel[[j]]=levels(as.factor(X[,Xcat[j]]))
      X1[, catMap] <- (matrix(X[,Xcat[j]],n,numCat[j])==matrix(catLabel[[j]],n,numCat[j],byrow = TRUE))+0
      col.idx <- col.idx + numCat[j]
    }
    X=cbind(X1,X[,-Xcat]) 
    rm(X1)
    p = ncol(X);
  }
  
  
  #Variable scaling.
  minCol = NULL;maxminCol = NULL;
  if(Xscale!="No"){
    indp=(sum(numCat)+1):p
    if(Xscale=="Min-max"){
      minCol=apply(X[,indp],2,min)
      maxminCol=apply(X[,indp],2,function(x){max(x)-min(x)})
    }
    if(Xscale=="Quantile"){
      minCol=apply(X[,indp],2,quantile,0.05)
      maxminCol=apply(X[,indp],2,function(x){quantile(x,0.95)-quantile(x,0.05)})
    }
    X[,indp]=(X[,indp]-matrix(minCol,n,length(indp),byrow = T))/matrix(maxminCol,n,length(indp),byrow = T)
  }
  
  # rotate the data?
  rotdims=NULL;rotmat=NULL
  if (TreeRandRotate) {
    # if p > 1000 then rotate only a random subset of 1000 of the dimensions
    if (p > 1000L) {
      rotmat <- RandRot(1000L)
      rotdims <- sample.int(p, 1000L)
      X[, rotdims] <- X[, rotdims] %*% rotmat
    } else {
      rotdims=1:p
      rotmat <- RandRot(p)
      X <- X %*% rotmat
    }
  }
  
  
  d=paramList$d
  paramList = ppRF:::defaults(paramList,x=X,y=y,method=method,catLabel=catLabel)
  #if(!FUNName%in%ppRFOtherFun){paramList$d=NULL}
  paramList$x = NULL;paramList$y = NULL;paramList$d=d
  
  #L = 2*ceiling(n);
  #L = 2*ceiling(n/MinLeaf) - 1;
  # L = 2*ceiling(n/MinLeaf) + 1;
  # L = factorial(ceiling(n/MinLeaf));
  # L = 2^ceiling(n/MinLeaf);
  #MaxDepth=min(ceiling(log2(n/MinLeaf)),MaxDepth)
  #MaxDepth=min(MaxDepth,n/MinLeaf-1) 
  #MaxDepth=n/MinLeaf-1
  #numNode=2*MaxDepth+1
  if(is.infinite(MaxDepth)) {
    numNode=min(numNode,sum(2^(0:ceiling(log2(n/MinLeaf))))) 
    index=c()
  }else{
    MaxDepth=min(MaxDepth,ceiling(log2(n/MinLeaf)),n/MinLeaf-1)
    numNode=min(numNode,sum(2^(0:MaxDepth)))
    index=2:numNode;
    totalN=0;
    totalNode = rep(0,numNode);
    totalNode[1] = 1;
  }
  
  nodeXIndx = vector("list", numNode+1);
  nodeXIndx[[1]] = 1:n;
  nodeRotaMat=NULL;
  nodeCutValue = rep(0,numNode);
  nodeLabel = nodeCutValue;
  childNode = nodeCutValue;
  
  
  #start create pptree nodeFlags
  ##############################################################################
  currentNode = 1;
  freeNode=2
  while(!is.null(nodeXIndx[[currentNode]])){
    if(!is.infinite(MaxDepth)){
      repeat{
        if((totalNode[currentNode]!=-1)||(length(index)==0)){break}
        totalNode[index[1:2]] = -1
        totalNode=totalNode[-currentNode]
        index=index[-c(1,2,numNode-totalN)]-1
        totalN=totalN+1
      }
    }
    
    if((length(unique(y[nodeXIndx[[currentNode]]]))==1)||
       (length(nodeXIndx[[currentNode]])<=(MinLeaf+1))||
       (!is.infinite(MaxDepth)&(length(index)==0))||
       (freeNode>=numNode)){
      if(method!="regression"){
        #leafLabel = which.max(table(factor(y[nodeXIndx[[currentNode]]],levels=Levels)));
        leafLabel = table(Levels[y[nodeXIndx[[currentNode]]]])
        nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
      }else{
        nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
      }
      if(length(index)>0){
        totalNode[index[c(1,2)]] = -1;
        index=index[-(1:2)]
      }
      
      nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
      nodeXIndx[currentNode]=NA
      currentNode = currentNode+1;
      next
    }
    
    
    if(length(Xweights)>1){
      Wcd = Xweights[nodeXIndx[[currentNode]]];
    }else{
      Wcd = 1;
    }
    
    
    paramList = ppRF:::defaults(paramList,x=X[nodeXIndx[[currentNode]],],y=y[nodeXIndx[[currentNode]]],method=method,catLabel=catLabel)
    sparseM <- do.call(FUN, paramList)
    paramList$x = NULL;paramList$y = NULL;paramList$d=d
    
    numDr=unique(sparseM[,2]);
    rotaX=matrix(0,length(numDr),p);
    for(i in 1:length(numDr)){
      lrows = which(sparseM[,2] == numDr[i]);
      rotaX[i,sparseM[lrows, 1]] = sparseM[lrows, 3];
    }
    rotaX=X[nodeXIndx[[currentNode]],,drop = FALSE] %*% t(rotaX);
    
    bestCut=ppRF:::best_cut_node(method0,rotaX,y[nodeXIndx[[currentNode]]],Wcd,MinLeaf,maxLabel);
    
    if(bestCut$BestCutVar==-1){
      TF=TRUE
    }else{
      Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
      TF=min(length(Lindex),length(nodeXIndx[[currentNode]])-length(Lindex))<=MinLeaf;
    }
    if(TF){
      if(method!="regression"){
        #leafLabel = which.max(table(factor(y[nodeXIndx[[currentNode]]],levels=Levels)));
        leafLabel = table(Levels[y[nodeXIndx[[currentNode]]]])
        nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
      }else{
        nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
      }
      if(length(index)>0){
        totalNode[index[c(1,2)]] = -1;
        index=index[-(1:2)]
      }
      
      nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
      nodeXIndx[currentNode]=NA
      currentNode = currentNode+1;
      next
    }
    
    sparseM=sparseM[sparseM[,2]==numDr[bestCut$BestCutVar],,drop = FALSE]
    sparseM[,2]=currentNode
    nodeRotaMat=rbind(nodeRotaMat,sparseM)
    rm(rotaX)
    rm(sparseM)
    
    nodeCutValue[currentNode] = bestCut$BestCutVal;
    
    #Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
    nodeXIndx[[freeNode]] = nodeXIndx[[currentNode]][Lindex];
    nodeXIndx[[freeNode+1]] = nodeXIndx[[currentNode]][setdiff(1:length(nodeXIndx[[currentNode]]),Lindex)];
    
    childNode[currentNode]=freeNode;
    if(!is.infinite(MaxDepth)){
      totalNode[index[c(1,2)]] = freeNode:(freeNode + 1);
      index=index[-(1:2)];
    }
    
    nodeXIndx[currentNode]=NA
    freeNode = freeNode+2;
    currentNode = currentNode+1;
  }
  
  colnames(nodeRotaMat)=c("var","node","coef")
  
  ppTree=list()
  ppTree$method=method
  ppTree$Levels=Levels
  ppTree$NodeRotateFun=NodeRotateFun
  ppTree$data=list(n=n,p=p,Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel,
                   TreeRandRotate=TreeRandRotate,rotdims=rotdims,rotmat=rotmat);
  ppTree$tree=list(FunDir=FunDir,paramList=paramList,MaxDepth=MaxDepth,MinLeaf=MinLeaf,numNode=numNode,Xweights=Xweights);
  ppTree$stuction=list(nodeRotaMat = nodeRotaMat,nodeCutValue =nodeCutValue[1:(currentNode-1)],
                       childNode = childNode[1:(currentNode-1)],nodeLabel = nodeLabel[1:(currentNode-1)])
  return(ppTree)
}