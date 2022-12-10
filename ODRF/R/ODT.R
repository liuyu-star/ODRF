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
#' @param weights      : a vector of values which weigh the samples when considering a split (default [])
#' @param nvartosample : the number of (randomly selected) variables to consider at each node (default all)
#' @return ppTree
#'
#' @useDynLib ODRF
#' @import Rcpp
#' @importFrom stats model.frame model.extract model.matrix na.fail
#' 
#' @export
#'
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
#X=X;y=y;method='g-classification';NodeRotateFun="RotMatPPO";paramList=NULL;catLabel=NULL;MaxDepth=Inf;
#MinLeaf=ifelse(method=='regression',5,1);weights=1;numNode=Inf;
#Levels=if(method=='regression')NULL else levels(as.factor(y));Xcat=NULL;
#Xscale=c("Min-max","Quantile","No")[1];FunDir=getwd();TreeRandRotate=FALSE
#Xcat=c(NULL,0)[1]
ODT=function(formula,data,subset=NULL,weights=NULL,na.action=na.fail,method='i-classification',
             NodeRotateFun="RotMatPPO",FunDir=getwd(),paramList=NULL,catLabel=NULL,
             Xcat=0,MaxDepth=Inf,numNode=Inf,MinLeaf=ifelse(method=='regression',5,1),
             Levels=NULL,Xscale=c("Min-max","Quantile","No")[1],TreeRandRotate=FALSE,...)
{
  # address na values.
  if (any(is.na(data))) {
    warning("NA values exist in data matrix")
  }
  
  Call<-match.call()
  indx<-match(c("formula", "data", "subset", "na.action"),names(Call),nomatch=0L)#, "weights"
  if (indx[[1]]==0)
    stop("a 'formula' argument is required")
  if (indx[[2]]==0){
    stop("a 'data' argument is required")
    #data <- environment(formula)
    #data <- data.frame(eval(formula[[2]]),eval(formula[[3]]))
    #varName=colnames(eval(formula[[3]]))
    #colnames(data)=c(as.character(formula[[2]]),varName)
    #Call$formula=
  }
    
  #Call=match.call(expand.dots = FALSE)
  #indx<-match(c("formula", "data", "subset", "weights", "na.action"),names(Call),nomatch=0L)
  temp<-Call[c(1L,indx)]
  temp[[1L]]<-quote(stats::model.frame)
  temp$drop.unused.levels <- TRUE
  temp <- eval(temp, parent.frame())
  Terms<-attr(temp, "terms")
  
  varName=setdiff(colnames(data),as.character(formula[[2]]))
  y <- model.extract(temp, "response")
  X <- model.matrix(Terms, temp)
  int <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  if (int > 0)
    X <- X[, -int, drop = FALSE]
  
  weights=c(weights,paramList$weights)
  if(!is.null(weights))
    X <- X * matrix(weights,length(y),ncol(X))
  colnames(X)=varName
  
  
  if(!NodeRotateFun%in%ls("package:ODRF")){
    source(paste0(FunDir,"/",NodeRotateFun,".R"))
  }
  
  FUN <- match.fun(NodeRotateFun, descend = TRUE)
  method0=strsplit(method,split = "")[[1]][1]
  
  
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
  
  paramList = ODRF:::defaults(paramList,p,catLabel,NodeRotateFun,method,weights)
  
  if(is.infinite(MaxDepth)) {
    numNode=min(numNode,sum(2^(0:ceiling(log2(n/MinLeaf))))) 
  }else{
    MaxDepth=min(MaxDepth,ceiling(log2(n/MinLeaf)),n/MinLeaf-1)
    numNode=min(numNode,sum(2^(0:MaxDepth)))
  }
  
  nodeXIndx = vector("list", numNode+1);
  nodeXIndx[[1]] = 1:n;
  
  if(method!="regression"){
    nodeNumLabel=matrix(0,0,maxLabel)
    colnames(nodeNumLabel)=Levels
    sl=seq(maxLabel)
  }else{
    nodeNumLabel=matrix(0,0,2)
    colnames(nodeNumLabel)=c("prediction","number")
  }

  nodeRotaMat=NULL;
  nodeDepth = rep(1,numNode);
  nodeCutValue = rep(0,numNode);
  nodeCutIndex = nodeCutValue;
  #nodeLabel = nodeCutValue;
  #nodeNumLabel = nodeCutValue;
  childNode = nodeCutValue;
  
  
  #start create pptree nodeFlags
  ##############################################################################
  currentNode = 1;
  freeNode=2
  while(!is.null(nodeXIndx[[currentNode]])){
    
    if((length(unique(y[nodeXIndx[[currentNode]]]))==1)||
       (length(nodeXIndx[[currentNode]])<=(MinLeaf+1))||
       (nodeDepth[currentNode]>=MaxDepth)||
       (freeNode>=numNode)){
      if(method!="regression"){
        leafLabel = table(Levels[c(sl,y[nodeXIndx[[currentNode]]])])-1
        #nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        #nodeNumLabel[currentNode]=max(leafLabel)
        #nodeNumLabel[currentNode,]=leafLabel
      }else{
        #nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
        #nodeNumLabel[currentNode]=length(y[nodeXIndx[[currentNode]]])
        leafLabel=c(mean(y[nodeXIndx[[currentNode]]]),length(y[nodeXIndx[[currentNode]]]))
      }
      nodeNumLabel=rbind(nodeNumLabel,leafLabel)
      
      nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
      nodeXIndx[currentNode]=NA
      currentNode = currentNode+1;
      next
    }
    
    if(!is.null(weights)){
      Wcd = weights[nodeXIndx[[currentNode]]];
    }else{
      Wcd = 1;
    }
    
    ##########################################
    if(NodeRotateFun=="RotMatMake"){
      sparseM <- MakeRotMat(X[nodeXIndx[[currentNode]],], y[nodeXIndx[[currentNode]]], 
                            paramList$RotMatFun, paramList$PPFun, FunDir, paramList)
    }
    
    if(!NodeRotateFun%in%ls("package:ODRF")){
      paramList$x =X[nodeXIndx[[currentNode]],];
      paramList$y = y[nodeXIndx[[currentNode]]];
      sparseM <- do.call(FUN, paramList)
      paramList$x = NULL;paramList$y = NULL;
    }
    
    if(NodeRotateFun%in%c('RotMatRF','RotMatRand')){
      sparseM <- do.call(FUN, paramList)
    }
    
    if(NodeRotateFun=="RotMatPPO"){
      sparseM=RotMatPPO(x=X[nodeXIndx[[currentNode]],],y=y[nodeXIndx[[currentNode]]],numProj=paramList$numProj,
                        dimProj=paramList$dimProj,catLabel = paramList$catLabel,weights=paramList$weights,
                        ppMethod = paramList$ppMethod,method=paramList$method)#"NNet"
    }
    
    if(NodeRotateFun=="PPO"){
      sparseM=ppRF:::PPO(x=X[nodeXIndx[[currentNode]],],y=y[nodeXIndx[[currentNode]]],q = paramList$q,
                         ppMethod = paramList$ppMethod,weight = paramList$weight,r = paramList$r,
                         lambda = paramList$lambda,energy = paramList$energy,cooling = paramList$cooling,
                         TOL = paramList$TOL,maxiter = paramList$maxiter)
      sparseM=sparseM$projMat
    }
    
    numDr=unique(sparseM[,2]);
    rotaX=matrix(0,p,length(numDr));
    for(i in 1:length(numDr)){
      lrows = which(sparseM[,2] == numDr[i]);
      rotaX[sparseM[lrows, 1],i] = sparseM[lrows, 3];
    }
    ###################################################################
    
    rotaX=X[nodeXIndx[[currentNode]],,drop = FALSE] %*% rotaX;
    
    bestCut=ODRF:::best_cut_node(method0,rotaX,y[nodeXIndx[[currentNode]]],Wcd,MinLeaf,maxLabel);
    
    if(bestCut$BestCutVar==-1){
      TF=TRUE
    }else{
      Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
      TF=min(length(Lindex),length(nodeXIndx[[currentNode]])-length(Lindex))<=MinLeaf;
    }
    if(TF){
      if(method!="regression"){
        leafLabel = table(Levels[c(sl,y[nodeXIndx[[currentNode]]])])-1
        #nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        #nodeNumLabel[currentNode]=max(leafLabel)
        #nodeNumLabel[currentNode,]=leafLabel
      }else{
        #nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
        #nodeNumLabel[currentNode]=length(y[nodeXIndx[[currentNode]]])
        leafLabel=c(mean(y[nodeXIndx[[currentNode]]]),length(y[nodeXIndx[[currentNode]]]))
      }
      nodeNumLabel=rbind(nodeNumLabel,leafLabel)
      
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
    nodeCutIndex[currentNode]=bestCut$BestIndex[bestCut$BestCutVar]
    childNode[currentNode]=freeNode;
    
    #Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
    nodeXIndx[[freeNode]] = nodeXIndx[[currentNode]][Lindex];
    nodeXIndx[[freeNode+1]] = nodeXIndx[[currentNode]][setdiff(1:length(nodeXIndx[[currentNode]]),Lindex)];
    nodeDepth[freeNode+c(0,1)]=nodeDepth[currentNode]+1
    
    nodeNumLabel=rbind(nodeNumLabel,0)
    
    nodeXIndx[currentNode]=NA
    freeNode = freeNode+2;
    currentNode = currentNode+1;
  }
  
  colnames(nodeRotaMat)=c("var","node","coef")
  
  ppTree=list()
  ppTree$call=Call
  ppTree$terms=Terms
  ppTree$method=method
  ppTree$Levels=Levels
  ppTree$NodeRotateFun=NodeRotateFun
  ppTree$paramList=paramList
  ppTree$data=list(subset=subset,weights=weights,na.action=na.action,n=n,p=p,varName=varName,
                   Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel,
                   TreeRandRotate=TreeRandRotate,rotdims=rotdims,rotmat=rotmat);
  ppTree$tree=list(FunDir=FunDir,MaxDepth=MaxDepth,MinLeaf=MinLeaf,numNode=numNode);
  ppTree$structure=list(nodeRotaMat = nodeRotaMat,nodeNumLabel = nodeNumLabel,nodeCutValue =nodeCutValue[1:(currentNode-1)],
                        nodeCutIndex =nodeCutIndex[1:(currentNode-1)],childNode = childNode[1:(currentNode-1)],
                        nodeDepth=nodeDepth[1:(currentNode-1)])
  class(ppTree) <- "ODT"
  return(ppTree)
}
