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
#' @useDynLib ODRF
#' @import Rcpp
#' 
#' @aliases prune.ODT
#' @rdname prune.ODT
#' @method prune ODT
#' @export
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
#ppTree=ppTree0
#Xnew=X1
#ynew=y1
#Levels=ppTree$Levels
#method=ppTree$method
prune.ODT=function(ppTree,data,weights=NULL,MaxDepth=NULL)
{
  structure=ppTree$structure
  if(!is.null(MaxDepth)){
    MaxDepth=min(MaxDepth,max(structure$nodeDepth))
  }else{
    MaxDepth=1
  }
  numNode=length(structure$nodeCutValue)
  
  # address na values.
  if (any(is.na(data))) {
    data=ppTree$data$na.action(data.frame(data))
    warning("NA values exist in data matrix")
  }
  ynew= data[,all.vars(ppTree$terms)[1]]
  Xnew= data[,all.vars(ppTree$terms)[-1]]
  Xnew=as.matrix(Xnew) 
  
  if(!is.null(ppTree$data$subset))
    Xnew=Xnew[ppTree$data$subset,]
  weights0=c(ppTree$data$weights,ppTree$paramList$weights)
  if(!is.null(weights0))
    Xnew <- Xnew * matrix(weights,length(y),ncol(Xnew))
  
  p=ncol(Xnew)
  n=nrow(Xnew)
  
  Xcat=ppTree$data$Xcat
  catLabel=ppTree$data$catLabel
  numCat=0
  if(sum(Xcat)>0){
    xj=1
    Xnew1 <- matrix(0, nrow = n, ncol = length(unlist(catLabel))) # initialize training data matrix X
    # one-of-K encode each categorical feature and store in X
    for (j in 1:length(Xcat)) {
      catMap= which(catLabel[[j]]%in%unique(Xnew[,Xcat[j]]))
      indC=catLabel[[j]][catMap]
      Xnewj <- (matrix(Xnew[,Xcat[j]],n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
      
      if(length(indC)>length(catLabel[[j]])){
        Xnewj=Xnewj[,1:length(catLabel[[j]])]
      }
      
      xj1=xj+length(catLabel[[j]])
      Xnew1[,(xj:(xj1-1))[catMap]]=Xnewj
      xj=xj1
    }
    
    Xnew=cbind(Xnew1,Xnew[,-Xcat]) 
    p=ncol(Xnew)
    numCat=length(unlist(catLabel))
    rm(Xnew1)
    rm(Xnewj)
  }
  
  #Variable scaling.
  if(ppTree$data$Xscale!="No"){
    indp=(numCat+1):p
    Xnew[,indp]=(Xnew[,indp]-matrix(ppTree$data$minCol,n,length(indp),byrow = T))/
      matrix(ppTree$data$maxminCol,n,length(indp),byrow = T)
  }
  
  if (ppTree$data$TreeRandRotate) {
    Xnew[,ppTree$data$rotdims] <- Xnew[,ppTree$data$rotdims,drop = FALSE] %*% ppTree$data$rotmat
  }
  
  
  if(ppTree$method!="regression"){
    #nodeLabel=rep(0,NROW(structure$nodeNumLabel))
    #leafid=which(rowSums(structure$nodeNumLabel)!=0)
    #leafLabel=structure$nodeNumLabel[leafid,,drop = FALSE]
    nodeLabel=colnames(structure$nodeNumLabel)[max.col(structure$nodeNumLabel)] ## "random"
    nodeLabel[which(rowSums(structure$nodeNumLabel)==0)]=0
  }else{
    nodeLabel=as.character(structure$nodeNumLabel[,1])
  }
  
  if(all(structure$nodeCutValue==0)){
    prediction=rep(nodeLabel,n)
  }else{
    prediction = .Call('_ODRF_predict_ODT', PACKAGE = 'ODRF',Xnew, structure$nodeRotaMat, 
                       structure$nodeCutValue, structure$childNode,nodeLabel)$prediction
  }
  
  if(ppTree$method!="regression"){
    err0 <- mean(prediction != ynew)
  }else{
    err0 <- mean((as.numeric(prediction)-ynew)^2)
    #structure$nodeLabel=as.numeric(structure$nodeLabel)
  }
  
  
  #start create pptree.
  ##############################################################################
  cutNode=which(structure$nodeCutValue!=0)
  ncut=length(cutNode)
  pruneError=matrix(0,ncut+1,4)
  colnames(pruneError)=c("cutNode","numNode","depth","error")
  #pruneError[1,]=c(1,1,1,err0)
  pruneError[ncut+1,]=c(ncut,numNode,max(structure$nodeDepth),err0)
  
  currentNode=cutNode[ncut]
  while((ncut>=1)&(structure$nodeDepth[currentNode]>=MaxDepth)){
    
    nodeRotaMat=structure$nodeRotaMat
    childNode=structure$childNode
    nodeCutValue=structure$nodeCutValue
    #nodeLabel=structure$nodeLabel
    nodeNumLabel=structure$nodeNumLabel
    
    freeNode=childNode[currentNode];
    idx=c(freeNode,freeNode+1)
    
    ni=1;#id=c()
    while(ni<=length(idx)){ 
      #ni=which(node==idx)
      node=idx[ni]
      if(nodeCutValue[node]!=0){
        cn=childNode[node]
        #idx=c(idx,cn,cn+1)
        #idx=setdiff(idx,node)
        idx=c(idx,cn,cn+1)
      }
      ni=ni+1
    }
    
    
    for (ni in seq(length(idx),1,-1))  {
      node=idx[ni]
      if(ni%%2==0){
        nodeRotaMat=nodeRotaMat[-which(nodeRotaMat[,2]%in%c(node-1,node)),,drop = FALSE]
        id=which(nodeRotaMat[,2]>node)
        nodeRotaMat[id,2]=nodeRotaMat[id,2]-2
      }
      
      if(nodeCutValue[node]!=0){
        id=seq(node)
        childNode[-id]=(childNode[-id]!=0)*(childNode[-id]-2)
      }
    }
    
    id=which(nodeRotaMat[,2]==currentNode)
    nodeRotaMat[id[1],]=c(0,currentNode,0)
    if(length(id[-1])>0){nodeRotaMat=nodeRotaMat[-id[-1],,drop = FALSE]}
    #id=min(which(nodeRotaMat[,2]==idx[length(idx)]))
    #nodeRotaMat[-seq(max(id)),2]=nodeRotaMat[-seq(max(id)),2]-2
    #nodeRotaMat=nodeRotaMat[id,,drop = FALSE]
    
    childNode[currentNode]=0;
    childNode[-seq(currentNode)]=(childNode[-seq(currentNode)]!=0)*(childNode[-seq(currentNode)]-2)
    childNode=childNode[-idx];
    
    id=idx[nodeCutValue[idx]==0]
    #id=idx[!idx%in%cutNode]
    #nodeLabel[currentNode]=ifelse(ppTree$method!="regression",nodeLabel[idx][which.max(structure$nodeNumLabel[idx])],
    #                              structure$nodeNumLabel[idx]*nodeLabel[idx]/sum(structure$nodeNumLabel[idx]))
    if(ppTree$method!="regression"){
      #nnl=rep(nodeLabel[id],nodeNumLabel[id]) 
      #nnl=table(nnl)
      #nodeLabel[currentNode]=names(nnl)[which.max(nnl)]
      #nodeNumLabel[currentNode]=nnl[which.max(nnl)]
      nodeNumLabel[currentNode,]=colSums(nodeNumLabel[id,,drop = FALSE])
    }else{
      #nodeLabel[currentNode]=nodeNumLabel[id]*nodeLabel[id]/sum(nodeNumLabel[id])
      #nodeNumLabel[currentNode]=sum(nodeNumLabel[id])
      nodeNumLabel[currentNode,1]=sum(nodeNumLabel[id,1]*nodeNumLabel[id,2])/sum(nodeNumLabel[id,2])
      nodeNumLabel[currentNode,2]=sum(nodeNumLabel[id,2])
    }
    #nodeLabel=nodeLabel[-idx]
    nodeNumLabel=nodeNumLabel[-idx,,drop = FALSE]
    
    
    ###################################################
    nodeCutValue[currentNode] = 0;
    nodeCutValue=nodeCutValue[-idx]
    
    if(ppTree$method!="regression"){
      nodeLabel=colnames(nodeNumLabel)[max.col(nodeNumLabel)] ## "random"
      #nodeLabel[which(rowSums(structure$nodeNumLabel)==0),]=0
    }else{
      nodeLabel=as.character(nodeNumLabel[,1])
    }
    
    if(all(nodeCutValue==0)){
      prediction=rep(nodeLabel,n)
    }else{
      prediction = .Call('_ODRF_predict_ODT', PACKAGE = 'ODRF',Xnew, 
                         nodeRotaMat,nodeCutValue,childNode,nodeLabel)$prediction
    }
    
    if(ppTree$method!="regression"){
      err <- mean(prediction != ynew)
    }else{
      err <- mean((as.numeric(prediction)-ynew)^2)
    }
    pruneError[ncut,]=c(currentNode-1,length(nodeCutValue),max(structure$nodeDepth[-idx]),err)
    
    if(err<err0){
      err0=err
      
      structure$nodeRotaMat=nodeRotaMat
      structure$nodeCutValue=nodeCutValue
      structure$childNode=childNode
      #structure$nodeLabel=nodeLabel
      structure$nodeNumLabel=nodeNumLabel
      
      structure$nodeDepth=structure$nodeDepth[-idx]
      structure$nodeCutIndex[currentNode]=0;
      structure$nodeCutIndex=structure$nodeCutIndex[-idx];
    }
    
    ncut=ncut-1
    currentNode=cutNode[max(ncut,1)]
  }
  
  #structure$nodeLabel=as.character(structure$nodeLabel)
  colnames(nodeRotaMat)=c("var","node","coef")
  
  ppTree$structure=structure
  pruneError=pruneError[(ncut+1):(length(cutNode)+1),]
  ppTree$pruneError=pruneError[order(pruneError[,1],decreasing = TRUE),]
  #ppTree$tree$MaxDepth=MaxDepth
  
  class(ppTree) <- c('ODT',"prune.ODT")
  return(ppTree)
}
