#' Pruning of class \code{ODT}.
#' 
#' Prune \code{ODT} from bottom to top with validation data based on prediction error.
#'
#' @param ppTree an object of class \code{ODT}.
#' @param X An n by d numeric matrix (preferable) or data frame is used to prune the object of class \code{ODT}.
#' @param y A response vector of length n.
#' @param MaxDepth The maximum depth of the tree after pruning. (Default 1)
#' @param ... optional parameters to be passed to the low level function.
#' 
#' @return \itemize{an object of class \code{ODT} and \code{prune.ODT}. 
#' \item{\code{ppTree} The same result as \code{ODT}.}
#' \item{\code{pruneError} Error of validation data after each pruning, misclassification rate (MR) for classification or mean square error (MSE) for regression.}
#'}
#' @seealso \code{\link{ODT}} \code{\link{plot.prune.ODT}} \code{\link{prune.ODRF}} \code{\link{online.ODT}}
#' 
#' @examples
#' #Classification with Oblique Decision Tree
#' data(seeds)
#' set.seed(221212)
#' train = sample(1:209,100)
#' train_data = data.frame(seeds[train,])
#' test_data = data.frame(seeds[-train,])
#' index=seq(floor(nrow(train_data)/2))
#' 
#' tree = ODT(varieties_of_wheat~.,train_data[index,],type='i-classification')
#' tree = prune(tree,train_data[-index,-8],train_data[-index,8])
#' pred <- predict(tree,test_data[,-8])
#' #classification error
#' (mean(pred!=test_data[,8]))
#' 
#' 
#' #Regression with Oblique Decision Tree
#' data(body_fat)
#' set.seed(221212)
#' train = sample(1:252,100)
#' train_data = data.frame(body_fat[train,])
#' test_data = data.frame(body_fat[-train,])
#' index=seq(floor(nrow(train_data)/2))
#' 
#' tree = ODT(Density~.,train_data[index,],type='regression')
#' tree = prune(tree,train_data[-index,-1],train_data[-index,1])
#' pred <- predict(tree,test_data[,-1])
#' #estimation error
#' mean((pred-test_data[,1])^2)
#'
#' @rdname prune.ODT
#' @aliases prune.ODT
#' @method prune ODT
#' @export
prune.ODT=function(ppTree,X,y,MaxDepth=1,...)
{
  structure=ppTree$structure
  if(!is.null(MaxDepth)){
    MaxDepth=min(MaxDepth,max(structure$nodeDepth))
  }
  numNode=length(structure$nodeCutValue)
  
  #vars=all.vars(ppTree$terms)
  Xcat=ppTree$data$Xcat
  catLabel=ppTree$data$catLabel
  if(sum(Xcat)>0&is.null(catLabel)){
    #vars=vars[-(1+seq(length(unlist(catLabel))))]
    #}else{
    stop("'Xcat!=0' however 'catLabel' does not exist!") 
  }
  # address na values.
  data=data.frame(y,X)
  if (any(is.na(data))) {
    data=ppTree$data$na.action(data.frame(data))
    warning("NA values exist in data matrix")
  }
  
  #ynew= data[,setdiff(colnames(data),vars[-1])]
  #Xnew= data[,vars[-1]]
  ynew= data[,1]
  Xnew= data[,-1]
  Xnew=as.matrix(Xnew) 
  rm(data);rm(X);rm(y)
  
  #weights0=c(ppTree$data$weights,ppTree$paramList$weights)
  #if(!is.null(ppTree$data$weights))
  #  Xnew <- Xnew * matrix(weights,length(y),ncol(Xnew))
  
  p=ncol(Xnew)
  n=nrow(Xnew)
  
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
    
    Xnew=cbind(Xnew1,apply(Xnew[,-Xcat], 2, as.numeric)) 
    p=ncol(Xnew)
    numCat=length(unlist(catLabel))
    rm(Xnew1)
    rm(Xnewj)
  }
  Xnew=as.matrix(Xnew)
  colnames(Xnew)=ppTree$data$varName
  
  if(!is.null(ppTree$data$subset))
    Xnew=Xnew[ppTree$data$subset,]
  
  #Variable scaling.
  if(ppTree$data$Xscale!="No"){
    indp=(numCat+1):p
    Xnew[,indp]=(Xnew[,indp,drop = FALSE]-matrix(ppTree$data$minCol,n,length(indp),byrow = T))/
      matrix(ppTree$data$maxminCol,n,length(indp),byrow = T)
  }
  
  if (ppTree$data$TreeRandRotate) {
    Xnew[,ppTree$data$rotdims] <- Xnew[,ppTree$data$rotdims,drop = FALSE] %*% ppTree$data$rotmat
  }
  
  
  if(ppTree$type!="regression"){
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
  
  if(ppTree$type!="regression"){
    err0 <- mean(prediction != ynew)
  }else{
    #e.0 = mean((ynew-mean(y))^2)
    err0 <- mean((as.numeric(prediction)-ynew)^2)#/e.0
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
    #nodeLabel[currentNode]=ifelse(ppTree$type!="regression",nodeLabel[idx][which.max(structure$nodeNumLabel[idx])],
    #                              structure$nodeNumLabel[idx]*nodeLabel[idx]/sum(structure$nodeNumLabel[idx]))
    if(ppTree$type!="regression"){
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
    
    if(ppTree$type!="regression"){
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
    
    if(ppTree$type!="regression"){
      err <- mean(prediction != ynew)
    }else{
      err <- mean((as.numeric(prediction)-ynew)^2)#/e.0
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
  #rownames(nodeRotaMat)=rep(structure$nodeDepth,table(nodeRotaMat[,2]))
  #rownames(nodeNumLabel)=structure$nodeDepth
  
  ppTree$structure=structure
  pruneError=pruneError[(ncut+1):(length(cutNode)+1),]
  ppTree$pruneError=pruneError[order(pruneError[,1],decreasing = TRUE),]
  #ppTree$tree$MaxDepth=MaxDepth
  
  class(ppTree) <- c('ODT',"prune.ODT")
  return(ppTree)
}
