#' A machine learning algorithm for training class \code{ODT}.
#' 
#' The \code{\link{ODT}}is constantly updated by multiple batches of data to optimize the model.
#'
#' @param obj an object of class \code{ODT}.
#' @param data Training data of class \code{data.frame} is used to update the object of class \code{ODT}.
#' @param weights A vector of length same as \code{data} that are positive weights. (default NULL)
#' 
#' @return The same result as \code{ODT}.
#'
#' @seealso \code{ODT} \code{prune.ODT}
#' 
#' @examples
#' library(ODRF)
#' 
#' #Classification with Oblique Decision Tree
#' data(seeds)
#' set.seed(221212)
#' train = sample(1:209,100)
#' train_data = data.frame(seeds[train,])
#' test_data = data.frame(seeds[-train,])
#' 
#' tree = ODT(varieties_of_wheat~.,train_data[seq(floor(nrow(train_data)/2)),],type='i-classification')
#' tree = online(tree,train_data[-seq(floor(nrow(train_data)/2)),])
#' pred <- predict(tree,test_data[,-8])
#' #estimation error
#' (mean(pred!=test_data[,8]))
#' 
#' #Regression with Oblique Decision Tree
#' data(body_fat)
#' set.seed(221212)
#' train = sample(1:252,100)
#' train_data = data.frame(body_fat[train,])
#' test_data = data.frame(body_fat[-train,])
#' 
#' tree = ODT(Density~.,train_data[seq(floor(nrow(train_data)/2)),],type='regression')
#' tree = online(tree,train_data[-seq(floor(nrow(train_data)/2)),])
#' pred <- predict(tree,test_data[,-8])
#' #estimation error
#' mean((pred-test_data[,1])^2)
#' 
#' @import Rcpp
#' @aliases online.ODRF
#' @rdname online.ODRF
#' @method online ODRF
#' @export
online.ODT=function(ppTree,data,weights=NULL)
{
  weights0=weights
  Call=ppTree$call
  Terms=ppTree$terms
  type=ppTree$type
  Levels=ppTree$Levels
  NodeRotateFun=ppTree$NodeRotateFun
  paramList=ppTree$paramList
  ppTree=ppTree[-seq(6)]
  ppTreeVar=c(names(ppTree$data),names(ppTree$tree),names(ppTree$structure))
  ppTree=do.call("c",ppTree)
  
  for(v in seq(length(ppTreeVar))){
    assign(ppTreeVar[v], ppTree[[v]])
  }
  rm(ppTree)
  
  vars=all.vars(Terms)
  # address na values.
  if (any(is.na(data))) {
    data=na.action(data.frame(data))
    warning("NA values exist in data matrix")
  }
  y= data[,setdiff(colnames(data),vars[-1])]
  X= data[,vars[-1]]
  X=as.matrix(X) 
  
  if(!is.null(subset))
    X=X[subset,]
  #weights=c(weights,paramList$weights)
  if(!is.null(weights))
    X <- X * matrix(weights0,length(y),ncol(X))
  weights=weights0
  
  n = length(y);
  p = ncol(X);

  
  if(!NodeRotateFun%in%ls("package:ODRF")){
    source(paste0(FunDir,"/",NodeRotateFun,".R"))
  }
  
  #get()
  FUN <- match.fun(NodeRotateFun, descend = TRUE)
  method0=strsplit(type,split = "")[[1]][1]
  
  
  if(type!="regression"){
    if (!is.integer(y)) {
      y <- as.integer(as.factor(y));
    }
    maxLabel = length(Levels); 
    
    nodeLabel=Levels[max.col(nodeNumLabel)] ## "random"
    nodeLabel[which(rowSums(nodeNumLabel)==0)]=0
    sl=seq(maxLabel)
  }else{
    y=c(y)
    maxLabel= 0;
    
    nodeLabel=as.character(nodeNumLabel[,1])
  }
  
  
  if(all(nodeCutValue==0)){
    Nodes=rep(1,n)
  }else{
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
    if(Xscale!="No"){
      indp=(numCat+1):p
      X[,indp]=(X[,indp]-matrix(minCol,n,length(indp),byrow = T))/
        matrix(maxminCol,n,length(indp),byrow = T)
    }
    
    if (TreeRandRotate) {
      X[, rotdims] <- X[,rotdims,drop = FALSE] %*% rotmat
    }
    
    #Nodes = ppCARTNode(X, nodeRotaMat, nodeCutValue, childNode)
    #Nodes = .Call(`_ppRF_ppCARTPredict`, X, nodeRotaMat, nodeCutValue, childNode, nodeLabel)$node
    Nodes = .Call('_ODRF_predict_ODT', PACKAGE = 'ODRF',X, nodeRotaMat, nodeCutValue, childNode, nodeLabel)$node
    #, as.character(ppTree$nodeLabel)
  }
  
  paramList = ODRF:::defaults(paramList,p,catLabel,NodeRotateFun,type,weights)
  
  if(is.infinite(MaxDepth)) {
    numNode=max(numNode,sum(2^(0:ceiling(log2(n/MinLeaf))))) 
  }else{
    MaxDepth=min(MaxDepth,ceiling(log2(n/MinLeaf)),n/MinLeaf-1)
    numNode=min(numNode,sum(2^(0:MaxDepth)))
  }
  
  numNode0=length(nodeCutValue)
  nodeX=sort(unique(Nodes))
  nodeXIndx = vector("list", numNode-numNode0+1);
  nodeXIndx[seq(numNode0)]=NA
  for (nx in nodeX) {
    nodeXIndx[[nx]]=which(Nodes==nx)
  }
  cutNode=which(is.na(nodeXIndx[seq(numNode0)]))
  
  rep0=rep(0,numNode-numNode0)
  childNode0 = childNode;
  #nodeNumLabel0=c(nodeNumLabel,rep0);
  nodeNumLabel0=nodeNumLabel
  #names(nodeNumLabel0)=paste0("v",seq(numNode-numNode0))
  #names(nodeNumLabel0)[seq(numNode0)]=nodeLabel

  nodeDepth = c(nodeDepth,rep0);
  nodeCutValue =c(nodeCutValue,rep0)
  nodeCutIndex =c(nodeCutIndex,rep0)
  childNode = c(childNode,rep0);
  #nodeLabel = c(nodeLabel,rep0);
  #nodeNumLabel=c(nodeNumLabel,rep0);
  
  
  #start create pptree nodeFlags
  ##############################################################################
  currentNode = nodeX[1];
  freeNode=ifelse(currentNode==1,2,childNode[max(which(childNode[seq(currentNode)]!=0))]+2)
  while(!is.null(nodeXIndx[[currentNode]])){
    if(is.na(nodeXIndx[currentNode])){
      if(nodeCutValue[currentNode]!=0){
        freeNode=freeNode+2 
      }
      currentNode = currentNode+1;
      next
    }
    
    if((length(unique(y[nodeXIndx[[currentNode]]]))==1)||
       (length(nodeXIndx[[currentNode]])<=(MinLeaf+1))||
       (nodeDepth[currentNode]>=MaxDepth)||
       (freeNode>=numNode)){
      nn=length(y[nodeXIndx[[currentNode]]])
      #r=ceiling(nn*nodeNumLabel0[currentNode])
      #parentLabel=nn*nodeNumLabel0[currentNode,]
      if(type!="regression"){
        leafLabel = table(Levels[c(sl,y[nodeXIndx[[currentNode]]])])-1+nn*nodeNumLabel0[currentNode,]
        #leafLabel = table(c(Levels[y[nodeXIndx[[currentNode]]]],rep(names(nodeNumLabel0[currentNode]),r)))
        #nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        #nodeNumLabel[currentNode]=max(leafLabel)
      }else{
        #nodeLabel[currentNode]=(r*as.numeric(names(nodeNumLabel0[currentNode]))+nn*mean(y[nodeXIndx[[currentNode]]]))/(r+nn);
        #nodeNumLabel[currentNode]=nn
        leafLabel=c((nn*nodeNumLabel0[currentNode,2]*nodeNumLabel0[currentNode,1]+
                     nn*mean(y[nodeXIndx[[currentNode]]]))/(nn*nodeNumLabel0[currentNode,2]+nn))
        leafLabel=c(leafLabel,(nn*nodeNumLabel0[currentNode,2]+nn))
      }
      #nodeNumLabel0=rbind(nodeNumLabel0,leafLabel)
      
      if(currentNode>numNode0){
        nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
        nodeNumLabel=rbind(nodeNumLabel,leafLabel)
      }else{
        nodeNumLabel[currentNode,]=leafLabel
      }
      
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
                        model = paramList$model,type=paramList$type)#"NNet"
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
      nn=length(y[nodeXIndx[[currentNode]]])
      #r=ceiling(nn*nodeNumLabel0[currentNode])
      #parentLabel=nn*nodeNumLabel0[currentNode,]
      if(type!="regression"){
        leafLabel = table(Levels[c(sl,y[nodeXIndx[[currentNode]]])])-1+nn*nodeNumLabel0[currentNode,]
        #leafLabel = table(c(Levels[y[nodeXIndx[[currentNode]]]],rep(names(nodeNumLabel0[currentNode]),r)))
        #nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        #nodeNumLabel[currentNode]=max(leafLabel)
      }else{
        #nodeLabel[currentNode]=(r*as.numeric(names(nodeNumLabel0[currentNode]))+nn*mean(y[nodeXIndx[[currentNode]]]))/(r+nn);
        #nodeNumLabel[currentNode]=nn
        leafLabel=c((nn*nodeNumLabel0[currentNode,2]*nodeNumLabel0[currentNode,1]+
                       nn*mean(y[nodeXIndx[[currentNode]]]))/(nn*nodeNumLabel0[currentNode,2]+nn))
        leafLabel=c(leafLabel,(nn*nodeNumLabel0[currentNode,2]+nn))
      }
      #nodeNumLabel0=rbind(nodeNumLabel0,leafLabel)
      
      if(currentNode>numNode0){
        nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
        nodeNumLabel=rbind(nodeNumLabel,leafLabel)
      }else{
        nodeNumLabel[currentNode,]=leafLabel
      }
      
      nodeXIndx[currentNode]=NA
      currentNode = currentNode+1;
      next
    }
    
    sparseM=sparseM[sparseM[,2]==numDr[bestCut$BestCutVar],,drop = FALSE]
    sparseM[,2]=currentNode
    if(currentNode<=numNode0){
      #r=rep(nodeNumLabel0[currentNode]/length(y[nodeXIndx[[currentNode]]]),2)
      #names(r)=rep(nodeLabel[currentNode],2)
      LRnode= rbind(nodeNumLabel0[currentNode,],nodeNumLabel0[currentNode,])
      if(!is.na(childNode0[currentNode])) {
        if(type!="regression"){
          #r=rep(nodeNumLabel0[currentNode],2)
          #names(r)=rep(names(nodeNumLabel0[currentNode]),2)
          LRnode= LRnode/length(y[nodeXIndx[[currentNode]]])
        }else{
          LRnode[,2]= LRnode[,2]/length(y[nodeXIndx[[currentNode]]])
        }
      }
      
      idx=which(nodeRotaMat[,2]==currentNode)
      nodeRotaMat=rbind(nodeRotaMat[seq(idx-1),,drop = FALSE],sparseM,nodeRotaMat[-seq(idx),,drop = FALSE])
      
      idx=which(childNode0!=0)
      idx=idx[idx>currentNode]
      childNode0[idx]=childNode0[idx]+2
      childNode0[currentNode]=freeNode;

      if(freeNode<=numNode0){
        numNode0=numNode0+2
        idx=seq(freeNode-1)
        nodeXIndx=c(nodeXIndx[idx],list(NA,NA),nodeXIndx[-idx])
        childNode0=c(childNode0[idx],NA,NA,childNode0[-idx])
        
        nodeDepth=c(nodeDepth[idx],0,0,nodeDepth[-idx])
        nodeCutValue=c(nodeCutValue[idx],0,0,nodeCutValue[-idx])
        nodeCutIndex=c(nodeCutIndex[idx],0,0,nodeCutIndex[-idx])
        
        #nodeLabel=c(nodeLabel[idx],0,0,nodeLabel[-idx])
        #nodeNumLabel=c(nodeNumLabel[idx],0,0,nodeNumLabel[-idx])
        #nodeNumLabel0=c(nodeNumLabel0[idx],r,nodeNumLabel0[-idx])
        nodeNumLabel=rbind(nodeNumLabel[idx,],0,0,nodeNumLabel[-idx,,drop = FALSE])
        nodeNumLabel0=rbind(nodeNumLabel0[idx,],LRnode,nodeNumLabel0[-idx,,drop = FALSE])

        idx=1:max(which(nodeRotaMat[,2]==(freeNode-1)))
        nodeRotaMat=rbind(nodeRotaMat[idx,,drop = FALSE],0,0,nodeRotaMat[-idx,,drop = FALSE])
        nodeRotaMat[-idx,2]=c(freeNode,freeNode+1,nodeRotaMat[-c(1,2,idx+2),2]+2)
      }else{
        #nodeNumLabel0[freeNode+c(0,1)]=r
        #names(nodeNumLabel0)[freeNode+c(0,1)]=names(r)
        nodeNumLabel0=rbind(nodeNumLabel0,LRnode)
        #nodeNumLabel=rbind(nodeNumLabel,0,0)
      }
    }else{
      #nodeNumLabel0[freeNode+c(0,1)]=rep(nodeNumLabel0[currentNode],2)
      #names(nodeNumLabel0)[freeNode+c(0,1)]=rep(names(nodeNumLabel0[currentNode]),2)
      nodeNumLabel=rbind(nodeNumLabel,0)
      nodeNumLabel0=rbind(nodeNumLabel0,nodeNumLabel0[currentNode,],nodeNumLabel0[currentNode,])
      nodeRotaMat=rbind(nodeRotaMat,sparseM)
      #childNode[currentNode]=freeNode;
    }
    rm(rotaX)
    rm(sparseM)
    
    nodeCutValue[currentNode] = bestCut$BestCutVal;
    nodeCutIndex[currentNode]=bestCut$BestIndex[bestCut$BestCutVar]
    childNode[currentNode]=freeNode;
    
    #Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
    nodeXIndx[[freeNode]] = nodeXIndx[[currentNode]][Lindex];
    nodeXIndx[[freeNode+1]] = nodeXIndx[[currentNode]][setdiff(1:length(nodeXIndx[[currentNode]]),Lindex)];
    nodeDepth[freeNode+c(0,1)]=nodeDepth[currentNode]+1
    
    nodeXIndx[currentNode]=NA
    currentNode = currentNode+1;
    freeNode =freeNode+2
  }
  childNode0[is.na(childNode0)]=0
  childNode[1:numNode0]=childNode0
  nodeDepth=nodeDepth[1:(currentNode-1)]
  colnames(nodeRotaMat)=c("var","node","coef")
  rownames(nodeRotaMat)=rep(nodeDepth,table(nodeRotaMat[,2]))
  rownames(nodeNumLabel)=nodeDepth
  
  ppTree=list(call=Call,terms=Terms,type=type,Levels=Levels,NodeRotateFun=NodeRotateFun,paramList=paramList)
  ppTree$data=list(subset=subset,weights=weights,na.action=na.action,n=n,p=p,varName=varName,
                   Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel,
                   TreeRandRotate=TreeRandRotate,rotdims=rotdims,rotmat=rotmat);
  ppTree$tree=list(FunDir=FunDir,MaxDepth=MaxDepth,MinLeaf=MinLeaf,numNode=numNode);
  ppTree$structure=list(nodeRotaMat = nodeRotaMat,nodeNumLabel = nodeNumLabel,nodeCutValue =nodeCutValue[1:(currentNode-1)],
                        nodeCutIndex =nodeCutIndex[1:(currentNode-1)],childNode = childNode[1:(currentNode-1)],
                        nodeDepth=nodeDepth)
  class(ppTree) <- "ODT"
  return(ppTree)
}
