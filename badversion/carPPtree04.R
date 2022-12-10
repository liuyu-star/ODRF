#' RerF Tree Generator
#' 
#' Grows a CARTree using X(samplesXfeatures) and one of the following criteria which can be set via the parameter 'method'
#'
#' @param       'g' : gini impurity index (classification)
#' @param       'c' : information gain (classification, default)
#' @param       'r' : squared error (regression)
#' 
#' Other parameters that can be set are:
#' 
#' @param paramList parameters in a named list to be used by FUN. If left unchanged,
#' @param catMap a list specifying which columns in X correspond to the same one-of-K encoded feature. Each element of catMap is a numeric vector specifying the K column indices of X corresponding to the same categorical feature after one-of-K encoding. All one-of-K encoded features in X must come after the numeric features. The K encoded columns corresponding to the same categorical feature must be placed contiguously within X. The reason for specifying catMap is to adjust for the fact that one-of-K encoding cateogorical features results in a dilution of numeric features, since a single categorical feature is expanded to K binary features. If catMap = NULL, then RerF assumes all features are numeric (i.e. none of the features have been one-of-K encoded).
#' @param FUN :FUN=RandMatPPR.
#' @param minparent    : the minimum amount of samples in an impure node for it to be considered for splitting (default 2)
#' @param minleaf      : the minimum amount of samples in a leaf (default 1)
#' @param weights      : a vector of values which weigh the samples when considering a split (default [])
#' @param nvartosample : the number of (randomly selected) variables to consider at each node (default all)
#' @return ppTree
#'
#' @useDynLib ppRF1
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
#X=X;y=y;method='r';FUN="RandMatPPR";paramList=NULL;catMap=NULL;MaxLayer=Inf;
#minleaf=ifelse(method=='r',5,1);weights=1;
#Levels=if(method=='r')NULL else levels(as.factor(y))
carPPtree=function(X,y,method='c',FUN="RandMatPPR",paramList=NULL,catMap=NULL,MaxLayer=Inf,numNode=Inf,
                   minleaf=ifelse(method=='r',5,1),weights=1,Levels=if(method=='r')NULL else levels(as.factor(y)))
{
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if(tolower(method)%in%c('c','g')){
    if (!is.integer(y)) {
      y <- as.integer(as.factor(y));
    }
    maxLabel = length(Levels); 
  }else{
    y=c(y)
    maxLabel= 0;
  }
  
  
  N = length(y);
  M = ncol(X);
  
  if(is.null(paramList)){
    paramList <- ppRF1:::defaults(ncolX = ncol(X), paramList = paramList, cat.map = catMap)
  }
  #FUN <- match.fun(FUN, descend = TRUE)
  funName <- FUN
  FUN <- match.fun(FUN, descend = TRUE)
  
  #L = 2*ceiling(N);
  #L = 2*ceiling(N/minleaf) - 1;
  # L = 2*ceiling(N/minleaf) + 1;
  # L = factorial(ceiling(N/minleaf));
  # L = 2^ceiling(N/minleaf);
  #MaxLayer=min(ceiling(log2(N/minleaf)),MaxLayer)
  #MaxLayer=min(MaxLayer,N/minleaf-1) 
  
  if(is.infinite(MaxLayer)) {
    L=sum(2^(0:ceiling(log2(N/minleaf))))
    if(!is.infinite(numNode)){
      L=numNode
    }
    index=c()
  }else{
    MaxLayer=min(MaxLayer,ceiling(log2(N/minleaf)),N/minleaf-1)
    L=sum(2^(0:MaxLayer))
    if(!is.infinite(numNode)){
      L=numNode
    }
    index=2:L;
    totalN=0;
    totalnode = rep(0,L);
    totalnode[1] = 1;
  }

  #MaxLayer=N/minleaf-1
  #L=2*MaxLayer+1
  
  nodeXIndx = vector("list", L);
  nodeXIndx[[1]] = 1 : N;
  
  # d = ceilinging(min(sqrt(N),M/3));
  # zeros(2*M,3,L)
  #s = c(0,round(quantile(1:M, (1:d)/d)))
  #max(s[-1] - s[-d1 - 1])
  d = ceiling(min(sqrt(N),M/3))
  nodesparseM = matrix(0,(M+2*d)*L,3);
  
  nodeCutVar = rep(0,L);
  nodeCutValue = nodeCutVar;
  nodelabel = nodeCutVar;
  childnode = nodeCutVar;
  nodeflags = rep(0,L+1);
  nodeflags[1] = 1;
  
  
  #start create pptree
  ##############################################################################
  current_node = 1;
  numPP=0;
  free_node = 2
  while(nodeflags[current_node] == 1){
    #free_node = min(which(nodeflags == 0));

    if(!is.infinite(MaxLayer)){
      #index=which(totalnode==0)
      #if(totalnode[current_node]==-1){
      repeat{
        if((totalnode[current_node]!=-1)||(length(index)==0)){break}
        totalnode[index[1:2]] = -1
        totalnode=totalnode[-current_node]
        index=index[-c(1,2,L-totalN)]-1
        totalN=totalN+1
      }
      #}
    }
    
    if((length(unique(y[nodeXIndx[[current_node]]]))==1)||
       (length(nodeXIndx[[current_node]])<=(minleaf+1))||
       (!is.infinite(MaxLayer)&(length(index)==0))||
       (sum(nodeflags)>=L)){
      if(tolower(method)%in%c('c','g')){
        #leaf_label = which.max(table(factor(y[nodeXIndx[[current_node]]],levels=Levels)));
        leaf_label = table(Levels[y[nodeXIndx[[current_node]]]])
        nodelabel[current_node]=names(leaf_label)[which.max(leaf_label)];
      }else{
        nodelabel[current_node] = mean(y[nodeXIndx[[current_node]]]);
      }
      if(length(index)>0){
        totalnode[index[c(1,2)]] = -1;
        index=index[-(1:2)]
      }
      
      current_node = current_node+1;
      next
    }
    
    
    if(length(weights)>1){
      Wcd = weights[nodeXIndx[[current_node]]];
    }else{
      Wcd = 1;
    }
    
    #  node_var = randperm(M);
    #  node_var = node_var(1:m);
    ##  node_var = (randsample(M,m,0));
    
    if(funName=="RandMatPPR"){
      sparseM <- RandMatPPR(x=X[nodeXIndx[[current_node]],],y=y[nodeXIndx[[current_node]]],method=method,catMap=paramList$catMap)
    }else{
      sparseM <- do.call(FUN, paramList)
    }
    
    numDr=unique(sparseM[,2]);
    rotaX=matrix(0,length(numDr),M);
    for(ndr in 1:length(numDr)){
      lrows = which(sparseM[,2] == ndr);
      rotaX[ndr,sparseM[lrows, 1]] = sparseM[lrows, 3];
    }
    
    rotaX=as.matrix(X[nodeXIndx[[current_node]],] %*% t(rotaX));
    bestCut =ppRF1:::best_cut_node(method,rotaX,y[nodeXIndx[[current_node]]],Wcd,minleaf,maxLabel);
    bestCutVar = bestCut$BestCutVar;
    bestCutValue = bestCut$BestCutVal;
    
    
    if(bestCutVar==-1){
      TF=TRUE
    }else{
      Lindex=which(rotaX[,bestCutVar]<bestCutValue);
      TF=min(length(Lindex),length(nodeXIndx[[current_node]])-length(Lindex))<=minleaf;
    }
    if(TF){##bestCutVar=1614907703L?(bestCutVar>ncol(rotaX))
      if(tolower(method)%in%c('c','g')){
        #leaf_label = which.max(table(factor(y[nodeXIndx[[current_node]]],levels=Levels)));
        leaf_label = table(Levels[y[nodeXIndx[[current_node]]]])
        nodelabel[current_node]=names(leaf_label)[which.max(leaf_label)];
      }else{
        nodelabel[current_node] = mean(y[nodeXIndx[[current_node]]]);
      }
      if(!is.infinite(MaxLayer)){
        totalnode[index[c(1,2)]] = -1;
        index=index[-(1:2)]
      }
      
      current_node = current_node+1;
      next
    }
    
    nodesparseM[numPP+(1:nrow(sparseM)),]=sparseM;
    numPP=numPP+nrow(sparseM)+1;
    rm(rotaX)
    rm(sparseM)
    
    nodeCutVar[current_node] = bestCutVar;            
    nodeCutValue[current_node] = bestCutValue;
    
    #Lindex=which(rotaX[,bestCutVar]<bestCutValue);
    nodeXIndx[[free_node]] = nodeXIndx[[current_node]][Lindex];
    nodeXIndx[[free_node+1]] = nodeXIndx[[current_node]][setdiff(1:length(nodeXIndx[[current_node]]),Lindex)];
    
    nodeflags[free_node:(free_node + 1)] = 1;
    childnode[current_node]=free_node;
    if(!is.infinite(MaxLayer)){
      totalnode[index[c(1,2)]] = free_node:(free_node + 1);
      index=index[-(1:2)]
    }
    free_node = free_node+2;
    
    current_node = current_node+1;
  }
  
  #childnode1 = childnode[1:(current_node-1)]
  #childnode1[which(childnode1[1:which.max(childnode1)]==0)]=2
  #totalnode= sum(childnode1[1:which.max(childnode1)])+1 
  #MaxLayer=min(which(totalnode<cumsum(L)))-1
  #MaxLayer=ceiling(log2(length(which(totalnode==1))+max(totalnode)))
  #totalnode=length(which(totalnode==1))+which.max(totalnode)
  #MaxLayer=min(which(totalnode<cumsum(L)))-1
  return(list(nodesparseM = nodesparseM[1:numPP,],
              nodeCutVar = nodeCutVar[1:(current_node-1)],
              nodeCutValue =nodeCutValue[1:(current_node-1)],
              childnode = childnode[1:(current_node-1)],
              nodelabel = nodelabel[1:(current_node-1)],
              method=method)) 
}
