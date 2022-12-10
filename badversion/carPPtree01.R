#' RerF Tree Generator
#' 
#' Grows a CARTree using X(samplesXfeatures) and one of the following criteria which can be set via the parameter 'method'
# 
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
#X=X;y=y;FUN="RandMatPPR";paramList=NULL;catMap=NULL; minparent=2;minleaf=1;nvartosample=ncol(X);method='c';weights=1
carPPtree=function(X,y,FUN="RandMatPPR",paramList=NULL,catMap=NULL, minparent=2,minleaf=1,
                   method='c',weights=1,Levels=levels(as.factor(y)))
{
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  y <- as.integer(as.factor(y));
  Levels=as.character(Levels);
  
  N = length(y);
  M = ncol(X);
  
  paramList <- defaults(ncolX = ncol(X), paramList = paramList, cat.map = catMap)
  #FUN <- match.fun(FUN, descend = TRUE)
  funName <- FUN
  FUN <- match.fun(FUN, descend = TRUE)
  
  # L = 2*ceiling(N);
  L = 2*ceiling(N/minleaf) - 1;
  # L = 2*ceiling(N/minleaf) + 1;
  # L = factorial(ceiling(N/minleaf));
  # L = 2^ceiling(N/minleaf);
  
  nodeXIndx = vector("list", L);
  nodeXIndx[[1]] = 1 : N;
  
  # d = ceilinging(min(sqrt(N),M/3));
  # zeros(2*M,3,L)
  nodesparseM = matrix(0,2*M*L,3);
  
  nodeCutVar = rep(0,L);
  nodeCutValue = rep(0,L);
  nodeflags = rep(0,L+1);
  nodelabel = rep(0,L);
  childnode = rep(0,L);
  nodeflags[1] = 1;
  
  if(tolower(method)%in%c('c','g')){
    max_label = length(Levels);  
  }else{
    max_label= 0;
  }
  
  #start create pptree
  ##############################################################################
  current_node = 1;
  numPP=0;
  
  while((nodeflags[current_node] == 1)&&(length(which(nodeflags == 0))!=0)){
    free_node = min(which(nodeflags == 0));
    currentXIndx = nodeXIndx[[current_node]];
    
    if(length(unique(y[currentXIndx]))==1){
      if(tolower(method)%in%c('c','g')){
        nodelabel[current_node]=Levels[y[currentXIndx[1]]];
      }else{
        nodelabel[current_node] = y[currentXIndx[1]];
      }
      #nodeCutVar[current_node] = 0;
      #nodeCutValue[current_node] = 0;
    }else{
      if(length(currentXIndx)>=minparent){
        if(length(weights)>1){
          Wcd = weights[currentXIndx];
        }else{
          Wcd = 1;
        }
        
        #  node_var = randperm(M);
        #  node_var = node_var(1:m);
        ##  node_var = (randsample(M,m,0));
        
        if(funName=="RandMatPPR"){
          sparseM <- RandMatPPR(x=X[currentXIndx,],y=y[currentXIndx],catMap=paramList$catMap)
        }else{
          sparseM <- do.call(FUN, paramList)
        }
        
        numDr=unique(sparseM[,2]);
        rotaMat=matrix(0,length(numDr),M);
        for(ndr in 1:length(numDr)){
          lrows = which(sparseM[,2] == ndr);
          rotaMat[ndr,sparseM[lrows, 1]] = sparseM[lrows, 3];
        }
        
        rotaX=as.matrix(X[currentXIndx,] %*% t(rotaMat));
        bestCut = best_cut_node(method,rotaX,y[currentXIndx],Wcd,minleaf,max_label);
        bestCutVar = bestCut$BestCutVar;
        bestCutValue = bestCut$BestCutVal;
        
        if((bestCutVar!=-1)&(bestCutVar<=ncol(rotaX))) {##bestCutVar=1614907703L?
          
          nodesparseM[numPP+(1:nrow(sparseM)),]=sparseM;
          numPP=numPP+nrow(sparseM)+1;
          
          nodeCutVar[current_node] = bestCutVar;            
          nodeCutValue[current_node] = bestCutValue;

          Lindex=which(rotaX[,bestCutVar]<bestCutValue);
          nodeXIndx[[free_node]] = currentXIndx[Lindex];
          nodeXIndx[[free_node+1]] = currentXIndx[setdiff(1:length(currentXIndx),Lindex)];
          
          nodeflags[free_node:(free_node + 1)] = 1;
          childnode[current_node]=free_node;
        }else{
          if(tolower(method)%in%c('c','g')){
            leaf_label = which.max(table(factor(y[currentXIndx],levels=Levels)));
            nodelabel[current_node]=Levels[leaf_label];
          }else{
            nodelabel[current_node] = mean(y[currentXIndx]);
          }
        }
        
      }else{
        if(tolower(method)%in%c('c','g')){
          leaf_label = which.max(table(factor(y[currentXIndx],levels=Levels)));
          nodelabel[current_node]=Levels[leaf_label];
        }else{
          nodelabel[current_node] = mean(y[currentXIndx]);
        }
      }
    }
      current_node = current_node+1;
  }
  
  return(list(nodesparseM = nodesparseM[1:numPP,],
              nodeCutVar = nodeCutVar[1:(current_node-1)],
              nodeCutValue =nodeCutValue[1:(current_node-1)],
              childnode = childnode[1:(current_node-1)],
              nodelabel = nodelabel[1:(current_node-1)])) 
}
