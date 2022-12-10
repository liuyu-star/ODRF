#' ppRF forest Generator
#'
#' Creates an ensemble of carPPtrees using X(samplesXfeatures).
#'
#' The following parameters can be set :
#'method='c';oobe='n';weights=1;

#' @param       ntrees       : number of trees in the ensemble (default 50).
#'
#' @param       oobe         : out-of-bag error calculation, values ('y'/'n' -> yes/no) (default 'n').
#'
#' @param       nsamtosample : number of randomly selected (withreplacement) samples to use to grow each tree (default num_samples)
#' 
#'
#' Furthermore the following parameters can be set regarding the trees themselves :
# 
#' @param       method       : the criterion used for splitting the nodes
#'                           'g' : gini impurity index (classification)
#'                           'c' : information gain (classification)
#'                           'r' : squared error (regression)
#'
#' @param       minparent    : the minimum amount of samples in an impure node for it to be considered for splitting
#'
#' @param       minleaf      : the minimum amount of samples in a leaf
#'
#' @param       weights      : a vector of values which weigh the samples when considering a split
#'
#' @param       nvartosample : the number of (randomly selected) variables to consider at each node 

#FUN="RandMatPPR";paramList=NULL;cat.map=NULL;minparent=2;minleaf=1;
#nvartosample=round(sqrt(ncol(X)));ntrees=500;nsamtosample=length(y);
#method='c';numOOB=c(0,ceiling(length(y)/3))[1];weights=1;
ppRF = function(X,y,FUN="RandMatPPR",paramList=NULL,cat.map=NULL,minparent=2,minleaf=1,
                nvartosample=round(sqrt(ncol(X))),ntrees=500,nsamtosample=length(y),
                method='c',numOOB=c(0,ceiling(length(y)/3))[1],weights=1)
{
  if(tolower(method)%in%c('c','g')){
    ppForest <- list(Levels = NULL, params = NULL,trees = NULL);
    
    # adjust y to go from 1 to num.class if needed
    if (is.factor(y)) {
      ppForest$Levels <- levels(y)
      y <- as.integer(y)
    } else if (is.numeric(y)) {
      ppForest$Levels <- sort(unique(y))
      y <- as.integer(as.factor(y))
    } else {
      stop("Incompatible X type. y must be of type factor or numeric.")
    }
  }else{
    ppForest <- list(params = NULL,trees = NULL);
    y=c(y)
  }
  ppForest$params=list(n=length(y),p=ncol(X),minparent=minparent,minleaf=minleaf,nvartosample=nvartosample,ntrees=ntrees,
                       nsamtosample=nsamtosample,method=method,numOOB=numOOB,weights=weights)
  
  # keep from making copies of X
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  # address na values.
  if (any(is.na(X))) {
    if (exists("na.action")) {
      stats::na.action(X, y)
    }
    if (any(is.na(X))) {
      warning("NA values exist in data matrix")
    }
  }
  
  N=length(y)  
  
  
  ppForest$trees=vector("list",ntrees);
  ppForestT=vector("list", 6);
  ppForestT[6]=c(1);
  names(ppForestT)=c("nodesparseM","nodeCutVar", "nodeCutValue","childnode" ,"nodelabel","oobErr")
  for(i in 1:ntrees){
    TDindx=sample.int(N-numOOB);
    
    ppForestT[1:5] = carPPtree(X[TDindx,],y[TDindx],FUN,paramList,cat.map, minparent,minleaf,
                               nvartosample,method,weights,ppForest$Levels);
    
    if (numOOB>0){
      NTD = setdiff(1:N,TDindx);
      tree_output = carPPtree_Prediction_Cpp(X[NTD,],ppForestT[1:5][[1]],ppForestT[1:5][[2]],
                                             ppForestT[1:5][[3]],ppForestT[1:5][[4]],ppForestT[1:5][[5]]);

      if(tolower(method)%in%c('c','g')){
        ppForestT[6]=mean(tree_output!=y[NTD]);
      }else{
        ppForestT[6]=mean((tree_output-y[NTD])^2);
      }
      
    }    
    
    ppForest$trees[[i]] = ppForestT;
  }
  
  return(ppForest)
}