#' ppRF forest Generator
#'
#' Creates an ensemble of carPPtrees using X(samplesXfeatures).
#'
#' The following parameters can be set :
#'
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
#'
#' @return ppForest
#'
#' @useDynLib ppRF1
#' @import Rcpp
#' @import doParallel
#' @importFrom foreach foreach `%dopar%`
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @importFrom stats na.action
#' 
#' @export
#'
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
##
#method='r';FUN="RandMatPPR";paramList=NULL;catMap=NULL;
#minleaf=ifelse(method=='r',5,1);
#ntrees=500;numOOB=c(0,ceiling(length(y)/5))[2];
#weights=1;numCores=0L;seed=220924;MaxLayer=Inf
ppRF = function(X,y,method='c',FUN="RandMatPPR",paramList=NULL,catMap=NULL,
                minleaf=ifelse(method=='r',5,1),storeOOB = FALSE,
                ntrees=500,numOOB=c(0,ceiling(length(y)/5))[2],
                weights=1,numCores=0L,seed=220924,MaxLayer=Inf)
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
    ppForest <- list(Levels = NULL, params = NULL,trees = NULL);
    y=c(y)
  }
  #minparent=minparent,
  ppForest$params=list(n=length(y),p=ncol(X),FUN=FUN,paramList=paramList,catMap=catMap,MaxLayer=MaxLayer,
                       minleaf=minleaf,ntrees=ntrees,method=method,numOOB=numOOB,storeOOB = storeOOB,
                       weights=weights,numCores=numCores,seed=seed)
  
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
  
  
  N=length(y);Levels=ppForest$Levels
  PPtree=function(...){
    ppForestT=list()
    
    TDindx=sample.int(N-numOOB);
    ppForestT=c(ppForestT,carPPtree(X[TDindx,],y[TDindx],method,FUN,paramList,catMap,MaxLayer,
                                   minleaf,weights,Levels));
    
    oobErr=1
    if ((numOOB>0)&storeOOB){
      NTD = setdiff(1:N,TDindx);
      tree_output = PPtreePredict(X[NTD,],ppForestT);
      
      if(tolower(method)%in%c('c','g')){
        oobErr=mean(tree_output!=y[NTD]);
      }else{
        oobErr=mean((tree_output-y[NTD])^2);
      }
    }
    ppForestT=c(ppForestT,oobErr=oobErr)
    
    return(ppForestT)
  }
  
  
  if (numCores != 1L) {
    #RNGkind("L'Ecuyer-CMRG")
    if (numCores == 0) {
      # Use all but 1 core if numCores=0.
      numCores <- parallel::detectCores()- 1L #logical = FALSE
    }
    numCores <- min(numCores, ntrees)
    gc()
    
    #cl <- parallel::makePSOCKcluster(num.cores)
    #library("ppRF1")
    library(foreach)
    #foreach::registerDoSEQ()
    cl <- parallel::makeCluster(numCores,type = ifelse(.Platform$OS.type == "windows","PSOCK","FORK"))
    chunks <- parallel::clusterSplit(cl,seq(ntrees))
    doParallel::registerDoParallel(cl,numCores)
    set.seed(seed)
    ppForestT = foreach::foreach(co=seq(length(chunks)),.combine=list,.multicombine = TRUE,
                                 .packages = "ppRF1",.noexport ="ppForest")%dopar%{
                                   E=vector("list",length(chunks[[co]]))
                                   for(ico in chunks[[co]]) {
                                     #set.seed(ico)
                                     #ppForestT = try(PPtree(), silent=TRUE)
                                     #if ('try-error' %in% class(ppForestT)){
                                     #  next
                                     #}
                                     E[[ico-(chunks[[co]][1]-1)]] = PPtree()#ppForestT
                                   }
                                   # return local results
                                   E
                                 }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
    for (i in 1:numCores) {
      ppForest$trees=c(ppForest$trees,ppForestT[[i]]) 
    }
  } else {
    # Use just one core.
    set.seed(seed)
    ppForest$trees <- lapply(1:ntrees, PPtree)
  }
  
  return(ppForest)
}