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
#method='c';FUN="RandMatPPR";paramList=NULL;catMap=NULL;
#minleaf=ifelse(method=='r',5,1);rotate = FALSE;replacement = TRUE;stratify = TRUE;
#ntrees=500;numOOB=c(0,ceiling(length(y)/5))[2];storeOOB = FALSE;
#weights=1;numCores=0L;seed=220924;MaxLayer=Inf;numNode=Inf
ppRF = function(X,y,method='c',FUN="RandMatPPR",paramList=NULL,catMap=NULL,
                minleaf=ifelse(method=='r',5,1),storeOOB = FALSE,
                rotate = FALSE,replacement = TRUE, stratify = TRUE,
                ntrees=500,numOOB=c(0,ceiling(.368*length(y)))[2],
                weights=1,numCores=0L,seed=220924,MaxLayer=Inf,numNode=Inf)
{
  if(tolower(method)%in%c('c','g')){
    ppForest <- list(Levels = NULL, params = NULL,trees = NULL);
    
    # adjust y to go from 1 to numClass if needed
    if (is.factor(y)) {
      ppForest$Levels <-levels(y)
      y <- as.integer(y)
    } else if (is.numeric(y)) {
      ppForest$Levels <- sort(unique(y))
      y <- as.integer(as.factor(y))
    } else {
      stop("Incompatible X type. y must be of type factor or numeric.")
    }
    
    numClass <- length(ppForest$Levels)
    classCt <- cumsum(table(y))
    if (stratify) {
      Cindex <- vector("list", numClass)
      for (m in 1L:numClass) {
        Cindex[[m]] <- which(y == m)
      }
    }
    
  }else{
    ppForest <- list(Levels = NULL, params = NULL,trees = NULL);
    y=c(y)
  }

  Levels=ppForest$Levels
  n=length(y);p=ncol(X);
  ppForest$params=list(n=n,p=p,FUN=FUN,paramList=paramList,catMap=catMap,MaxLayer=MaxLayer,
                       minleaf=minleaf,ntrees=ntrees,method=method,numOOB=numOOB,storeOOB = storeOOB,
                       rotate=rotate,replacement=replacement,stratify=stratify,numNode=numNode,
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
  
  # rotate the data?
  if (rotate) {
    # if p > 1000 then rotate only a random subset of 1000 of the dimensions
    if (p > 1000L) {
      rotmat <- RandRot(1000L)
      rotdims <- sample.int(p, 1000L)
      X[, rotdims] <- X[, rotdims] %*% rotmat
    } else {
      rotmat <- RandRot(p)
      X[] <- X %*% rotmat
    }
  }

  
  PPtree=function(seed,...){
    set.seed(seed)
    ppForestT=list()

    TDindx <- seq(n)
    if (replacement) {
      go <- TRUE
      while (go) {
        # make sure each class is represented in proportion to classes in initial dataset
        if (stratify&(method!='r')) {
          if (classCt[1L] != 0L) {
            TDindx[1:classCt[1L]] <- sample(Cindex[[1L]], classCt[1L], replace = TRUE)
          }
          for (z in 2:numClass) {
            if (classCt[z - 1L] != classCt[z]) {
              TDindx[(classCt[z - 1L] + 1L):classCt[z]] <- sample(Cindex[[z]], classCt[z] - classCt[z - 1L], replace = TRUE)
            }
          }
        } else {
          TDindx <- sample(1L:n, n, replace = TRUE)
        }
        go <- all(1L:n %in% TDindx)
      }
    } else {
      TDindx <- sample.int(1:n, n-numOOB, replace = FALSE)
    }
    
    
    ppForestT=c(ppForestT,ppCART(X[TDindx,],y[TDindx],method,FUN,paramList,catMap,MaxLayer,
                                    numNode,minleaf,weights,Levels));
    
    oobErr=1
    if ((numOOB>0)&storeOOB){
      NTD = setdiff(1:n,TDindx);
      tree_output = ppCARTPredict(X[NTD,],ppForestT);
      
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
    #set.seed(seed)
    ppForestT = foreach::foreach(co=seq(length(chunks)),.combine=list,.multicombine = TRUE,
                                 .packages = "ppRF1",.noexport ="ppForest")%dopar%{
                                   lapply(seed+chunks[[co]], PPtree)
                                 }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
    for (i in 1:numCores) {
      ppForest$trees=c(ppForest$trees,ppForestT[[i]]) 
    }
  } else {
    # Use just one core.
    ppForest$trees <- lapply(seed+(1:ntrees), PPtree)
  }
  
  return(ppForest)
}