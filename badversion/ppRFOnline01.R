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
#'                           'g-classification' : gini impurity index (classification)
#'                           'i-classification' : information gain (classification)
#'                           'regression' : squared error (regression)
#'
#' @param       minparent    : the minimum amount of samples in an impure node for it to be considered for splitting
#'
#' @param       minleaf      : the minimum amount of samples in a leaf
#'
#' @param       Xweights      : a vector of values which weigh the samples when considering a split
#'
#' @param       nvartosample : the number of (randomly selected) variables to consider at each node 
#'
#' @return ppForest
#'
#' @useDynLib ppRF
#' @import Rcpp
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @importFrom stats na.action
#' 
#' @export
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
## ##@importFrom foreach foreach `%dopar%`
#method='g-classification';NodeRotateFun="RandMatPPR02";FunDir=getwd();paramList=NULL;catLabel=NULL;
#MinLeaf=ifelse(method=='regression',5,1);rotate = FALSE;replacement = TRUE;stratify = TRUE;
#ntrees=500;numOOB=c(0,ceiling(length(y)/5))[2];storeOOB = FALSE;
#Xweights=1;numCores=0L;seed=220924;MaxDepth=Inf;numNode=Inf;Xcat=0;
#Xscale=c("Min-max","Quantile","No")[1];FUNDir=getwd();TreeRandRotate=F
#Xcat=c(NULL,0)[1]
ppRFOnline = function(ppForest,X,y)
{
  #ppForest=ppForest0
  #X=Xnew
  #y=ynew
  ppTrees=ppForest$ppTrees
  Levels=ppForest$Levels
  ppForest=ppForest[-c(1,3,8)]
  ppForestVar=c("method","NodeRotateFun",names(ppForest$data),names(ppForest$tree),names(ppForest$forest))
  ppForest=do.call("c",ppForest)
  
  for(v in seq(length(ppForestVar))){
    assign(ppForestVar[v], ppForest[[v]])
  }
  rm(ppForest)

  
  ppForest <- list(call=match.call(),method=method,Levels = NULL,NodeRotateFun=NodeRotateFun)
  if(method!="regression"){
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
  }
  
  Levels=ppForest$Levels
  n=length(y);p=ncol(X);
  
  
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
  
  ppForest$data=list(n=n,p=p,Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel);
  ppForest$tree=list(FunDir=FunDir,paramList=paramList,MaxDepth=MaxDepth,MinLeaf=MinLeaf,ntrees=ntrees,
                     numNode=numNode,TreeRandRotate=TreeRandRotate);
  ppForest$forest=list(numOOB=numOOB,storeOOB = storeOOB,replacement=replacement,stratify=stratify,
                       Xweights=Xweights,numCores=numCores,seed=seed)
  
  PPtree=function(itree,...){
    ppTree=ppTrees[[itree]]
    set.seed(seed+itree)
    ppForestT=list()
    
    TDindx0 <- seq(n)
    TDindx<-TDindx0
    if (replacement) {
      go <- TRUE
      while (go) {
        # make sure each class is represented in proportion to classes in initial dataset
        if (stratify&(method!='regression')) {
          if (classCt[1L] != 0L) {
            TDindx[1:classCt[1L]] <- sample(Cindex[[1L]], classCt[1L], replace = TRUE)
          }
          for (z in 2:numClass) {
            if (classCt[z - 1L] != classCt[z]) {
              TDindx[(classCt[z - 1L] + 1L):classCt[z]] <- sample(Cindex[[z]], classCt[z] - classCt[z - 1L], replace = TRUE)
            }
          }
        } else {
          TDindx <- sample(TDindx0, n, replace = TRUE)
        }
        go <- all(TDindx0 %in% TDindx)
      }
    } else {
      TDindx <- sample.int(TDindx0, n-numOOB, replace = FALSE)
    }
    
    
    ppForestT=c(ppForestT,ppCARTOnline(ppTree[-length(ppTree)+c(0,1)],X[TDindx,],y[TDindx]));
    
    oobErr=1
    NTD = setdiff(TDindx0,TDindx);
    if ((numOOB>0)&storeOOB){
      pred = predict_ppCART(ppForestT,X[NTD,]);
      
      if(method!="regression"){
        oobErr=mean(pred!=Levels[y[NTD]]);
      }else{
        oobErr=mean((pred-y[NTD])^2);
      }
    }
    ppForestT=c(ppForestT,list(oobErr=oobErr,oobIndex=NTD))
    
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
    #library(foreach)
    #foreach::registerDoSEQ()
    cl <- parallel::makeCluster(numCores,type = ifelse(.Platform$OS.type == "windows","PSOCK","FORK"))
    chunks <- parallel::clusterSplit(cl,seq(ntrees))
    doParallel::registerDoParallel(cl,numCores)
    #set.seed(seed)
    ppForestT = foreach::foreach(icore=seq(length(chunks)),.combine=list,.multicombine = TRUE,
                                 .packages = "ppRF",.noexport ="ppForest")%dopar%{
                                   lapply(chunks[[icore]], PPtree)
                                 }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
    #do.call(rbind.fill,list1)
    ppForest$ppTrees=do.call("c",ppForestT)
    #ppForest$ppTrees=NULL
    #for (i in 1:numCores) {
    #  ppForest$ppTrees=c(ppForest$ppTrees,ppForestT[[i]]) 
    #}
  } else {
    # Use just one core.
    ppForest$ppTrees <- lapply(1:ntrees, PPtree)
  }
  
  return(ppForest)
}
