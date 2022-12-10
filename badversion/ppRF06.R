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
#' @param       weights      : a vector of values which weigh the samples when considering a split
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
#weights=1;numCores=0L;seed=220924;MaxLayer=Inf;numNode=Inf;Xcat=0;
#Xscale=c("Min-max","Quantile","No")[1];FUNDir=getwd();TreeRandRotate=F
#Xcat=c(NULL,0)[1]
ppRF = function(X,y,method='g-classification',NodeRotateFun="RandMatPPR",FunDir=getwd(),paramList=NULL,
                  catLabel=NULL,Xcat=0L,MinLeaf=ifelse(method=='regression',5,1),storeOOB = FALSE,
                  TreeRandRotate = FALSE,replacement = TRUE, stratify = TRUE,
                  ntrees=500,numOOB=c(0,ceiling(.368*length(y)))[2],
                  weights=1,numCores=0L,seed=220924,MaxLayer=Inf,numNode=Inf,
                  Xscale=c("Min-max","Quantile","No")[1])
{
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
  
  #ppForest$params=list(data=list(n=n,p=p,Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel),
  #                     tree=list(MaxLayer=MaxLayer,minleaf=minleaf,ntrees=ntrees,numNode=numNode),
  #                     rotation=list(rotate=rotate,FUN=FUN,FUNDir=FUNDir,paramList=paramList),
  #                     forest=list(numOOB=numOOB,storeOOB = storeOOB,replacement=replacement,stratify=stratify,weights=weights),
  #                     numCores=numCores,seed=seed)
  ppForest$data=list(n=n,p=p,Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel);
  ppForest$tree=list(MaxLayer=MaxLayer,MinLeaf=MinLeaf,ntrees=ntrees,numNode=numNode,TreeRandRotate=TreeRandRotate);
  ppForest$forest=list(numOOB=numOOB,storeOOB = storeOOB,replacement=replacement,stratify=stratify,
                       weights=weights,numCores=numCores)
  
  PPtree=function(seed,...){
    set.seed(seed)
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
    
    
    ppForestT=c(ppForestT,ppCART(X[TDindx,],y[TDindx],method,NodeRotateFun,FunDir,paramList,catLabel,Xcat=0L,MaxLayer,
                                   numNode,MinLeaf,weights,Levels,Xscale="No",TreeRandRotate));
    
    oobErr=1
    NTD = setdiff(TDindx0,TDindx);
    if ((numOOB>0)&storeOOB){
      pred = predict_ppCART(X[NTD,],ppForestT);
      
      if(method!="regression"){
        oobErr=mean(pred!=Levels[y[NTD]]);
      }else{
        oobErr=mean((pred-y[NTD])^2);
      }
    }
    ppForestT=c(ppForestT,list(oobErr=oobErr,OOBindx=NTD))
    
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
    ppForestT = foreach::foreach(co=seq(length(chunks)),.combine=list,.multicombine = TRUE,
                                 .packages = "ppRF",.noexport ="ppForest")%dopar%{
                                   lapply(seed+chunks[[co]], PPtree)
                                 }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
    #do.call(rbind.fill,list1)
    #do.call("c",ppTree)
    ppForest$ppTrees=NULL
    for (i in 1:numCores) {
      ppForest$ppTrees=c(ppForest$ppTrees,ppForestT[[i]]) 
    }
  } else {
    # Use just one core.
    ppForest$ppTrees <- lapply(seed+(1:ntrees), PPtree)
  }
  
  return(ppForest)
}
