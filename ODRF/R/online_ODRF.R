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
#' @import Rcpp
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @importFrom stats na.action
#' 
#' @aliases online.ODRF
#' @rdname online.ODRF
#' @method online ODRF
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
#ppForest=ppForest0
#X=Xnew
#y=ynew
online.ODRF = function(ppForest,data,weights=NULL)
{
  weights0=weights
  ppTrees=ppForest$ppTrees
  Levels=ppForest$Levels
  method=ppForest$method
  NodeRotateFun=ppForest$NodeRotateFun
  Call=ppForest$call
  Terms=ppForest$terms
  paramList=ppForest$paramList
  
  ppForest=ppForest[c(9,10,11)]
  ppForestVar=c(names(ppForest$data),names(ppForest$tree),names(ppForest$forest))
  ppForest=do.call("c",ppForest)
  
  for(v in seq(length(ppForestVar))){
    assign(ppForestVar[v], ppForest[[v]])
  }
  rm(ppForest)
  
  vars=all.vars(Terms)
  # address na values.
  if (any(is.na(data))) {
    data=na.action(data.frame(data))
    warning("NA values exist in data matrix")
  }
  y= data[,vars[1]]
  X= data[,vars[-1]]
  X=as.matrix(X) 
  
  if(!is.null(subset))
    X=X[subset,]
  #if(!is.null(weights))
  #  X <- X * matrix(weights0,length(y),ncol(X))
  weights=weights0
  
  
  ppForest <- list(call=Call,terms=Terms,method=method,Levels = NULL,
                   NodeRotateFun=NodeRotateFun,paramList=paramList,oobErr=NULL,oobConfusionMat=NULL)
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
  
  ppForest$data=list(subset=subset,weights=weights,na.action=na.action,n=n,p=p,varName=varName,
                     Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel);
  ppForest$tree=list(FunDir=FunDir,MaxDepth=MaxDepth,MinLeaf=MinLeaf,numNode=numNode,TreeRandRotate=TreeRandRotate);
  ppForest$forest=list(ntrees=ntrees,numOOB=numOOB,storeOOB = storeOOB,replacement=replacement,stratify=stratify,
                       numCores=numCores,seed=seed)
  
  
  PPtree=function(itree,...){
    ppTree=ppTrees[[itree]] #[1:7]
    if ((numOOB>0)&storeOOB){ppTree=ppTree[-(length(ppTree)-c(2,1,0))]}
    class(ppTree)="ODT"
    set.seed(seed+itree)
    
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

    data=data.frame(y[TDindx],X[TDindx,])
    colnames(data)=vars
    weights1=weights[TDindx]
    ppForestT=online(ppTree,data,weights1);
    
    if ((numOOB>0)&storeOOB){
      oobErr=1
      NTD = setdiff(TDindx0,TDindx);
      pred = predict(ppForestT,X[NTD,]);
      
      if(method!="regression"){
        oobErr=mean(pred!=Levels[y[NTD]]);
      }else{
        oobErr=mean((pred-y[NTD])^2);
      }
      
      ppForestT=c(ppForestT,list(oobErr=oobErr,oobIndex=NTD,oobPred=pred))
    }
    
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
    #library("ODRF1")
    #library(foreach)
    #foreach::registerDoSEQ()
    cl <- parallel::makeCluster(numCores,type = ifelse(.Platform$OS.type == "windows","PSOCK","FORK"))
    chunks <- parallel::clusterSplit(cl,seq(ntrees))
    doParallel::registerDoParallel(cl,numCores)
    #set.seed(seed)
    ppForestT = foreach::foreach(icore=seq(length(chunks)),.combine=list,.multicombine = TRUE,
                                 .packages = "ODRF",.noexport ="ppForest")%dopar%{
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
  
  ####################################
  if ((numOOB>0)&storeOOB){
    oobVotes=matrix(NA,n,ntrees)
    for (t in 1:ntrees) {
      oobVotes[ppForest$ppTrees[[t]]$oobIndex,t]=ppForest$ppTrees[[t]]$oobPred
    }
    idx=which(rowSums(is.na(oobVotes))<ntrees) 
    oobVotes=oobVotes[idx,,drop = FALSE]
    yy=y[idx]
    
    if(method!="regression"){
      ny=length(yy)
      nC=numClass
      weights=rep(1,ny*ntrees)
      Votes=factor(c(t(oobVotes)),levels =Levels)
      Votes=as.integer(Votes)+nC*rep(0:(ny-1),rep(ntrees,ny));
      Votes=aggregate(c(rep(0,ny*nC),weights), by=list(c(1:(ny*nC),Votes)),sum)[,2];
      
      prob=matrix(Votes,ny,nC,byrow = TRUE);
      pred=max.col(prob) ## "random"
      oobPred=Levels[pred]
      ppForest$oobErr=mean(Levels[yy]!=oobPred)
      
      #oobPred=rep(NA,noob)
      #for (i in 1:noob) {
      #  oobTable = table(oobVotes[i,])
      #  oobPred[i]=names(oobTable)[which.max(oobTable)];
      #}
      
      #oobErr=mean(oobPred!=Levels[y[idx]]);
      XConfusionMat=table(oobPred,yy)
      class_error=(rowSums(XConfusionMat)-diag(XConfusionMat))/rowSums(XConfusionMat)
      XConfusionMat=cbind(XConfusionMat,class_error)
      ppForest$oobConfusionMat=XConfusionMat
    }else{
      oobPred=rowMeans(oobVotes,na.rm = TRUE);
      ppForest$oobErr=mean((oobPred-yy)^2)/mean((yy-mean(y))^2);
    }
  }
  
  #class(ppTree) <- append(class(ppTree),"ODRF")
  class(ppForest) <- "ODRF"
  return(ppForest)
}
