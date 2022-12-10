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
#' @aliases prune.ODRF
#' @rdname prune.ODRF
#' @method prune ODRF
#' @export
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
## ##@importFrom foreach foreach `%dopar%` X y
#method='g-classification';NodeRotateFun="RandMatPPR02";FunDir=getwd();paramList=NULL;catLabel=NULL;
#MinLeaf=ifelse(method=='regression',5,1);rotate = FALSE;replacement = TRUE;stratify = TRUE;
#ntrees=500;numOOB=c(0,ceiling(length(y)/5))[2];storeOOB = FALSE;
#Xweights=1;numCores=0L;seed=220924;MaxDepth=Inf;numNode=Inf;Xcat=0;
#Xscale=c("Min-max","Quantile","No")[1];FUNDir=getwd();TreeRandRotate=F
#Xcat=c(NULL,0)[1]
#ppForest=ppForest0
#Xnew=X1
#ynew=y1
prune.ODRF = function(ppForest,data,weights=NULL,MaxDepth=NULL,useOOB=TRUE)
{
  ppTrees=ppForest$ppTrees
  method=ppForest$method
  numOOB=ppForest$forest$numOOB
  storeOOB=ppForest$forest$storeOOB
  seed=ppForest$forest$seed
  replacement=ppForest$forest$replacement
  stratify=ppForest$forest$stratify
  numCores=ppForest$forest$numCores
  
  if ((!storeOOB)&useOOB) {
    stop("out-of-bag indices for each tree are not stored. ODRF must be called with storeOOB = TRUE.")
  }
  
  vars=all.vars(ppForest$terms)
  # address na values.
  if (any(is.na(data))) {
    data=ppForest$data$na.action(data.frame(data))
    warning("NA values exist in data matrix")
  }
  ynew= data[,vars[1]]
  Xnew= data[,vars[-1]]
  Xnew=as.matrix(Xnew) 
  
  if(!is.null(ppForest$data$subset))
    Xnew=Xnew[ppForest$data$subset,]
  
  p=ncol(Xnew)
  n=nrow(Xnew)
  nC=length(ppForest$Levels)
  numClass=nC
  ntrees=length(ppTrees);
  
  classCt <- cumsum(table(ynew))
  if (stratify) {
    Cindex <- vector("list", numClass)
    for (m in 1L:numClass) {
      Cindex[[m]] <- which(ynew == ppForest$Levels[m])
    }
  }

  Xcat=ppForest$data$Xcat
  catLabel=ppForest$data$catLabel
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
    
    Xnew=cbind(Xnew1,Xnew[,-Xcat]) 
    p=ncol(Xnew)
    numCat=length(unlist(catLabel))
    rm(Xnew1)
    rm(Xnewj)
  }
  
  #Variable scaling.
  if(ppForest$data$Xscale!="No"){
    indp=(sum(numCat)+1):p
    Xnew[,indp]=(Xnew[,indp]-matrix(ppForest$data$minCol,n,length(indp),byrow = T))/
      matrix(ppForest$data$maxminCol,n,length(indp),byrow = T)
  }

  
  PPtree=function(itree,...){
    ppTree=ppTrees[[itree]] #[1:7]
    class(ppTree)="ODT"
    set.seed(seed+itree)

    if(useOOB){
      data=data.frame(y=ynew[ppTree$oobIndex],Xnew[ppTree$oobIndex,])
      colnames(data)=vars
      weights1=weights[ppTree$oobIndex]
      ppForestT=prune(ppTree,data,weights1,MaxDepth)#[seq(7)]
      ppTree=ppForestT[-length(ppForestT)]
    }else{
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
      
      if ((numOOB>0)&storeOOB){ppTree=ppTree[-(length(ppTree)-c(2,1,0))]}
      data=data.frame(y=ynew[TDindx],Xnew[TDindx,])
      colnames(data)=vars
      weights1=weights[TDindx]
      ppTree=prune(ppTree,data,weights1,MaxDepth)
      ppTree=ppTree[-length(ppTree)]
      class(ppTree)="ODT"
  
      if ((numOOB>0)&storeOOB){
        oobErr=1
        #if(useOOB){
        #  NTD = ppTree$oobIndex
        #}else{
          NTD = setdiff(TDindx0,TDindx);
        #}
        pred = predict(ppTree,Xnew[NTD,]);

        if(method!="regression"){
          oobErr=mean(pred!=ynew[NTD]);
        }else{
          oobErr=mean((pred-ynew[NTD])^2);
        }
        
        ppTree=c(ppTree,list(oobErr=oobErr,oobIndex=NTD,oobPred=pred))
      }
    }
    
    return(ppTree)
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
  if ((numOOB>0)&storeOOB&(!useOOB)){
    oobVotes=matrix(NA,n,ntrees)
    for (t in 1:ntrees) {
      oobVotes[ppForest$ppTrees[[t]]$oobIndex,t]=ppForest$ppTrees[[t]]$oobPred
    }
    idx=which(rowSums(is.na(oobVotes))<ntrees) 
    oobVotes=oobVotes[idx,,drop = FALSE]
    
    if(method!="regression"){
      yy=ynew[idx]
      ny=length(yy)
      nC=numClass
      weights=rep(1,ny*ntrees)
      Votes=factor(c(t(oobVotes)),levels =ppForest$Levels)
      Votes=as.integer(Votes)+nC*rep(0:(ny-1),rep(ntrees,ny));
      Votes=aggregate(c(rep(0,ny*nC),weights), by=list(c(1:(ny*nC),Votes)),sum)[,2];
      
      prob=matrix(Votes,ny,nC,byrow = TRUE);
      pred=max.col(prob) ## "random"
      oobPred=ppForest$Levels[pred]
      ppForest$oobErr=mean(yy!=oobPred)
      
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
  class(ppForest) <- c('ODRF',"prune.ODRF")
  return(ppForest)
}
