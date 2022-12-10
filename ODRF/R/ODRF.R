#' ODRF forest Generator
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
#' @useDynLib ODRF
#' @import Rcpp
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @importFrom stats model.frame model.extract model.matrix na.fail
#' 
#' @export
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
## ##@importFrom foreach foreach `%dopar%`
#formula=y~.;data=data.frame(X,y=y);subset=NULL;weights=NULL;na.action=na.fail;method='i-classification';
#NodeRotateFun="RotMatPPO";FunDir=getwd();paramList=NULL;
#catLabel=NULL;Xcat=0L;MinLeaf=ifelse(method=='regression',5,1);storeOOB = TRUE;
#TreeRandRotate = FALSE;replacement = TRUE; stratify = TRUE;
#ntrees=100;numOOB=c(0,ceiling(length(y)/3))[2];
#numCores=0L;seed=220924;MaxDepth=Inf;numNode=Inf;
#Xscale=c("Min-max","Quantile","No")[1]
#Xcat=c(NULL,0)[1]
ODRF = function(formula,data,subset=NULL,weights=NULL,na.action=na.fail,method='i-classification',
                NodeRotateFun="RotMatPPO",FunDir=getwd(),paramList=NULL,catLabel=NULL,Xcat=0L,
                MinLeaf=ifelse(method=='regression',5,1),storeOOB = TRUE,TreeRandRotate = FALSE,
                replacement = TRUE, stratify = TRUE,ntrees=100,numOOB=c(0,ceiling(length(y)/3))[2],
                numCores=0L,seed=220924,MaxDepth=Inf,numNode=Inf,Xscale=c("Min-max","Quantile","No")[1])
{
  # address na values.
  if (any(is.na(data))) {
    warning("NA values exist in data matrix")
  }
  
  Call<-match.call()
  indx<-match(c("formula", "data", "subset", "na.action"),names(Call),nomatch=0L)#"weights",
  if (indx[[1]]==0)
    stop("a 'formula' argument is required")
  if (indx[[2]]==0){
    stop("a 'data' argument is required")
    #data <- environment(formula)
    #data <- data.frame(eval(formula[[2]]),eval(formula[[3]]))
    #varName=colnames(eval(formula[[3]]))
    #colnames(data)=c(as.character(formula[[2]]),varName)
    #Call$formula=
  }
  
  #Call=match.call(expand.dots = FALSE)
  #indx<-match(c("formula", "data", "subset", "weights", "na.action"),names(Call),nomatch=0L)
  temp<-Call[c(1L,indx)]
  temp[[1L]]<-quote(stats::model.frame)
  temp$drop.unused.levels <- TRUE
  temp <- eval(temp, parent.frame())
  Terms<-attr(temp, "terms")
  
  varName=setdiff(colnames(data),as.character(formula[[2]]))
  y <- model.extract(temp, "response")
  X <- model.matrix(Terms, temp)
  int <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  if (int > 0)
    X <- X[, -int, drop = FALSE]
  #if(!is.null(weights))
  #  X <- X * matrix(weights,length(y),ncol(X))
  colnames(X)=varName
  
  
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
  
  ppForest$data=list(subset=subset,weights=weights,na.action=na.action,n=n,p=p,varName=varName,
                     Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel);
  ppForest$tree=list(FunDir=FunDir,MaxDepth=MaxDepth,MinLeaf=MinLeaf,numNode=numNode,TreeRandRotate=TreeRandRotate);
  ppForest$forest=list(ntrees=ntrees,numOOB=numOOB,storeOOB = storeOOB,replacement=replacement,stratify=stratify,
                       numCores=numCores,seed=seed)
  
  #Weights=weights
  vars=all.vars(Terms)
  PPtree=function(itree,...){
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
      TDindx <- sample(TDindx0, n-numOOB, replace = FALSE)
    }
    
    
    data=data.frame(y[TDindx],X[TDindx,])
    colnames(data)=vars
    weights1=weights[TDindx]
    #weights1=rep(1,nrow(data))
    #if(is.null(Call$data)){
    #  formula=y~.
    #}else{
    #  colnames(data)=all.vars(Terms)
    #}
    paramList$weights=weights1
    ppForestT=ODT(formula,data,subset=NULL,weights=NULL,na.action=NULL,method,NodeRotateFun,FunDir,
                  paramList,catLabel,Xcat=0L,MaxDepth,numNode,MinLeaf,Levels,Xscale="No",TreeRandRotate);
    
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
      tree_weights=rep(1,ny*ntrees)
      Votes=factor(c(t(oobVotes)),levels =Levels)
      Votes=as.integer(Votes)+nC*rep(0:(ny-1),rep(ntrees,ny));
      Votes=aggregate(c(rep(0,ny*nC),tree_weights), by=list(c(1:(ny*nC),Votes)),sum)[,2];
      
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
