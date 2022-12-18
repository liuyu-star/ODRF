#' Classification and Regression with Oblique Decision Random Forest
#' 
#' ODRF is the assembly of multiple ODTs, which can effectively reduce the overfitting of individual ODTs and improve the accuracy of classification and regression.
#'
#' @param formula Object of class \code{formula} with a response but no interaction terms describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} in which to interpret the variables named in the formula.If data is missing it is obtained from the current environment by \code{formula}.
#' @param type The criterion used for splitting the nodes. g-classification': gini impurity index(default) and i-classification': information gain for classification; 'regression': mean square error for regression.
#' y is a factor then is a classification.
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictor variables in the split node. Default is "RotMatPPO" with model="PPR". (see \code{\link{RotMatPPO}}) 
#' Users can define this function, for details see \code{\link{RotMatMake}}.
#' @param FunDir The path to the \code{function} of the user-defined \code{NodeRotateFun}. (default current Workspace)
#' @param paramList Parameters in a named list to be used by \code{NodeRotateFun}.If left unchanged, default values will be populated, for details see \code{\link[ODRF]{defaults}}.
#' @param ntrees The number of trees in the forest. (default 100)
#' @param storeOOB if TRUE then the samples omitted during the creation of a tree are stored as part of the tree.
#' @param replacement if TRUE then n samples are chosen, with replacement, from training data. (default TRUE)
#' @param stratify if TRUE then class sample proportions are maintained during the random sampling. Ignored if replacement = FALSE. (default TRUE)
#' @param numOOB  Ratio of 'out-of-bag'.
#' @param parallel Parallel computing or not. (default TRUE)
#' @param numCores Number of cores to be used for parallel computing. (default Inf)
#' @param seed Random seeds in order to reproduce results.
#' @param MaxDepth The maximum depth of the tree (default \code{Inf}).
#' @param numNode The number of nodes used by the tree (default \code{Inf}).
#' @param MinLeaf Minimal node size. Default 1 for classification, 5 for regression.
#' @param subset An index vector indicating which rows should be used. (NOTE: If given, this argument must be named.)
#' @param weights A vector of length same as \code{data} that are positive weights.(default NULL)
#' @param na.action A function to specify the action to be taken if NAs are found. (NOTE: If given, this argument must be named.)
#' @param catLabel A category labels of class \code{list} in prediction variables. (for details see Details and Examples)
#' @param Xcat A class \code{vector} is used to indicate which variables are class variables. The default Xcat=0 means that no special treatment is given to category variables. 
#' When Xcat=NULL, the variable x that satisfies the condition \code{(length(unique(x))<10) & (n>20)} is judged to be a category variable #' @param Xscale Predictor variable standardization methods." Min-max", "Quantile", "No" denote Min-max transformation, Quantile transformation and No transformation (default "Min-max"), respectively.
#' @param TreeRandRotate If or not to randomly rotate the Training data before building the tree. (default FALSE)
#' @param ... optional parameters to be passed to the low level function.
#' 
#' @return An object of class ODT, which is a list with the following components:
#' \itemize{
#' \item{\code{call}: The original call to ODT.}
#' \item{\code{terms}: An object of class \code{c("terms", "formula")} (see \code{\link{terms.object}}) summarizing the formula. Used by various methods, but typically not of direct relevance to users.}
#' \item{\code{ppTrees}: Each tree used to build the forest. \itemize{
#' \item{\code{oobErr}: 'out-of-bag' error for tree, classification error rate for classification or mean square error for regression.}
#' \item{\code{oobIndex}: Which training data to use as 'out-of-bag'?}
#' \item{\code{oobPred}: Predicted value for 'out-of-bag'.}
#' \item{\code{other}: For other tree-related values see \code{\link{ODT}}.}
#' }}
#' \item{\code{oobErr}: 'out-of-bag' error for forest, classification error rate for classification or mean square error for regression.}
#' \item{\code{oobConfusionMat}: 'out-of-bag' confusion matrix for forest.}
#' \item{\code{type}, \code{Levels} and \code{NodeRotateFun} are important parameters for building the tree.}
#' \item{\code{paramList}: Parameters in a named list to be used by \code{NodeRotateFun}.}
#' \item{\code{data}: The list of data related parameters used to build the tree.}
#' \item{\code{tree}: The list of tree related parameters used to build the tree.}
#' \item{\code{forest}: The list of forest related parameters used to build the tree.}
#' }
#'
#' @seealso \code{\link{predict.OORF}} \code{\link{print.OORF}} \code{\link{ODRF.error}} \code{\link{VarImp}}
#' 
#' @author YU Liu and Yingcun Xia
#' @references \itemize{
#' \item{Zhan H, Liu Y, Xia Y. Consistency of The Oblique Decision Tree and Its Random Forest[J]. arXiv preprint arXiv:2211.12653, 2022.}
#' \item{Tomita T M, Browne J, Shen C, et al. Sparse projection oblique randomer forests[J]. Journal of machine learning research, 2020, 21(104).}
#' }
#' @examples
#' #Classification with Oblique Decision Tree
#' data(seeds)
#' set.seed(221212)
#' train = sample(1:209,100)
#' train_data = data.frame(seeds[train,])
#' test_data = data.frame(seeds[-train,])
#' 
#' tree = ODRF(varieties_of_wheat~.,train_data,type='i-classification')
#' pred <- predict(tree,test_data[,-8])
#' #estimation error
#' (mean(pred!=test_data[,8]))
#' 
#' #Regression with Oblique Decision Tree
#' data(body_fat)
#' set.seed(221212)
#' train = sample(1:252,100)
#' train_data = data.frame(body_fat[train,])
#' test_data = data.frame(body_fat[-train,])
#' 
#' tree = ODRF(Density~.,train_data,type='regression')
#' pred <- predict(tree,test_data[,-1])
#' #estimation error
#' mean((pred-test_data[,1])^2)
#' 
#' 
#' ### Train ODT on one-of-K encoded categorical data ###
#' Xcol1=sample(c("A","B","C"),100,replace = TRUE)
#' Xcol2=sample(c("1","2","3","4","5"),100,replace = TRUE)
#' Xcon=matrix(rnorm(100*3),100,3)
#' X=data.frame(Xcol1,Xcol2,Xcon)
#' Xcat=c(1,2)
#' catLabel=NULL
#' y=as.factor(sample(c(0,1),100,replace = TRUE))
#' tree = ODRF(y~X,type='g-classification')
#' 
#' numCat <- apply(X[,Xcat,drop = FALSE], 2, function(x) length(unique(x)))
#' X1 <- matrix(0, nrow = nrow(X), ncol = sum(numCat)) # initialize training data matrix X
#' catLabel <- vector("list", length(Xcat))
#' names(catLabel)<- colnames(X)[Xcat]
#' col.idx <- 0L
#' # one-of-K encode each categorical feature and store in X
#' for (j in 1:length(Xcat)) {
#'   catMap <- (col.idx + 1L):(col.idx + numCat[j])
#'   # convert categorical feature to K dummy variables
#'   catLabel[[j]]=levels(as.factor(X[,Xcat[j]]))
#'   X1[, catMap] <- (matrix(X[,Xcat[j]],nrow(X),numCat[j])==matrix(catLabel[[j]],nrow(X),numCat[j],byrow = TRUE))+0
#'   col.idx <- col.idx + numCat[j]
#' }
#' X=cbind(X1,X[,-Xcat])
#' 
#' #Print the result after processing of category variables
#' X
#' #  1 2 3 4 5 6 7 8          X1         X2          X3
#' #1 0 1 0 0 1 0 0 0 -0.81003483  0.7900958 -1.94504333
#' #2 0 0 1 0 0 0 0 1 -0.02528851 -0.5143964 -0.18628226
#' #3 1 0 0 1 0 0 0 0  1.15532067  2.0236020  1.02942500
#' #4 1 0 0 0 0 1 0 0  1.18598589  1.0594630  0.42990019
#' #5 1 0 0 1 0 0 0 0 -0.21695438  1.5145973  0.09316665
#' #6 0 0 1 0 0 0 0 1 -1.11507717 -0.5775602  0.09918911
#' catLabel
#' #$Xcol1
#' #[1] "A" "B" "C"
#' #$Xcol2
#' #[1] "1" "2" "3" "4" "5"
#' 
#' @useDynLib ODRF
#' @import Rcpp doParallel foreach
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @export
ODRF = function(formula,data=NULL,type=NULL,NodeRotateFun="RotMatPPO",FunDir=getwd(),paramList=NULL,
                ntrees=100,storeOOB = TRUE,replacement = TRUE, stratify = TRUE,numOOB=1/3,parallel=TRUE,
                numCores=Inf,seed=220924, MaxDepth=Inf,numNode=Inf,MinLeaf=5,subset=NULL,weights=NULL,
                na.action=na.fail,catLabel=NULL,Xcat=0,Xscale="Min-max",TreeRandRotate=FALSE,...)
{
  Call<-match.call()
  indx<-match(c("formula", "data", "subset", "na.action"),names(Call),nomatch=0L)#, "weights"
  if (indx[[1]]==0)
    stop("A 'formula' argument is required")
  if (indx[[2]]==0){
    #stop("a 'data' argument is required")
    #data <- environment(formula)
    data <- data.frame(y=eval(formula[[2]]),eval(formula[[3]]))
    formula=y~.
    Call$formula=formula
    Call$data=quote(data)
  }
  # address na values.
  if (any(is.na(data))) {
    warning("NA values exist in data frame")
  }
  
  #Call=match.call(expand.dots = FALSE)
  indx<-match(c("formula", "data", "subset", "na.action"),names(Call),nomatch=0L)
  temp<-Call[c(1L,indx)]
  temp[[1L]]<-quote(stats::model.frame)
  temp$drop.unused.levels <- TRUE
  temp <- eval(temp)#, parent.frame()
  Terms<-attr(temp, "terms")
  #Terms<-stats::model.frame(temp)
  
  varName=setdiff(colnames(data),as.character(formula[[2]]))
  y <- model.extract(temp, "response")
  X <- model.matrix(Terms, temp)
  int <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  if (int > 0)
    X <- X[, -int, drop = FALSE]
  #if(!is.null(weights))
  #  X <- X * matrix(weights,length(y),ncol(X))
  colnames(X)=varName
  
  if(is.factor(y)&is.null(type)){
    type='i-classification'
    warning("You are creating a forest for classification")
  }
  if(is.numeric(y)&is.null(type)){
    type='regression'
    warning("You are creating a forest for regression")
  }
  if(is.factor(y)&type=='regression')
    stop(paste0("When ",formula[[2]]," is a factor type, type cannot take 'regression'."))
  if(MinLeaf==5)
    MinLeaf=ifelse(type=='regression',5,1)
  
  ppForest <- list(call=Call,terms=Terms,type=type,Levels = NA,
                   NodeRotateFun=NodeRotateFun,paramList=paramList,oobErr=NULL,oobConfusionMat=NULL)
  if(type!="regression"){
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
                       parallel=parallel,numCores=numCores,seed=seed)
  
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
        if (stratify&(type!='regression')) {
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
      TDindx <- sample(TDindx0, n*(1-numOOB), replace = FALSE)
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
    #paramList$weights=weights1
    ppForestT=ODT(formula,data,type,NodeRotateFun,FunDir,paramList,MaxDepth,numNode,MinLeaf,
                  Levels,subset=NULL,weights=weights1,na.action=NULL,catLabel,Xcat=0L,Xscale="No",TreeRandRotate);
    
    if ((numOOB>0)&storeOOB){
      oobErr=1
      NTD = setdiff(TDindx0,TDindx);
      pred = predict(ppForestT,X[NTD,]);
      
      if(type!="regression"){
        oobErr=mean(pred!=Levels[y[NTD]]);
      }else{
        oobErr=mean((pred-y[NTD])^2);
      }
      
      ppForestT=c(ppForestT,list(oobErr=oobErr,oobIndex=NTD,oobPred=pred))
    }
    
    return(ppForestT)
  }
  
  
  if (parallel) {
    #RNGkind("L'Ecuyer-CMRG")
    if (is.infinite(numCores)) {
      # Use all but 1 core if numCores=0.
      numCores <- parallel::detectCores()- 1L #logical = FALSE
    }
    numCores <- min(numCores, ntrees)
    gc()
    
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
    
    if(type!="regression"){
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
      XConfusionMat=table(oobPred,Levels[yy])
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
