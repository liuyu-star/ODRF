#' Classification and Regression with Oblique Decision Tree
#' 
#' ODT uses a linear combination of predictors as partitioning variables to create classification and regression tree. 
#'
#' @param formula Object of class \code{formula} with a response but no interaction terms describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} in which to interpret the variables named in the formula.If data is missing it is obtained from the current environment by \code{formula}.
#' @param type The criterion used for splitting the nodes. g-classification': gini impurity index(default) and i-classification': information gain for classification; 'regression': mean square error for regression.
#' y is a factor then is a classification.
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictor variables in the split node. Default is "RotMatPPO" with model="PPR". (see \code{\link{RotMatPPO}}) 
#' Users can define this function, for details see \code{\link{RotMatMake}}.
#' @param FunDir The path to the \code{function} of the user-defined \code{NodeRotateFun}. (default current Workspace)
#' @param paramList Parameters in a named list to be used by \code{NodeRotateFun}.If left unchanged, default values will be populated, for details see \code{\link[ODRF]{defaults}}.
#' @param MaxDepth The maximum depth of the tree (default \code{Inf}).
#' @param numNode The number of nodes used by the tree (default \code{Inf}).
#' @param MinLeaf Minimal node size. Default 1 for classification, 5 for regression.
#' @param Levels The category label of the response variable when \code{type} is not equal to 'regression'.
#' @param subset An index vector indicating which rows should be used. (NOTE: If given, this argument must be named.)
#' @param weights A vector of length same as \code{data} that are positive weights.
#' @param na.action A function to specify the action to be taken if NAs are found. (NOTE: If given, this argument must be named.)
#' @param catLabel A category labels of class \code{list} in prediction variables. (for details see Details and Examples)
#' @param Xcat A class \code{vector} is used to indicate which variables are class variables. The default Xcat=0 means that no special treatment is given to class variables.
#' @param Xscale Predictor variable standardization methods." Min-max", "Quantile", "No" denote Min-max transformation, Quantile transformation and No transformation, respectively. (default "Min-max")
#' @param TreeRandRotate If or not to randomly rotate the Training data before building the tree. (default FALSE)
#' @param ... optional parameters to be passed to the low level function.
#' 
#' @return An object of class ODT, which is a list with the following components:
#' \itemize{
#' \item{\code{call}: The original call to ODT.}
#' \item{\code{terms}: An object of class \code{c("terms", "formula")} (see \code{\link{terms.object}}) summarizing the formula. Used by various methods, but typically not of direct relevance to users.}
#' \item{\itemize{\code{structure}: A set of tree structure data records.
#' \item{\code{nodeRotaMat}: Record the split variables (first column), split node serial number (second column) and rotation direction (third column) for each node. (The first column and the third column are 0 means leaf nodes)}
#' \item{\code{nodeNumLabel}: Record each leaf node's category for classification or predicted value for regression (second column is data size). (Each column is 0 means it is not a leaf node)}
#' \item{\code{nodeCutValue}: Record the split point of each node. (0 means leaf nodes)}
#' \item{\code{nodeCutIndex}: Record the index values of the partitioning variables selected based on the partition criterion \code{type}.}
#' \item{\code{childNode}: Record the number of child nodes after each splitting.}
#' \item{\code{nodeDepth}: Record the depth of the tree where each node is located.}
#' }}
#' \item{\code{type}, \code{Levels} and \code{NodeRotateFun} are important parameters for building the tree.}
#' \item{\code{paramList}: Parameters in a named list to be used by \code{NodeRotateFun}.}
#' \item{\code{data}: The list of data related parameters used to build the tree.}
#' \item{\code{tree}: The list of tree related parameters used to build the tree.}
#' }
#'
#' @seealso \code{\link{predict.ODT}} \code{\link{print.ODT}} \code{\link{plot.ODT}} \code{\link{plot_ODT_depth}}
#' 
#' @author YU Liu and Yingcun Xia
#' @references Zhan H, Liu Y, Xia Y. Consistency of The Oblique Decision Tree and Its Random Forest[J]. arXiv preprint arXiv:2211.12653, 2022.
#'
#' @examples
#' library(ODRF)
#' 
#' #Classification with Oblique Decision Tree
#' data(seeds)
#' set.seed(221212)
#' train = sample(1:209,100)
#' train_data = data.frame(seeds[train,])
#' test_data = data.frame(seeds[-train,])
#' 
#' tree = ODT(varieties_of_wheat~.,train_data,type='i-classification')
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
#' tree = ODT(Density~.,train_data,type='regression')
#' pred <- predict(tree,test_data[,-1])
#' #estimation error
#' mean((pred-test_data[,1])^2)
#' 
#' @useDynLib ODRF
#' @import Rcpp
#' @importFrom stats model.frame model.extract model.matrix na.fail
#' @export
ODT=function(formula,data=NULL,type=NULL,NodeRotateFun="RotMatPPO",FunDir=getwd(),paramList=NULL,
             MaxDepth=Inf,numNode=Inf,MinLeaf=5,Levels=NULL,subset=NULL,weights=NULL,na.action=na.fail,
             catLabel=NULL,Xcat=0,Xscale="Min-max",TreeRandRotate=FALSE,...)
{
  Call<-match.call()
  indx<-match(c("formula", "data", "subset", "na.action"),names(Call),nomatch=0L)#, "weights"
  if (indx[[1]]==0)
    stop("a 'formula' argument is required")
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
  
  #weights=c(weights,paramList$weights)
  if(!is.null(weights))
    X <- X * matrix(weights,length(y),ncol(X))
  colnames(X)=varName
  
  
  if(is.factor(y)&is.null(type)){
    type='i-classification'
    warning("You are creating a tree for classification")
  }
  if(is.numeric(y)&is.null(type)){
    type='regression'
    warning("You are creating a tree for regression")
  }
  if(is.factor(y)&type=='regression')
    stop(paste0("When ",formula[[2]]," is a factor type, type cannot take 'regression'."))
  if(MinLeaf==5)
    MinLeaf=ifelse(type=='regression',5,1)
  
  if(!NodeRotateFun%in%ls("package:ODRF")){
    source(paste0(FunDir,"/",NodeRotateFun,".R"))
  }
  
  FUN <- match.fun(NodeRotateFun, descend = TRUE)
  type0=strsplit(type,split = "")[[1]][1]
  
  
  n = length(y);
  p = ncol(X);
  
  if(type!="regression"){
    if(is.null(Levels)){
      Levels=levels(as.factor(y))
    }
    if (!is.integer(y)) {
      y <- as.integer(as.factor(y));
    }
    maxLabel = length(Levels); 
  }else{
    y=c(y)
    maxLabel= 0;
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
  
  # rotate the data?
  rotdims=NULL;rotmat=NULL
  if (TreeRandRotate) {
    # if p > 1000 then rotate only a random subset of 1000 of the dimensions
    if (p > 1000L) {
      rotmat <- RandRot(1000L)
      rotdims <- sample.int(p, 1000L)
      X[, rotdims] <- X[, rotdims] %*% rotmat
    } else {
      rotdims=1:p
      rotmat <- RandRot(p)
      X <- X %*% rotmat
    }
  }
  
  paramList = ODRF:::defaults(paramList,p,catLabel,NodeRotateFun,type,weights)
  
  if(is.infinite(MaxDepth)) {
    numNode=min(numNode,sum(2^(0:ceiling(log2(n/MinLeaf))))) 
  }else{
    MaxDepth=min(MaxDepth,ceiling(log2(n/MinLeaf)),n/MinLeaf-1)
    numNode=min(numNode,sum(2^(0:MaxDepth)))
  }
  
  nodeXIndx = vector("list", numNode+1);
  nodeXIndx[[1]] = 1:n;
  
  if(type!="regression"){
    nodeNumLabel=matrix(0,0,maxLabel)
    colnames(nodeNumLabel)=Levels
    sl=seq(maxLabel)
  }else{
    nodeNumLabel=matrix(0,0,2)
    colnames(nodeNumLabel)=c("prediction","number")
  }

  nodeRotaMat=NULL;
  nodeDepth = rep(1,numNode);
  nodeCutValue = rep(0,numNode);
  nodeCutIndex = nodeCutValue;
  #nodeLabel = nodeCutValue;
  #nodeNumLabel = nodeCutValue;
  childNode = nodeCutValue;
  
  
  #start create pptree nodeFlags
  ##############################################################################
  currentNode = 1;
  freeNode=2
  while(!is.null(nodeXIndx[[currentNode]])){
    
    if((length(unique(y[nodeXIndx[[currentNode]]]))==1)||
       (length(nodeXIndx[[currentNode]])<=(MinLeaf+1))||
       (nodeDepth[currentNode]>=MaxDepth)||
       (freeNode>=numNode)){
      if(type!="regression"){
        leafLabel = table(Levels[c(sl,y[nodeXIndx[[currentNode]]])])-1
        #nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        #nodeNumLabel[currentNode]=max(leafLabel)
        #nodeNumLabel[currentNode,]=leafLabel
      }else{
        #nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
        #nodeNumLabel[currentNode]=length(y[nodeXIndx[[currentNode]]])
        leafLabel=c(mean(y[nodeXIndx[[currentNode]]]),length(y[nodeXIndx[[currentNode]]]))
      }
      nodeNumLabel=rbind(nodeNumLabel,leafLabel)
      
      nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
      nodeXIndx[currentNode]=NA
      currentNode = currentNode+1;
      next
    }
    
    if(!is.null(weights)){
      Wcd = weights[nodeXIndx[[currentNode]]];
    }else{
      Wcd = 1;
    }
    
    ##########################################
    if(NodeRotateFun=="RotMatMake"){
      sparseM <- RotMatMake(X[nodeXIndx[[currentNode]],], y[nodeXIndx[[currentNode]]], 
                            paramList$RotMatFun, paramList$PPFun, FunDir, paramList)
    }
    
    if(!NodeRotateFun%in%ls("package:ODRF")){
      paramList$x =X[nodeXIndx[[currentNode]],];
      paramList$y = y[nodeXIndx[[currentNode]]];
      sparseM <- do.call(FUN, paramList)
      paramList$x = NULL;paramList$y = NULL;
    }
    
    if(NodeRotateFun%in%c('RotMatRF','RotMatRand')){
      sparseM <- do.call(FUN, paramList)
    }
    
    if(NodeRotateFun=="RotMatPPO"){
      sparseM=RotMatPPO(x=X[nodeXIndx[[currentNode]],],y=y[nodeXIndx[[currentNode]]],numProj=paramList$numProj,
                        dimProj=paramList$dimProj,catLabel = paramList$catLabel,weights=paramList$weights,
                        model = paramList$model,type=paramList$type)#"NNet"
    }
    
    if(NodeRotateFun=="PPO"){
      sparseM=ODRF:::PPO(x=X[nodeXIndx[[currentNode]],],y=y[nodeXIndx[[currentNode]]],q = paramList$q,
                         ppMethod = paramList$ppMethod,weight = paramList$weight,r = paramList$r,
                         lambda = paramList$lambda,energy = paramList$energy,cooling = paramList$cooling,
                         TOL = paramList$TOL,maxiter = paramList$maxiter)
      sparseM=sparseM$projMat
    }
    
    numDr=unique(sparseM[,2]);
    rotaX=matrix(0,p,length(numDr));
    for(i in 1:length(numDr)){
      lrows = which(sparseM[,2] == numDr[i]);
      rotaX[sparseM[lrows, 1],i] = sparseM[lrows, 3];
    }
    ###################################################################
    
    rotaX=X[nodeXIndx[[currentNode]],,drop = FALSE] %*% rotaX;
    
    bestCut=ODRF:::best_cut_node(type0,rotaX,y[nodeXIndx[[currentNode]]],Wcd,MinLeaf,maxLabel);
    
    if(bestCut$BestCutVar==-1){
      TF=TRUE
    }else{
      Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
      TF=min(length(Lindex),length(nodeXIndx[[currentNode]])-length(Lindex))<=MinLeaf;
    }
    if(TF){
      if(type!="regression"){
        leafLabel = table(Levels[c(sl,y[nodeXIndx[[currentNode]]])])-1
        #nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        #nodeNumLabel[currentNode]=max(leafLabel)
        #nodeNumLabel[currentNode,]=leafLabel
      }else{
        #nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
        #nodeNumLabel[currentNode]=length(y[nodeXIndx[[currentNode]]])
        leafLabel=c(mean(y[nodeXIndx[[currentNode]]]),length(y[nodeXIndx[[currentNode]]]))
      }
      nodeNumLabel=rbind(nodeNumLabel,leafLabel)
      
      nodeRotaMat=rbind(nodeRotaMat,c(0,currentNode,0))
      nodeXIndx[currentNode]=NA
      currentNode = currentNode+1;
      next
    }
    
    sparseM=sparseM[sparseM[,2]==numDr[bestCut$BestCutVar],,drop = FALSE]
    sparseM[,2]=currentNode
    nodeRotaMat=rbind(nodeRotaMat,sparseM)
    rm(rotaX)
    rm(sparseM)
    
    nodeCutValue[currentNode] = bestCut$BestCutVal;
    nodeCutIndex[currentNode]=bestCut$BestIndex[bestCut$BestCutVar]
    childNode[currentNode]=freeNode;
    
    #Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
    nodeXIndx[[freeNode]] = nodeXIndx[[currentNode]][Lindex];
    nodeXIndx[[freeNode+1]] = nodeXIndx[[currentNode]][setdiff(1:length(nodeXIndx[[currentNode]]),Lindex)];
    nodeDepth[freeNode+c(0,1)]=nodeDepth[currentNode]+1
    
    nodeNumLabel=rbind(nodeNumLabel,0)
    
    nodeXIndx[currentNode]=NA
    freeNode = freeNode+2;
    currentNode = currentNode+1;
  }
  
  nodeDepth=nodeDepth[1:(currentNode-1)]
  colnames(nodeRotaMat)=c("var","node","coef")
  rownames(nodeRotaMat)=rep(nodeDepth,table(nodeRotaMat[,2]))
  rownames(nodeNumLabel)=nodeDepth
  
  ppTree=list(call=Call,terms=Terms,type=type,Levels=Levels,NodeRotateFun=NodeRotateFun,paramList=paramList)
  ppTree$data=list(subset=subset,weights=weights,na.action=na.action,n=n,p=p,varName=varName,
                   Xscale=Xscale,minCol=minCol,maxminCol=maxminCol,Xcat=Xcat,catLabel=catLabel,
                   TreeRandRotate=TreeRandRotate,rotdims=rotdims,rotmat=rotmat);
  ppTree$tree=list(FunDir=FunDir,MaxDepth=MaxDepth,MinLeaf=MinLeaf,numNode=numNode);
  ppTree$structure=list(nodeRotaMat = nodeRotaMat,nodeNumLabel = nodeNumLabel,nodeCutValue =nodeCutValue[1:(currentNode-1)],
                        nodeCutIndex =nodeCutIndex[1:(currentNode-1)],childNode = childNode[1:(currentNode-1)],
                        nodeDepth=nodeDepth)
  class(ppTree) <- "ODT"
  return(ppTree)
}
