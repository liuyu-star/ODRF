#' predict method for ODT objects
#'
#' Prediction a oblique decision tree based on an input matrix or data frame using \code{\link{ODT}} function.
#' 
#' @param ppTree an object of class ODT, as that created by the function ODT.
#' @param Xnew an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @param leafnode If or not output the leaf node sequence number that \code{Xnew} is partitioned. (default FALSE)
#' 
#' @return A set of vectors in the following list:
#' \itemize{
#' \item prediction: the prediced response of the new data.
#' \item leafnode: the leaf node sequence number that the new data is partitioned.
#' }
#' 
#' @references \itemize{
#' \item{Zhan H, Liu Y, Xia Y. Consistency of The Oblique Decision Tree and Its Random Forest[J]. arXiv preprint arXiv:2211.12653, 2022.}
#' }
#' 
#' @seealso \code{\link{ODT}}
#' 
#' @examples
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
#' @import Rcpp
#' @aliases predict.ODT
#' @rdname predict.ODT
#' @method predict ODT
#' @export
predict.ODT <- function(ppTree,Xnew,leafnode=FALSE){
  #ppTreeVar=c("type",names(ppTree$structure),names(ppTree$data),names(ppTree$tree))
  #ppTree=do.call("c",ppTree)
  #assign("Levels", as.vector(unlist(ppTree[c(2,3)])))
  #ppTree=ppTree[-c(2,3)]
  #for(v in seq(length(ppTreeVar))){
  #  assign(ppTreeVar[v], ppTree[[v]])
  #}
  #rm(ppTree)
  # address na values.
  if (any(is.na(Xnew))) {
    Xnew=ppTree$data$na.action(data.frame(Xnew))
    warning("NA values exist in data matrix")
  }
  Xnew=as.matrix(Xnew) 
  
  if(!is.null(ppTree$data$subset))
    Xnew=Xnew[ppTree$data$subset,]
  #weights0=c(ppTree$data$weights,ppTree$paramList$weights)
  #if(!is.null(ppTree$data$weights))
  #  Xnew <- Xnew * matrix(weights,length(y),ncol(Xnew))
  
  p=ncol(Xnew)
  n=nrow(Xnew)
  
  if(ppTree$type!="regression"){
    nodeLabel=colnames(ppTree$structure$nodeNumLabel)[max.col(ppTree$structure$nodeNumLabel)] ## "random"
    nodeLabel[which(rowSums(ppTree$structure$nodeNumLabel)==0)]=0
  }else{
    nodeLabel=as.character(ppTree$structure$nodeNumLabel[,1])
  }
  
  if(all(ppTree$structure$nodeCutValue==0)){
    pred=rep(nodeLabel,n)
  }else{
    Xcat=ppTree$data$Xcat
    catLabel=ppTree$data$catLabel
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
    if(ppTree$data$Xscale!="No"){
      indp=(numCat+1):p
      Xnew[,indp]=(Xnew[,indp]-matrix(ppTree$data$minCol,n,length(indp),byrow = T))/
        matrix(ppTree$data$maxminCol,n,length(indp),byrow = T)
    }
    
    if (ppTree$data$TreeRandRotate) {
      Xnew[,ppTree$data$rotdims] <- Xnew[,ppTree$data$rotdims,drop = FALSE] %*% ppTree$data$rotmat
    }
    
    pred = .Call('_ODRF_predict_ODT', PACKAGE = 'ODRF',Xnew, ppTree$structure$nodeRotaMat, 
                       ppTree$structure$nodeCutValue, ppTree$structure$childNode,nodeLabel)
    
    if(leafnode){
      pred=pred$node
    }else{
      pred=pred$prediction
    }
  }
  
  if(ppTree$type=='regression'){pred=as.numeric(pred)}
  
  return(pred)
}
