#' carPPtree Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#' @param Xnew an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction
#'
#' @import Rcpp
#' @import grid
#' 
#' @aliases predict.ODT
#' @rdname predict.ODT
#' @method predict ODT
#' @export
#' 
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#'
#' ### Train RerF on one-of-K encoded categorical data ###
#' df1 <- as.data.frame(Titanic)
#' nc <- ncol(df1)
#' df2 <- df1[NULL, -nc]
#' for (i in which(df1$Freq != 0L)) {
#'   df2 <- rbind(df2, df1[rep(i, df1$Freq[i]), -nc])
#' }
#' 
## @rdname predict#'  #' @aliases predict_ODT #' @method predict ODT #' @useDynLib ppRF
#@rdname predict
#@method predict ODT
predict.ODT <- function(ppTree,Xnew,type=c("prediction","node")[1]){
  #ppTreeVar=c("method",names(ppTree$structure),names(ppTree$data),names(ppTree$tree))
  #ppTree=do.call("c",ppTree)
  #assign("Levels", as.vector(unlist(ppTree[c(2,3)])))
  #ppTree=ppTree[-c(2,3)]
  #for(v in seq(length(ppTreeVar))){
  #  assign(ppTreeVar[v], ppTree[[v]])
  #}
  #rm(ppTree)
  Xnew=as.matrix(Xnew)
  p=ncol(Xnew)
  n=nrow(Xnew)
  
  if(ppTree$method!="regression"){
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
    
    if(type=="node"){
      pred=pred$node
    }else{
      pred=pred$prediction
    }
  }
  
  if((ppTree$method=='regression')&(type=="prediction")){pred=as.numeric(pred)}
  
  #class(pred) <- append(class(pred),"ppCART")
  
  return(pred)
}

#predict <- function(x, ...)  UseMethod("predict")

