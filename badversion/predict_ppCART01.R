#' carPPtree Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param Xnew an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction
#'
#' @useDynLib ppRF
#' @import Rcpp
#' 
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
## @rdname predict#'  #' @aliases predict_ppCART #' @method predict ppCART
predict_ppCART <- function(Xnew, ppTree){
  if(all(ppTree$nodeCutValue==0)){
    prediction=rep(ppTree$nodeLabel,nrow(Xnew))
  }else{
    p=ncol(Xnew)
    n=nrow(Xnew)
    
    Xcat=ppTree$Xcat
    catLabel=ppTree$catLabel
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
    if(ppTree$Xscale!="No"){
      indp=(numCat+1):p
      Xnew[,indp]=(Xnew[,indp]-matrix(ppTree$minCol,n,length(indp),byrow = T))/
        matrix(ppTree$maxminCol,n,length(indp),byrow = T)
    }
    
    if (ppTree$TreeRandRotate) {
      Xnew[, ppTree$rotdims] <- Xnew[,ppTree$rotdims,drop = FALSE] %*% ppTree$rotmat
    }
    
    prediction = .Call('_ppRF_ppCARTPredict', PACKAGE = 'ppRF',as.matrix(Xnew), ppTree$nodeRotaMat, 
                       ppTree$nodeCutValue, ppTree$childNode, as.character(ppTree$nodeLabel))
  }
  
  if(ppTree$method=='regression'){prediction=as.numeric(prediction)}
  
  return(prediction)
}