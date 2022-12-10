#' carPPtree Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param X an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
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
ppCARTPredict01 <- function(Xnew, PPtree){
  if(all(PPtree$nodeCutValue==0)){
    predict=rep(PPtree$nodeLabel,nrow(Xnew))
  }else{
    p=ncol(Xnew)
    n=nrow(Xnew)
    
    Xcat=PPtree$Xcat
    numCat=0
    if(sum(Xcat)>0){
      numCat <- apply(Xnew[,Xcat,drop = FALSE], 2, function(x) length(unique(x)))
      Xnew1 <- matrix(0, nrow = n, ncol = sum(numCat)) # initialize training data matrix X
      catMap <- vector("list", length(Xcat))
      names(catMap)<- colnames(Xnew)[Xcat]
      col.idx <- 0L
      # one-of-K encode each categorical feature and store in X
      for (j in 1:length(Xcat)) {
        catMap[[j]] <- (col.idx + 1L):(col.idx + numCat[j])
        # convert categorical feature to K dummy variables
        indC=unique(Xnew[,Xcat[j]])
        Xnew1[, catMap[[j]]] <- (matrix(Xnew[,Xcat[j]],n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
        col.idx <- col.idx + numCat[j]
      }
      Xnew=cbind(Xnew1,Xnew[,-Xcat]) 
      p=ncol(Xnew)
      rm(Xnew1)
    }
    
    #Variable scaling.
    if(PPtree$Xscale!="No"){
      indp=(sum(numCat)+1):p
      Xnew[,indp]=(Xnew[,indp]-matrix(PPtree$minCol,n,length(indp),byrow = T))/
        matrix(PPtree$maxminCol,n,length(indp),byrow = T)
    }
    
    if (PPtree$rotate) {
      Xnew[, rotdims] <- Xnew[,PPtree$rotdims,drop = FALSE] %*% PPtree$rotmat
    }
    
    predict = .Call('_ppRF_ppCARTPredict', PACKAGE = 'ppRF',as.matrix(Xnew), PPtree$nodeRotaMat, 
                    PPtree$nodeCutValue, PPtree$childNode, as.character(PPtree$nodeLabel))
  }
  
  if(PPtree$method=='r'){predict=as.numeric(predict)}
  
  return(predict)
}