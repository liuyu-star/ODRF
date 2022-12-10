#' carPPtree Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param X an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction
#'
#' @useDynLib ppRF1
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
  PPtreePredict <- function(Xnew, PPtree){
    if(all(PPtree$nodeCutVar==0)){
        predict=rep(PPtree$nodelabel,nrow(Xnew))
    }else{
      predict = .Call('_ppRF1_PPtreePredict', PACKAGE = 'ppRF1',as.matrix(Xnew), PPtree$nodesparseM, PPtree$nodeCutVar, 
                     PPtree$nodeCutValue, PPtree$childnode, as.character(PPtree$nodelabel))
    }
    if(PPtree$method=='r'){predict=as.numeric(predict)}
    return(predict)
}