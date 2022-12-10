#' find  best split variable and node.
#'
#' @param       X       : an n by d numeric matrix (preferable) or data frame.
#' 
#' @param       y       : a n vector.
#' 
#' @param       method       : the criterion used for splitting the nodes
#'                           'g-classification' : gini impurity index (classification, default)
#'                           'i-classification' : information gain (classification)
#'                           'regression' : mean square error (regression)
#'
#' @param       MinLeaf      : the minimum amount of samples in a leaf
#'
#' @param       weights      : a vector of values which weigh the samples when considering a split
#'
#' @param       numLabels    : the number of categories
#' 
#' 
#' @return a list which contains:
#' \itemize{
#' \item BestCutVar: the best split variable
#' \item BestCutVal: the best split point for the best split variable.
#' \item BestIndex: Each variable corresponds to the min gini impurity index(method='g-classification'),
#' the max information gain(method='i-classification') and the min squared error(method='regression')
#' }
#'
#' @import Rcpp
#' 
#' @export
#' 
#' @examples
#' ### Find the best split variable ###
#' library(ppRF)
#' X=as.matrix(iris[, 1:4])
#' y=iris[[5L]]
#'
#' bestcut=BestCutNode(X,y,method='g-classification')
#' print(bestcut)
best_split_node <- function(X, y, method='g-classification', weights=1, MinLeaf=ifelse(method=='regression',5,1),
                          numLabels=ifelse(method=='regression',0,length(unique(y)))) {

  X <- as.matrix(X)
  if(method!="regression"){
    y <- as.integer(as.factor(y));
  }else{
    y=c(y)
  }
  
  
  .Call(`_ODRF_best_cut_node`, strsplit(method,split = "")[[1]][1], X, y, weights, MinLeaf, numLabels)
}

