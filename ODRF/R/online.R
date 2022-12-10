#' RerF Tree Generator
#' 
#' Grows a CARTree using X(samplesXfeatures) and one of the following criteria which can be set via the parameter 'method'
#'
#' @param       'g-classification' : gini impurity index (classification)
#' @param       'i-classification' : information gain (classification, default)
#' @param       'regression' : squared error (regression)
#' 
#' Other parameters that can be set are:
#' 
#' @param paramList parameters in a named list to be used by FUN. If left unchanged,
#' @param catMap a list specifying which columns in X correspond to the same one-of-K encoded feature. Each element of catMap is a numeric vector specifying the K column indices of X corresponding to the same categorical feature after one-of-K encoding. All one-of-K encoded features in X must come after the numeric features. The K encoded columns corresponding to the same categorical feature must be placed contiguously within X. The reason for specifying catMap is to adjust for the fact that one-of-K encoding cateogorical features results in a dilution of numeric features, since a single categorical feature is expanded to K binary features. If catMap = NULL, then RerF assumes all features are numeric (i.e. none of the features have been one-of-K encoded).
#' @param FUN :FUN=RandMatPPR.
#' @param minparent    : the minimum amount of samples in an impure node for it to be considered for splitting (default 2)
#' @param MinLeaf      : the minimum amount of samples in a leaf (default 1)
#' @param Xweights      : a vector of values which weigh the samples when considering a split (default [])
#' @param nvartosample : the number of (randomly selected) variables to consider at each node (default all)
#' @return ppTree
#'
#' @useDynLib ODRF
#' @import Rcpp
#' 
#' @export
#'
#'
#' @examples
#' ### Train RerF on numeric data ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#' 
#X=X;y=y;method='g-classification';NodeRotateFun="RotMatPPO";paramList=NULL;catLabel=NULL;MaxDepth=Inf;
#MinLeaf=ifelse(method=='regression',5,1);Xweights=1;numNode=Inf;
#Levels=if(method=='regression')NULL else levels(as.factor(y));Xcat=NULL;
#Xscale=c("Min-max","Quantile","No")[1];FunDir=getwd();TreeRandRotate=FALSE
#Xcat=c(NULL,0)[1]
online <- function(obj, ...)
  UseMethod("online",obj)
