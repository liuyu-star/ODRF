#' find best splitting variable and node
#'
#' A function to select the splitting variables and nodes using one of four criteria.
#'
#' @param X An n by d numeric matrix (preferable) or data frame.
#' @param y A response vector of length n.
#' @param Xsplit Splitting variables used to construct linear model trees. The default value is NULL and is only valid when split="linear".
#' @param split The criterion used for splitting the nodes. "entropy": information gain and "gini": gini impurity index for classification; "": mean square error for regression; "linear": mean square error for multiple linear regression.
#' @param lambda The argument of \code{split} is used to determine the penalty level of the partition criterion. Three options are provided including, \code{lambda=0}: no penalty; \code{lambda=2}: AIC penalty; \code{lambda='log'} (Default): BIC penalty. In Addition, lambda can be any value from 0 to n (training set size).
#' @param MinLeaf Minimal node size (Default 10).
#' @param weights A vector of values which weigh the samples when considering a split.
#' @param numLabels The number of categories.
#' @param glmnetParList List of parameters used by the functions \code{glmnet} and \code{cv.glmnet} in package \code{glmnet}. If left unchanged, default values will be used, for details see \code{\link[glmnet]{glmnet}} and \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return A list which contains:
#' \itemize{
#' \item BestCutVar: The best split variable.
#' \item BestCutVal: The best split points for the best split variable.
#' \item BestIndex: Each variable corresponds to maximum decrease in gini impurity index, information gain, and mean square error.
#' \item fitL and fitR: The multivariate linear models for the left and right nodes after splitting are trained using the function \code{\link{glmnet}}.
#' }
#'
#' @examples
#' ### Find the best split variable ###
#' #Classification
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' y <- iris[[5]]
#' (bestcut <- best.cut.node(X, y, split = "gini"))
#' (bestcut <- best.cut.node(X, y, split = "entropy"))
#'
#' #Regression
#' data(body_fat)
#' X=body_fat[, -1]
#' y=body_fat[, 1]
#' (bestcut <- best.cut.node(X, y, split = "mse"))
#'
#' set.seed(10)
#' cutpoint=50
#' X=matrix(rnorm(100*10),100,10)
#' age=sample(seq(20,80),100,replace = TRUE)
#' height=sample(seq(50,200),100,replace = TRUE)
#' weight=sample(seq(5,150),100,replace = TRUE)
#' Xsplit=cbind(age=age,height=height,weight=weight)
#' mu=rep(0,100)
#' mu[age<=cutpoint]=X[age<=cutpoint,1]+X[age<=cutpoint,2]
#' mu[age>cutpoint]=X[age>cutpoint,1]+X[age>cutpoint,3]
#' y=mu+rnorm(100)
#' bestcut <- best.cut.node(X, y, Xsplit, split = "linear",
#'            glmnetParList=list(lambda = 0))
#'
#' @export
best.cut.node <- function(X, y, Xsplit=X, split, lambda = "log", weights = 1, MinLeaf = 10,
                          numLabels = ifelse(split %in% c("gini","entropy"), length(unique(y)), 0),glmnetParList=NULL) {
  if (any(is.na(X))) {
    stop("data 'X' has Missing value, NA or NaN")
  }

  X <- as.matrix(X)
  Xsplit <- as.matrix(Xsplit)
  if (split %in% c("gini","entropy")) {
    y <- as.integer(as.factor(y))
  } else {
    y <- c(y)
  }

  if (lambda == "log") {
    lambda <- length(y)
  }

  if(split == "linear"){
    glmnetParList$weights=if(length(weights)==1)NULL
    bestcut=linear_split(X, Xsplit, y, MinLeaf, glmnetParList)
  }else{
    #if (split == "") method <- "r"
    #if (split == "entropy") method <- "i"
    #if (split == "gini") method <- "g"
    method <- switch(split,
                     "mse" = "r",
                     "entropy" = "i",
                     "gini" = "g")
    # strsplit(split, split = "")[[1]][1]
    bestcut=.Call("_ODRF_best_cut_node", PACKAGE = "ODRF", method, lambda, Xsplit, y, weights, MinLeaf, numLabels)
  }

  return(bestcut)
}
