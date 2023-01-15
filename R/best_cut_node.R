#' find best split variable and node
#'
#' A function to select the splitting variables and nodes using one of three criteria.
#'
#' @param X An n by d numeric matrix (preferable) or data frame.
#' @param y A response vector of length n.
#' @param type One of three criteria, 'i-classification': information gain (classification, default),
#' 'g-classification': gini impurity index (classification) or 'regression': mean square error (regression).
#' @param MinLeaf The minimum amount of samples in a leaf.
#' @param weights A vector of values which weigh the samples when considering a split.
#' @param numLabels The number of categories.
#'
#' @return A list which contains:
#' \itemize{
#' \item BestCutVar: The best split variable.
#' \item BestCutVal: The best split point for the best split variable.
#' \item BestIndex: Each variable corresponds to the min gini impurity index(method='g-classification'),
#' the max information gain(method='i-classification') or the min squared error(method='regression').
#' }
#'
#' @examples
#' ### Find the best split variable ###
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' y <- iris[[5]]
#' bestcut <- best.cut.node(X, y, type = "i-classification")
#' print(bestcut)
#'
#' @export
best.cut.node <- function(X, y, type = "i-classification", weights = 1, MinLeaf = ifelse(type == "regression", 5, 1),
                          numLabels = ifelse(type == "regression", 0, length(unique(y)))) {
  if (any(is.na(X))) {
    stop("data 'X' has Missing value, NA or NaN")
  }

  X <- as.matrix(X)
  if (type != "regression") {
    y <- as.integer(as.factor(y))
  } else {
    y <- c(y)
  }

  .Call("_ODRF_best_cut_node", PACKAGE = "ODRF", strsplit(type, split = "")[[1]][1], X, y, weights, MinLeaf, numLabels)
}
