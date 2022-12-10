#' Default values passed to RandMat*
#'
#' Given the parameter list and the categorical map this function
#' populates the values of the parameter list accoding to our 'best'
#' known general use case parameters.
#'
#' @param ncolX an integer denoting the number of columns in the design
#' matrix X.
#' @param paramList a list (possibly empty), to be populated with a set
#' of default values to be passed to a RandMat* function.
#' @param cat.map a list specifying which columns in X correspond to the
#' same one-of-K encoded feature. Each element of cat.map is a numeric
#' vector specifying the K column indices of X corresponding to the same
#' categorical feature after one-of-K encoding. All one-of-K encoded
#' features in X must come after the numeric features. The K encoded
#' columns corresponding to the same categorical feature must be placed
#' contiguously within X. The reason for specifying cat.map is to adjust
#' for the fact that one-of-K encoding cateogorical features results in
#' a dilution of numeric features, since a single categorical feature is
#' expanded to K binary features. If cat.map = NULL, then RerF assumes
#' all features are numeric (i.e. none of the features have been
#' one-of-K encoded).
#'
#' @return If \code{cat.map} is NULL, then
#' \itemize{
#' \item \code{p} is set to the number of columns of \code{X}
#' \item \code{d} is set to the ceiling of the square root of the number of columns of \code{X}
#' \item \code{sparsity}: if \eqn{\code{ncol(X)} \ge 10}, then sparsity is set
#' to 3 / \code{ncol{X}}, otherwise it is set to 1 / \code{ncol(X)}.
#' \item \code{prob} defaults to 0.5.
#' }
#'
#' @keywords internal
#'
defaults <- function(paramList,dimX,catLabel,NodeRotateFun,method,weights) {
  #public parameter.
  if (is.null(paramList[["dimProj"]])) {
    paramList$dimProj = NULL
  }
  if (is.null(paramList[["dimX"]])) {
    paramList$dimX <- dimX
  }
  if (is.null(paramList[["catLabel"]])) {
    paramList$catLabel = catLabel
  }
  if (is.null(paramList[["numProj"]])) {
    #  q<- min(ceiling(length(y)^0.4),ceiling(paramList$p*2/3))
    #  paramList$d <-min(max(5, ceiling(paramList$p/q)),paramList$p)
    paramList$numProj <-ifelse(NodeRotateFun=="RotMatPPO",NULL,ceiling(sqrt(paramList$dimX)))
  }
  if (is.null(paramList[["ppMethod"]])) {
    paramList$ppMethod =ifelse(NodeRotateFun=="PPO","LDA","PPR")
  }
  if (is.null(paramList[["weights"]])) {
    paramList$weights = weights
  }
  if (is.null(paramList[["lambda"]])) {
    paramList$lambda = ifelse(NodeRotateFun=="PPO",0.1,1)
  }
  
  #RandRotMat parameter.
  #############################################
  if (is.null(paramList[["sparsity"]])) {
    paramList$sparsity <- ifelse(paramList$dimX >=10, 3/paramList$dimX, 1/paramList$dimX)
  }
  if (is.null(paramList[["prob"]])) {
    paramList$prob <- 0.5
  }
  if (is.null(paramList[["RandDist"]])) {
    paramList$RandDist = c("Binary","Norm","Uniform")[1]
  }
  if (is.null(paramList[["RandDist"]])) {
    paramList$RandDist = "Binary"
  }
  
  #RotMatPPO parameter.
  ###########################################
  if (is.null(paramList[["method"]])) {
    paramList$method = method#method='g-classification'
  }
  
  #################################################
  #pp optimization parameter.
  if (is.null(paramList[["q"]])) {
    paramList$q=1
  }
  if (is.null(paramList[["weight"]])) {
    paramList$weight = TRUE
  }
  if (is.null(paramList[["r"]])) {
    paramList$r=1
  }
  if (is.null(paramList[["energy"]])) {
    paramList$energy=0
  }
  if (is.null(paramList[["cooling"]])) {
    paramList$cooling=0.9
  }
  if (is.null(paramList[["TOL"]])) {
    paramList$TOL=0.0001
  }
  if (is.null(paramList[["maxiter"]])) {
    paramList$maxiter=1000
  }

  
  return(paramList)
}
