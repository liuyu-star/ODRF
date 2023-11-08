#' Extract variable importance measure
#'
#' This is the extractor function for variable importance measures as produced by \code{\link{ODT}} and \code{\link{ODRF}}.
#'
#' @param obj An object of class \code{\link{ODT}} and \code{\link{ODRF}}.
#' @param type specifying the type of importance measure. "impurity": mean decrease in node impurity, "permutation" (default): mean decrease in accuracy.
#' @param X An n by d numerical matrix (preferably) or data frame is used in the \code{ODRF}.
#' @param y A response vector of length n is used in the \code{ODRF}.
#'
#' @return A matrix of importance measure, first column is the predictors and second column is Increased error. Misclassification rate (MR) for classification or mean square error (MSE) for regression.
#' The larger the increased error the more important the variable is.
#'
#' @details A note from \code{randomForest} package, here are the definitions of the variable importance measures.
#' \itemize{
#' \item The first measure is the total decrease in node impurities from splitting on the variable, averaged over all trees. For classification, the node impurity is measured by the Gini index. For regression, it is measured by residual sum of squares.
#' \item The second measure is computed from permuting OOB data: For each tree, the prediction error on the out-of-bag portion of the data is recorded.
#' Then the same is done after permuting each predictor variable. The difference between the two are then averaged over all trees.
#' }
#'
#' @seealso \code{\link{ODRF}} \code{\link{Accuracy}} \code{\link{plot.VarImp}}
#'
#' @examples
#' data(body_fat)
#' y=body_fat[,1]
#' X=body_fat[,-1]
#'
#' tree <- ODT(X, y, split = "mse")
#' (varimp <- VarImp(tree, type="impurity"))
#'
#' forest <- ODRF(X, y, split = "mse", parallel = FALSE, ntrees=50)
#' (varimp <- VarImp(forest, type="impurity"))
#' (varimp <- VarImp(forest, X, y, type="permutation"))
#'
#' @keywords forest tree
#' @export
VarImp <- function(obj, X=NULL, y=NULL, type="permutation") {

  p=obj$data$p
  varName <- obj$data$varName
  pp <- obj$data$p
  Xcat <- obj$data$Xcat
  catLabel <- obj$data$catLabel
  Rep=rep(seq(p),rep(1,pp))
  sp=seq(pp)
  if (!is.null(catLabel) && (sum(Xcat) > 0)) {
    pp <- pp - length(unlist(catLabel)) + length(Xcat)
    varName <-c(names(catLabel),varName[-seq(length(unlist(catLabel)))])
    Rep=rep(seq(pp),c(unlist(lapply(catLabel, length)),rep(1,pp-length(Xcat))))
    sp=seq(pp)
    sp=c(Xcat,sp[-Xcat])
  }

  if(type=="impurity"){
    VarImp.impurity <- function(structure){
      DecImpurity=rep(0,p)

      nodeRotaMat=structure[["nodeRotaMat"]]
      if(length(structure[["nodeCutIndex"]])==1){
        stop("No tree structure to measure the importance of variables!")
      }else{
        cutNodes=unique(nodeRotaMat[nodeRotaMat[,1]!=0,2])
      }

      for (node in  cutNodes) {
        idx=which(nodeRotaMat[,2]==node)
        DecImpurity[nodeRotaMat[idx,1]]=DecImpurity[nodeRotaMat[idx,1]]+
          (nodeRotaMat[idx,3])^2*structure[["nodeCutIndex"]][node]
      }

      DecImpurity <- DecImpurity/length(cutNodes)

      return(DecImpurity)
    }

    if("ODT"%in%class(obj)){
      DecImpurity=VarImp.impurity(obj$structure)
    }

    if("ODRF"%in%class(obj)){
      DecImpurity <- vapply(obj$structure, function(structure){
        if(length(structure[["nodeCutIndex"]])==1){
          rep(0,p)
        }else{
          VarImp.impurity(structure)
        }
        } , rep(0,p))

      DecImpurity <- rowMeans(DecImpurity)
    }

    DecImpurity <- aggregate(DecImpurity, list(Rep), mean)[-1]
    varimport <- cbind(sp,DecImpurity)
    colnames(varimport) <- c("varible","decrease_accuracy")
    rownames(varimport) <- varName
  }


  #####################################################################
  if(type=="permutation"){
    X <- as.matrix(X)
    colnames(X) <- varName
    if (ncol(X) != pp) {
      stop("The dimensions of 'Xnew' and training data do not match")
    }

    if("ODT"%in%class(obj)){
      stop("Tree structure can't use permutation method to measure the importance of variables!")
    }

    if("ODRF"%in%class(obj)){
      if(is.null(X)&&is.null(y)){
        stop("training data 'X' and 'y' argument is required")
      }

      if (!obj$forest$storeOOB) {
        stop("out-of-bag indices for each tree are not stored, so can't use permutation method!")
      }

      if (!is.null(obj$data$subset)) {
        X <- X[obj$data$subset, ]
      }
      # if(!is.null(weights))
      #  X <- X * matrix(weights0,length(y),ncol(X))
      # weights=weights0


      if (obj$split != "mse") {
        #y <- factor(y, levels = obj$Levels)
        y <- as.character(y)
      }
      # X=obj$data$X
      # y=obj$data$y
      # ntrees=obj$tree$ntrees
      n <- length(y)
      p <- ncol(X)

      numCat <- 0
      if (sum(Xcat) > 0) {
        xj <- 1
        X1 <- matrix(0, nrow = n, ncol = length(unlist(catLabel))) # initialize training data matrix X
        # one-of-K encode each categorical feature and store in X
        for (j in seq_along(Xcat)) {
          catMap <- which(catLabel[[j]] %in% unique(X[, Xcat[j]]))
          indC <- catLabel[[j]][catMap]
          Xj <- (matrix(X[, Xcat[j]], n, length(indC)) == matrix(indC, n, length(indC), byrow = TRUE)) + 0

          if (length(indC) > length(catLabel[[j]])) {
            Xj <- Xj[, seq_along(catLabel[[j]])]
          }

          xj1 <- xj + length(catLabel[[j]])
          X1[, (xj:(xj1 - 1))[catMap]] <- Xj
          xj <- xj1
        }

        X <- cbind(X1, X[, -Xcat])
        p <- ncol(X)
        numCat <- length(unlist(catLabel))
        rm(X1)
        rm(Xj)
      }
      if (!is.numeric(X)){
        X=apply(X, 2, as.numeric)
      }

      # Variable scaling.
      if (obj$data$Xscale != "No") {
        indp <- (sum(numCat) + 1):p
        X[, indp] <- (X[, indp] - matrix(obj$data$minCol, n, length(indp), byrow = T)) /
          matrix(obj$data$maxminCol, n, length(indp), byrow = T)
      }

      split=obj$split
      Levels=obj$Levels
      runOOBErr <- function(structure, ...) {
        oobErrs <- rep(0, p)
        oobIndex <- structure$oobIndex
        X0 <- X[oobIndex, ]
        y0 <- y[oobIndex]
        n0 <- length(y0)
        # if(obj$split=="regression"){
        #  e.0=mean((yi-mean(y[-oobIndex]))^2)
        # }

        pred <- predictTree(structure,X0,split,Levels)$prediction
        if (split != "mse") {
          oobErr0 <- mean(pred != y0)
        } else {
          oobErr0 <- mean((pred - y0)^2) # /e.0
        }

        for (j in seq(p)) {
          Xi=X0
          Xi[, j] <- Xi[sample(n0), j] #+rnorm(length(oobIndex))
          pred <- predictTree(structure,Xi,split,Levels)$prediction
          if (split != "mse") {
            oobErr <- mean(pred != y0)
          } else {
            oobErr <- mean((pred - y0)^2) # /e.0
          }

          oobErrs[j] <- oobErr - oobErr0#abs()
        }

        return(oobErrs)
      }

      IncErr <- vapply(obj$structure, runOOBErr, rep(0,p))
      #IncErr <- sapply(obj$ppTrees, runOOBErr)
      IncErr <- rowMeans(IncErr)
      IncErr <- aggregate(IncErr, list(Rep), mean)[-1]
      varimport <- cbind(sp,IncErr)
      colnames(varimport) <- c("varible","decrease_accuracy")
      rownames(varimport) <- varName
    }
  }

  varimport <- list(varImp = varimport[order(varimport[,2], decreasing = T), ], split = obj$split)
  class(varimport) <- "VarImp"

  return(varimport)
}
