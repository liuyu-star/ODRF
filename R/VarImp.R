#' variable importance of oblique decision random forest
#'
#' Variable importance is computed from permuting OOB data.
#'
#' @param obj An object of class \code{\link{ODRF}}.
#' @param X An n by d numerical matrix (preferably) or data frame is used in the \code{ODRF}.
#' @param y A response vector of length n is used in the \code{ODRF}.
#'
#' @return A matrix of importance measure, first column is the predictors and second column is Increased error. Misclassification rate (MR) for classification or mean square error (MSE) for regression.
#' The larger the increased error the more important the variable is.
#'
#' @details A note from \code{randomForest} package, here are the definitions of the variable importance measures. The measure is computed from permuting OOB data: For each tree, the prediction error on the out-of-bag portion of the data is recorded.
#' Then the same is done after permuting each predictor variable. The difference between the two are then averaged over all trees.
#'
#' @seealso \code{\link{ODRF}} \code{\link{Accuracy}} \code{\link{plot.VarImp}}
#'
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train <- sample(1:569, 200)
#' train_data <- data.frame(breast_cancer[train, -1])
#' forest <- ODRF(train_data[, -1], train_data[, 1], split = "gini",
#'   parallel = FALSE)
#' varimp <- VarImp(forest, train_data[, -1], train_data[, 1])
#' varimp
#' @keywords forest
#' @export
VarImp <- function(obj, X, y) {
  # vars=all.vars(forest$terms)
  # address na values.
  # if (any(is.na(data))) {
  #  data=forest$data$na.action(data.frame(data))
  #  warning("NA values exist in data matrix")
  # }
  # y= data[,setdiff(colnames(data),vars[-1])]
  # X= data[,vars[-1]]

  forest=obj
  X <- as.matrix(X)
  colnames(X)=forest$data$varName

  pp <- forest$data$p
  if (!is.null(forest$data$catLabel) && (sum(forest$data$Xcat) > 0)) {
    pp <- pp - length(unlist(forest$data$catLabel)) + length(forest$data$Xcat)
  }
  if (ncol(X) != pp) {
    stop("The dimensions of 'Xnew' and training data do not match")
  }


  if (!is.null(forest$data$subset)) {
    X <- X[forest$data$subset, ]
  }
  # if(!is.null(weights))
  #  X <- X * matrix(weights0,length(y),ncol(X))
  # weights=weights0


  if (forest$split != "mse") {
    y <- factor(y, levels = forest$Levels)
  }
  # X=forest$data$X
  # y=forest$data$y
  # ntrees=forest$tree$ntrees
  n <- length(y)
  p <- ncol(X)

  Xcat <- forest$data$Xcat
  catLabel <- forest$data$catLabel
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
  if (forest$data$Xscale != "No") {
    indp <- (sum(numCat) + 1):p
    X[, indp] <- (X[, indp] - matrix(forest$data$minCol, n, length(indp), byrow = T)) /
      matrix(forest$data$maxminCol, n, length(indp), byrow = T)
  }

  runOOBErr <- function(tree, ...) {
    class(tree) <- "ODT"
    oobErrs <- rep(0, p)
    oobIndex <- tree$oobIndex
    X0 <- X[oobIndex, ]
    y0 <- y[oobIndex]
    n0 <- length(y0)
    # if(forest$split=="regression"){
    #  e.0=mean((yi-mean(y[-oobIndex]))^2)
    # }

    pred <- predict(tree, X0)
    if (forest$split != "mse") {
      oobErr0 <- mean(pred != y0)
    } else {
      oobErr0 <- mean((pred - y0)^2) # /e.0
    }

    for (j in seq(p)) {
      Xi=X0
      Xi[, j] <- Xi[sample(n0), j] #+rnorm(length(oobIndex))
      pred <- predict(tree, Xi)
      if (forest$split != "mse") {
        oobErr <- mean(pred != y0)
      } else {
        oobErr <- mean((pred - y0)^2) # /e.0
      }

      oobErrs[j] <- abs(oobErr - oobErr0)
    }

    return(oobErrs)
  }

  IncErr <- vapply(forest$ppTrees, runOOBErr, rep(0,p))
  #IncErr <- sapply(forest$ppTrees, runOOBErr)
  IncErr <- rowMeans(IncErr)
  varimport <- cbind(varible = seq(p), increased_error = IncErr)
  rownames(varimport) <- colnames(X)

  varimport <- list(varImp = varimport[order(IncErr, decreasing = T), ], split = forest$split)
  class(varimport) <- "VarImp"

  return(varimport)
}
