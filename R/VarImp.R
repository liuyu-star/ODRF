#' variable importance of oblique decision random forest
#'
#' Variable importance is computed from permuting OOB data.
#'
#' @param forest An object of class \code{\link{ODRF}}.
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
#' test_data <- data.frame(breast_cancer[-train, -1])
#'
#' forest <- ODRF(diagnosis ~ ., train_data, type = "i-classification", parallel = FALSE)
#' (varimp <- VarImp(forest, train_data[, -1], train_data[, 1]))
#'
#' @keywords forest
#' @export
VarImp <- function(forest, X, y) {
  # vars=all.vars(forest$terms)
  # address na values.
  # if (any(is.na(data))) {
  #  data=forest$data$na.action(data.frame(data))
  #  warning("NA values exist in data matrix")
  # }
  # y= data[,setdiff(colnames(data),vars[-1])]
  # X= data[,vars[-1]]
  X <- as.matrix(X)

  if (!is.null(forest$data$subset)) {
    X <- X[forest$data$subset, ]
  }
  # if(!is.null(weights))
  #  X <- X * matrix(weights0,length(y),ncol(X))
  # weights=weights0


  if (forest$type != "regression") {
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

  # Variable scaling.
  if (forest$data$Xscale != "No") {
    indp <- (sum(numCat) + 1):p
    X[, indp] <- (X[, indp] - matrix(forest$data$minCol, n, length(indp), byrow = T)) /
      matrix(forest$data$maxminCol, n, length(indp), byrow = T)
  }

  runOOBErr <- function(tree, ...) {
    class(tree) <- "ODT"
    oobErrs <- rep(0, p + 1)
    oobIndex <- tree$oobIndex

    Xi <- X[oobIndex, ]
    yi <- y[oobIndex]
    yn <- length(yi)
    # if(forest$type=="regression"){
    #  e.0=mean((yi-mean(y[-oobIndex]))^2)
    # }
    for (j in 1:(p + 1)) {
      if (j != 1) {
        Xi[, j - 1] <- Xi[sample(yn), j - 1] #+rnorm(length(oobIndex))
      }
      pred <- predict(tree, Xi)
      if (forest$type != "regression") {
        oobErr <- mean(pred != yi)
      } else {
        oobErr <- mean((pred - yi)^2) # /e.0
      }

      if (j == 1) {
        oobErrs[1] <- oobErr
      } else {
        oobErrs[j] <- abs(oobErr - oobErrs[1])
      }
    }

    return(oobErrs)
  }

  oobErrVar <- sapply(forest$ppTrees, runOOBErr)[-1, ]
  oobErrVar <- rowMeans(oobErrVar)
  varimport <- cbind(varible = seq(p), increased_error = oobErrVar)
  rownames(varimport) <- colnames(X)

  varimport <- list(varImp = varimport[order(oobErrVar, decreasing = T), ], type = forest$type)
  class(varimport) <- "VarImp"

  return(varimport)
}
