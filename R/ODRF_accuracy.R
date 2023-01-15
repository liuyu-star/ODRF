#' accuracy of oblique decision random forest
#'
#' Prediction accuracy of ODRF at different tree sizes.
#'
#' @param ppForest an object of class \code{ODRF}, as that created by the function \code{\link{ODRF}}.
#' @param data Training data of class \code{data.frame} in \code{\link{ODRF}} is used to calculate the OOB error.
#' @param newdata A data frame or matrix containing new data is used to calculate the test error. If it is missing, let it be \code{data}.
#'
#' @return OOB error and test error, misclassification rate (MR) for classification or mean square error (MSE) for regression.
#'
#' @seealso \code{\link{ODRF}} \code{\link{plot.ODRF_accuracy}}
#'
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train <- sample(1:569, 200)
#' train_data <- data.frame(breast_cancer[train, -1])
#' test_data <- data.frame(breast_cancer[-train, -1])
#'
#' forest <- ODRF(diagnosis ~ ., train_data, type = "i-classification", parallel = FALSE)
#' (error <- ODRF_accuracy(forest, train_data, test_data))
#'
#' @export
ODRF_accuracy <- function(ppForest, data, newdata = NULL) {
  vars <- all.vars(ppForest$terms)
  if (!all(vars[-1] %in% colnames(data))) {
    stop("The column name of 'data' does not match the training data.")
  }

  y <- data[, setdiff(colnames(data), vars[-1])]
  if (is.null(newdata)) newdata <- data
  ynew <- newdata[, setdiff(colnames(newdata), vars[-1])]
  Xnew <- newdata[, vars[-1]]
  Xnew <- as.matrix(Xnew)
  if (ppForest$type != "regression") {
    y <- factor(y, levels = ppForest$Levels)
  }

  n <- length(y)
  nt <- ntrees <- ppForest$forest$ntrees
  nC <- length(ppForest$Levels)
  ny <- length(ynew)

  treeVotes <- predict(ppForest, Xnew, type = "tree")
  err.test <- rep(0, ntrees)
  if (ppForest$type == "regression") {
    # e.0=mean((ynew-mean(y))^2)
    pred <- rowSums(treeVotes)
    err.test[nt] <- mean((ynew - pred / nt)^2) # /e.0;
    for (t in seq(nt - 1, 1)) {
      pred <- pred - treeVotes[, t + 1]
      err.test[t] <- mean((ynew - pred / t)^2) # /e.0;
    }
  } else {
    weights <- rep(1, ny * nt)
    Votes <- factor(c(t(treeVotes)), levels = ppForest$Levels)
    treeVotes <- matrix(as.integer(Votes), nt, ny)

    Votes <- c(treeVotes) + nC * rep(0:(ny - 1), rep(nt, ny))
    Votes <- aggregate(c(rep(0, ny * nC), weights), by = list(c(1:(ny * nC), Votes)), sum)[, 2]
    # Votes=aggregate(c(rep(0,ny*nC),weights), by=list(c(1:(ny*nC),Votes)),cumsum)[,2];

    # prob=matrix(Votes,ny,nC,byrow = TRUE);
    Votes <- matrix(Votes, ny, nC, byrow = TRUE)
    pred <- ppForest$Levels[max.col(Votes)] ## "random"
    err.test[nt] <- mean(ynew != pred)
    treeC <- matrix(seq(nC), ny, nC, byrow = TRUE)
    for (t in seq(nt - 1, 1)) {
      Votes <- Votes - (treeC == treeVotes[t + 1, ]) * 1
      # pred=apply(prob,1,which.max);
      pred <- ppForest$Levels[max.col(Votes)] ## "random"
      err.test[t] <- mean(ynew != pred)
    }
  }


  err.oob <- rep(0, ntrees)
  for (tt in 1:ntrees) {
    oobVotes <- matrix(NA, n, tt)
    for (t in 1:tt) {
      oobVotes[ppForest$ppTrees[[t]]$oobIndex, t] <- ppForest$ppTrees[[t]]$oobPred
    }
    idx <- which(rowSums(is.na(oobVotes)) < tt)
    oobVotes <- oobVotes[idx, , drop = FALSE]

    if (ppForest$type == "regression") {
      pred <- rowMeans(oobVotes, na.rm = TRUE)
      err <- mean((y[idx] - pred)^2) / mean((y[idx] - mean(y))^2)
    } else {
      ny <- length(y[idx])
      nt <- ncol(oobVotes)
      weights <- rep(1, ny * nt)
      Votes <- factor(c(t(oobVotes)), levels = ppForest$Levels)
      Votes <- as.integer(Votes) + nC * rep(0:(ny - 1), rep(nt, ny))
      Votes <- aggregate(c(rep(0, ny * nC), weights), by = list(c(1:(ny * nC), Votes)), sum)[, 2]

      prob <- matrix(Votes, ny, nC, byrow = TRUE)
      # pred=apply(prob,1,which.max);
      pred <- max.col(prob) ## "random"
      pred <- ppForest$Levels[pred]
      err <- mean(y[idx] != pred)
    }
    err.oob[tt] <- err
  }

  error <- list(err.oob = err.oob, err.test = err.test, type = ppForest$type)

  class(error) <- "ODRF_accuracy"
  return(error)
}
