#' predict based on an ODRF object
#'
#' Prediction of ODRF for an input matrix or data frame.
#'
#' @param object An object of class ODRF, the same created by the function \code{\link{ODRF}}.
#' @param Xnew An n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' Note that if there are NA values in the data 'Xnew', which will be replaced with the average value.
#' @param type One of \code{response}, \code{prob} or \code{tree}, indicating the type of output: predicted values, matrix of class probabilities or predicted value for each tree.
#' @param weight.tree Whether to weight the tree, if \code{TRUE} then use the out-of-bag error of the tree as the weight. (default \code{FALSE})
#' @param ... Arguments to be passed to methods.
#'
#' @return A set of vectors in the following list:
#' \itemize{
#' \item \code{response}: the predicted values of the new data.
#' \item \code{prob}: matrix of class probabilities (one column for each class and one row for each input). If \code{object$split} is \code{mse}, a vector of tree weights is returned.
#' \item \code{tree}: It is a matrix where each column is a prediction for each tree.
#' }
#'
#' @seealso \code{\link{ODRF}} \code{\link{predict.ODT}}
#'
#' @references Zhan, H., Liu, Y., & Xia, Y. (2022). Consistency of The Oblique Decision Tree and Its Random Forest. arXiv preprint arXiv:2211.12653.
#'
#' @examples
#' # Classification with Oblique Decision Random Forest
#' data(seeds)
#' set.seed(221212)
#' train <- sample(1:209, 80)
#' train_data <- data.frame(seeds[train, ])
#' test_data <- data.frame(seeds[-train, ])
#' forest <- ODRF(varieties_of_wheat ~ ., train_data,
#'   split = "entropy", parallel = FALSE,ntrees = 50
#' )
#' pred <- predict(forest, test_data[, -8], weight.tree = TRUE)
#' # classification error
#' (mean(pred != test_data[, 8]))
#'
#' # Regression with Oblique Decision Random Forest
#' \donttest{
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 80)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#' forest <- ODRF(Density ~ ., train_data, split = "mse", parallel = FALSE,
#' ntrees = 50, TreeRandRotate=TRUE)
#' pred <- predict(forest, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#' }
#'
#' @keywords forest predict
#' @rdname predict.ODRF
#' @aliases predict.ODRF
#' @method predict ODRF
#' @export
predict.ODRF <- function(object, Xnew, type = "response", weight.tree = FALSE, ...) {

  pp <- object$data$p
  if (!is.null(object$data$catLabel) && (sum(object$data$Xcat) > 0)) {
    pp <- pp - length(unlist(object$data$catLabel)) + length(object$data$Xcat)
  }
  if (ncol(Xnew) != pp) {
    stop("The dimensions of 'Xnew' and training data do not match")
  }

  Xna <- is.na(Xnew)
  if (any(Xna)) {
    # Xnew <- ppTree$data$na.action(data.frame(Xnew))
    xj <- which(colSums(Xna) > 0)
    warning("There are NA values in columns ", paste(xj, collapse = ", "), " of the data 'Xnew', which will be replaced with the average value.")
    for (j in xj) {
      Xnew[Xna[, j], j] <- mean(Xnew[, j], na.rm = TRUE)
    }
  }
  Xnew <- as.matrix(Xnew)

  # if (!is.null(object$data$subset)) {
  #  Xnew <- Xnew[object$data$subset, ]
  # }

  p <- ncol(Xnew)
  n <- nrow(Xnew)
  nC <- length(object$Levels)
  ntrees <- length(object$structure)


  Xcat <- object$data$Xcat
  catLabel <- object$data$catLabel
  numCat <- 0
  if (sum(Xcat) > 0) {
    xj <- 1
    Xnew1 <- matrix(0, nrow = n, ncol = length(unlist(catLabel))) # initialize training data matrix X
    # one-of-K encode each categorical feature and store in X
    for (j in seq_along(Xcat)) {
      catMap <- which(catLabel[[j]] %in% unique(Xnew[, Xcat[j]]))
      indC <- catLabel[[j]][catMap]
      Xnewj <- (matrix(Xnew[, Xcat[j]], n, length(indC)) == matrix(indC, n, length(indC), byrow = TRUE)) + 0

      if (length(indC) > length(catLabel[[j]])) {
        Xnewj <- Xnewj[, seq_along(catLabel[[j]])]
      }

      xj1 <- xj + length(catLabel[[j]])
      Xnew1[, (xj:(xj1 - 1))[catMap]] <- Xnewj
      xj <- xj1
    }

    Xnew <- cbind(Xnew1, Xnew[, -Xcat])
    p <- ncol(Xnew)
    numCat <- length(unlist(catLabel))
    rm(Xnew1)
    rm(Xnewj)
  }
  if (!is.numeric(Xnew)){
    Xnew=apply(Xnew, 2, as.numeric)
  }

  # Variable scaling.
  if (object$data$Xscale != "No") {
    indp <- (sum(numCat) + 1):p
    Xnew[, indp] <- (Xnew[, indp] - matrix(object$data$minCol, n, length(indp), byrow = T)) /
      matrix(object$data$maxminCol, n, length(indp), byrow = T)
  }

  # Votes = matrix(0,ntrees,n);
  # oobErr = rep(0,ntrees);
  # for(i in 1:ntrees){
  # Votes[i,] = PPtreePredict(Xnew,object$trees[[i]]);
  #  oobErr[i] = object$trees[[i]]$oobErr;
  # }
  # Votes=t(sapply(seq(ntrees), function(i)PPtreePredict(Xnew,object$trees[[i]])))

  #VALUE <- rep(ifelse(split == "mse", 0, "0"), n)
  #TreePrediction <- vapply(object$structure,function(tree){predictTree(tree,Xnew,split,Levels)$prediction}, VALUE)
  #Votes <- t(TreePrediction)
  split=object$split
  Levels=object$Levels
  Rotate=object$data$TreeRandRotate
  VALUE <- rep(ifelse(split == "mse", 0, "0"), n)
  TreePrediction <- vapply(object$structure,function(tree){
    XXnew=Xnew
    if (Rotate) {
      XXnew[, tree$rotdims] <- XXnew[, tree$rotdims, drop = FALSE] %*% tree$rotmat
    }
    predictTree(tree,XXnew,split,Levels)$prediction}, VALUE)
  Votes <- t(TreePrediction)


  weights <- rep(1, ntrees)
  if (weight.tree) {
    if (object$forest$ratOOB == 0){
      warning("ratOOB=0, weight.tree = TRUE invalid")
      #stop("out-of-bag indices for each tree are not stored. ODRF must be called with storeOOB = TRUE.")
    }else{
      oobErr <- sapply(object$structure, function(trees) trees$oobErr)
      weights <- 1/(oobErr+1e-5)
    }
  }

  #weights <- weight.tree * oobErr + (!weight.tree)
  weights <- weights / sum(weights)
  if (split != "mse") {
    # prob=matrix(0,n,nC)
    # for (i in 1:n) {
    #    prob[i,]=aggregate(c(rep(0,nC),weights[,i]), by=list(c(1:nC, f_votes[,i])),sum)[,2];
    #  }
    # f_votes[i,] =sapply(f_votes[i,],function(pred)which(object$Levels%in%pred));
    # f_votes=as.numeric(as.factor(c(f_votes)))

    # Levels=levels(as.factor(c(object$Levels[1:nC],Votes[,1])))
    # Levels=max.col(matrix(Levels,nC,nC)==matrix(object$Levels,nC,nC,byrow = T))
    # prob=matrix(0,n,nC);
    # for (i in seq(n)) {
    #  prob[i,]=aggregate(c(rep(0,nC),weights), by=list(c(1:nC,as.integer(as.factor(Votes[,i])))),sum)[,2];
    # }
    # prob[,Levels]=prob

    weights <- rep(weights, n)
    Votes <- factor(c(Votes), levels = Levels)
    Votes <- as.integer(Votes) + nC * rep(0:(n - 1), rep(ntrees, n))
    Votes <- aggregate(c(rep(0, n * nC), weights), by = list(c(1:(n * nC), Votes)), sum)[, 2]

    prob <- matrix(Votes, n, nC, byrow = TRUE)
    prob <- prob / matrix(rowSums(prob), n, nC)
    colnames(prob) <- Levels

    # pred=apply(prob,1,which.max);
    pred <- max.col(prob) ## "random"
    pred <- Levels[pred]
  } else {
    prob <- weights #/ sum(weights)
    pred <- t(Votes) %*% prob
    # pred=colMeans(Votes);
    # prob=NULL
  }

  if (type == "response") {
    prediction <- pred
  }
  if (type == "prob") {
    prediction <- prob
  }
  if (type == "tree") {
    prediction <- TreePrediction
  }

  return(prediction)
}
