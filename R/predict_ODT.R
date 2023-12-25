#' making predict based on ODT objects
#'
#' Prediction of ODT for an input matrix or data frame.
#'
#' @param object An object of class ODT, the same as that created by the function \code{\link{ODT}}.
#' @param Xnew An n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' Note that if there are NA values in the data 'Xnew', which will be replaced with the average value.
#' @param Xsplit Splitting variables used to construct linear model trees. The default value is NULL and is only valid when \code{split="linear"}.
#' @param type Type of prediction required. Choosing \code{"pred"} (default) gives the prediction result, and choosing \code{"leafnode"} gives the leaf node sequence number that \code{Xnew} is partitioned into.
#' For classification tasks, including classification trees (\code{split= "gini" or "entropy"}) and linear classification models (\code{split= "linear" and glmnetParList= list(family="binomial" or "multinomial")}). Setting \code{type="prob"} gives the prediction probabilities.
#' @param ... Arguments to be passed to methods.
#'
#' @return A vector of the following:
#' \itemize{
#' \item pred: the prediced response of the new data.
#' \item leafnode: the leaf node sequence number that the new data is partitioned.
#' \item prob: the prediction probabilities for classification tasks.
#' }
#'
#' @references Zhan, H., Liu, Y., & Xia, Y. (2022). Consistency of The Oblique Decision Tree and Its Random Forest. arXiv preprint arXiv:2211.12653.
#'
#' @seealso \code{\link{ODT}} \code{\link{predict.ODRF}}
#'
#' @examples
#' # Classification with Oblique Decision Tree.
#' data(seeds)
#' set.seed(221212)
#' train <- sample(1:209, 100)
#' train_data <- data.frame(seeds[train, ])
#' test_data <- data.frame(seeds[-train, ])
#' tree <- ODT(varieties_of_wheat ~ ., train_data, split = "entropy")
#' pred <- predict(tree, test_data[, -8])
#' # classification error
#' (mean(pred != test_data[, 8]))
#'
#' # Regression with Oblique Decision Tree.
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 100)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#' tree <- ODT(Density ~ ., train_data, split = "mse")
#' pred <- predict(tree, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#'
#' # Use "Z" as the splitting variable to build a linear model tree for "X" and "y".
#' set.seed(1)
#' n = 200
#' p = 10
#' q = 5
#' X = matrix(rnorm(n*p), n, p)
#' Z = matrix(rnorm(n*q), n, q)
#' y = (Z[,1] > 1)*(X[,1] - X[,2] + 2)  +
#' (Z[,1] < 1)*(Z[,2] > 0)*(X[,1] + X[,2] + 0) +
#' (Z[,1] < 1)*(Z[,2] < 0)*(X[,3] - 2)
#' my.tree <- ODT(X=X, y=y, Xsplit=Z, split = "linear",
#'                NodeRotateFun = "RotMatRF",MinLeaf = 10, MaxDepth = 5,
#'                glmnetParList=list(lambda = 0.1,family = "gaussian"))
#' (leafnode <- predict(my.tree, X, Xsplit=Z, type="leafnode"))
#'
#' y1 = (y>0)*1
#' my.tree <- ODT(X=X, y=y1, Xsplit=Z, split = "linear",
#'                NodeRotateFun = "RotMatRF",MinLeaf = 10, MaxDepth = 5,
#'                glmnetParList=list(family = "binomial"))
#' (class <- predict(my.tree, X, Xsplit=Z, type="pred"))
#' (prob <- predict(my.tree, X, Xsplit=Z, type="prob"))
#'
#' y2 = (y < -2.5)*1+(y>=-2.5&y<2.5)*2+(y>=2.5)*3
#' my.tree <- ODT(X=X, y=y2, Xsplit=Z, split = "linear",
#'                NodeRotateFun = "RotMatRF",MinLeaf = 10, MaxDepth = 5,
#'                glmnetParList=list(family = "multinomial"))
#' (prob <- predict(my.tree, X, Xsplit=Z, type="prob"))
#'
#'
#' @importFrom stats aggregate as.formula na.action predict quantile runif
#' @importFrom glmnet predict.glmnet
#' @keywords tree predict
#' @rdname predict.ODT
#' @aliases predict.ODT
#' @method predict ODT
#' @export
predict.ODT <- function(object, Xnew, Xsplit=NULL, type = c("pred","leafnode","prob")[1], ...) {

  pp <- object$data$p
  if (!is.null(object$data$catLabel) && (sum(object$data$Xcat) > 0)) {
    pp <- pp - length(unlist(object$data$catLabel)) + length(object$data$Xcat)
  }
  if (ncol(Xnew) != pp) {
    stop("The dimensions of 'Xnew' and training data do not match")
  }

  Xna <- is.na(Xnew)
  if (any(Xna)) {
    # Xnew <- object$data$na.action(data.frame(Xnew))
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
  # weights0=c(object$data$weights,object$paramList$weights)
  # if(!is.null(object$data$weights))
  #  Xnew <- Xnew * matrix(weights,length(y),ncol(Xnew))

  p <- ncol(Xnew)
  n <- nrow(Xnew)

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

    #Xnew <- cbind(Xnew1, apply(Xnew[, -Xcat], 2, as.numeric))
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
    indp <- (numCat + 1):p
    Xnew[, indp] <- (Xnew[, indp] - matrix(object$data$minCol, n, length(indp), byrow = T)) /
      matrix(object$data$maxminCol, n, length(indp), byrow = T)
  }

  if (object$data$TreeRandRotate) {
    Xnew[, object$data$rotdims] <- Xnew[, object$data$rotdims, drop = FALSE] %*% object$data$rotmat
  }

  if(is.null(Xsplit))Xsplit=Xnew
  predict_tree <- predictTree(object$structure, Xsplit, object$split, object$Levels)

  LeafNode <- predict_tree$leafnode
  nodeSeq=which(object$structure$nodeNumLabel[,2]!=0)
  if (type =="leafnode") {
    pred <- LeafNode
  }
  if (type =="pred") {
    pred <- predict_tree$prediction
    # if (!object$split %in% c("mse" ,"gini","entropy")){
    #    pred[pred<min(object$predicted)]=min(object$predicted)
    #    pred[pred>max(object$predicted)]=max(object$predicted)
    # }
    if(object$split=="linear"){
      type0 =ifelse(object$glmnetParList$family%in%c("binomial","multinomial"),"class","link")
      #nodeSeq=which(object$structure$nodeNumLabel[,2]!=0)
      pred=rep(0,n)
      for (node in nodeSeq) {
        idx=which(LeafNode==node)
        #do.call(glmnet, glmnetParList)
        pred[idx]=predict(object$structure$glmnetFit[[node]],Xnew[idx,],type=type0)
      }
    }
  }
  if (type =="prob") {
    pred=matrix(0,n,length(object$Levels))
    for (node in nodeSeq) {
      idx=which(LeafNode==node)
      #do.call(glmnet, glmnetParList)
      if(object$split=="linear"){
        pred[idx,]=predict(object$structure$glmnetFit[[node]], Xnew[idx,], type ="response")
      }else{
        pred[idx,]=object$structure$nodeNumLabel[idx,]
      }
    }
    colnames(pred)=object$Levels
  }

  return(pred)
}
