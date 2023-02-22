#' predict method for ODT objects
#'
#' Prediction of ODT for an input matrix or data frame.
#'
#' @param object An object of class ODT, as that created by the function \code{\link{ODT}}.
#' @param Xnew An n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' Note that if there are NA values in the data 'Xnew', which will be replaced with the average value.
#' @param leafnode If or not output the leaf node sequence number that \code{Xnew} is partitioned. (default FALSE)
#' @param ... Arguments to be passed to methods.
#'
#' @return A vector of the following:
#' \itemize{
#' \item prediction: the prediced response of the new data.
#' \item leafnode: the leaf node sequence number that the new data is partitioned.
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
#'
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
#'
#' tree <- ODT(Density ~ ., train_data, split = "mse")
#' pred <- predict(tree, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#'
#' @importFrom stats aggregate as.formula na.action predict quantile runif
#' @keywords tree predict
#' @rdname predict.ODT
#' @aliases predict.ODT
#' @method predict ODT
#' @export
predict.ODT <- function(object, Xnew, leafnode = FALSE, ...) {

  ppTree <- object
  rm(object)

  pp <- ppTree$data$p
  if (!is.null(ppTree$data$catLabel) && (sum(ppTree$data$Xcat) > 0)) {
    pp <- pp - length(unlist(ppTree$data$catLabel)) + length(ppTree$data$Xcat)
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

  # if (!is.null(ppTree$data$subset)) {
  #  Xnew <- Xnew[ppTree$data$subset, ]
  # }
  # weights0=c(ppTree$data$weights,ppTree$paramList$weights)
  # if(!is.null(ppTree$data$weights))
  #  Xnew <- Xnew * matrix(weights,length(y),ncol(Xnew))

  p <- ncol(Xnew)
  n <- nrow(Xnew)

  if (ppTree$split != "mse") {
    nodeLabel <- colnames(ppTree$structure$nodeNumLabel)[max.col(ppTree$structure$nodeNumLabel)] ## "random"
    nodeLabel[which(rowSums(ppTree$structure$nodeNumLabel) == 0)] <- 0
  } else {
    nodeLabel <- as.character(ppTree$structure$nodeNumLabel[, 1])
  }

  if (all(ppTree$structure$nodeCutValue == 0)) {
    pred <- rep(nodeLabel, n)
  } else {
    Xcat <- ppTree$data$Xcat
    catLabel <- ppTree$data$catLabel
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

      Xnew <- cbind(Xnew1, apply(Xnew[, -Xcat], 2, as.numeric))
      p <- ncol(Xnew)
      numCat <- length(unlist(catLabel))
      rm(Xnew1)
      rm(Xnewj)
    }

    # Variable scaling.
    if (ppTree$data$Xscale != "No") {
      indp <- (numCat + 1):p
      Xnew[, indp] <- (Xnew[, indp] - matrix(ppTree$data$minCol, n, length(indp), byrow = T)) /
        matrix(ppTree$data$maxminCol, n, length(indp), byrow = T)
    }

    if (ppTree$data$TreeRandRotate) {
      Xnew[, ppTree$data$rotdims] <- Xnew[, ppTree$data$rotdims, drop = FALSE] %*% ppTree$data$rotmat
    }

    pred <- .Call("_ODRF_predict_ODT",
      PACKAGE = "ODRF", Xnew, ppTree$structure$nodeRotaMat,
      ppTree$structure$nodeCutValue, ppTree$structure$childNode, nodeLabel
    )

    if (leafnode) {
      pred <- as.integer(pred$node)
    } else {
      pred <- pred$prediction
    }
  }

  if ((ppTree$split == "mse") && (!leafnode)) {
    pred <- as.numeric(pred)
  }

  return(pred)
}
