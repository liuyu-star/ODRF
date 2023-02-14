#' using training data to update an existing \code{ODT}.
#'
#' Update existing \code{\link{ODT}} using batches of data to improve the model.
#'
#' @param obj an object of class \code{ODT}.
#' @param X An n by d numeric matrix (preferable) or data frame is used to update the object of class \code{ODT}.
#' @param y A response vector of length n is used to update the object of class \code{ODT}.
#' @param weights Vector of non-negative observational weights; fractional weights are allowed (default NULL).
#' @param ... optional parameters to be passed to the low level function.
#'
#' @return The same result as \code{ODT}.
#'
#' @seealso \code{\link{ODT}} \code{\link{prune.ODT}} \code{\link{online.ODRF}}
#'
#' @examples
#' # Classification with Oblique Decision Tree
#' data(seeds)
#' set.seed(221212)
#' train <- sample(1:209, 100)
#' train_data <- data.frame(seeds[train, ])
#' test_data <- data.frame(seeds[-train, ])
#' index <- seq(floor(nrow(train_data) / 2))
#' tree <- ODT(varieties_of_wheat ~ ., train_data[index, ], split = "gini")
#' online_tree <- online(tree, train_data[-index, -8], train_data[-index, 8])
#' pred <- predict(online_tree, test_data[, -8])
#' # classification error
#' (mean(pred != test_data[, 8]))
#'
#' # Regression with Oblique Decision Tree
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 100)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#' index <- seq(floor(nrow(train_data) / 2))
#' tree <- ODT(Density ~ ., train_data[index, ], split = "mse")
#' online_tree <- online(tree, train_data[-index, -1], train_data[-index, 1])
#' pred <- predict(online_tree, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#'
#' @keywords tree online
#' @rdname online.ODT
#' @aliases online.ODT
#' @method online ODT
#' @export
online.ODT <- function(obj, X = NULL, y = NULL, weights = NULL, ...) {
  ppTree <- obj
  rm(obj)
  if(length(ppTree[["structure"]][["nodeDepth"]])==1)
    stop("No tree structure to use 'online'!")
  weights0 <- weights
  Call <- ppTree$call
  Terms <- ppTree$terms
  split <- ppTree$split
  Levels <- ppTree$Levels
  NodeRotateFun <- ppTree$NodeRotateFun
  paramList <- ppTree$paramList
  ppTree <- ppTree[-seq(5)]
  if ("projections" %in% names(ppTree)) {
    projections <- ppTree$projections
    ppTree <- ppTree[-1]
  }
  ppTree <- ppTree[-1]

  subset <- weights <- na.action <- n <- p <- varName <- Xscale <- minCol <- maxminCol <- NULL
  Xcat <- catLabel <- TreeRandRotate <- rotdims <- rotmat <- NULL
  lambda <- FunDir <- MaxDepth <- MinLeaf <- numNode <- NULL
  nodeRotaMat <- nodeNumLabel <- nodeCutValue <- nodeCutIndex <- childNode <- nodeDepth <- NULL

  ppTreeVar <- c(names(ppTree$data), names(ppTree$tree), names(ppTree$structure))
  ppTree <- do.call("c", ppTree)
  # nameTree <- names(ppTree)

  # env <- new.env()
  for (v in seq_along(ppTreeVar)) {
    assign(ppTreeVar[v], ppTree[[v]]) # ,envir = env)
  }
  rm(ppTree)

  # vars=all.vars(Terms)
  if (sum(Xcat) > 0 && is.null(catLabel)) {
    # vars=vars[-(1+seq(length(unlist(catLabel))))]
    # }else{
    stop("'Xcat!=0' however 'catLabel' does not exist!")
  }

  # if(is.null(data)){
  data <- data.frame(y, X)
  #  colnames(data)=vars
  # }
  # address na values.
  if (any(is.na(data))) {
    data <- na.action(data.frame(data))
    warning("NA values exist in data matrix 'X'")
  }
  # y= data[,setdiff(colnames(data),vars[-1])]
  # X= data[,vars[-1]]
  # X=as.matrix(X)

  y <- data[, 1]
  X <- data[, -1]
  rm(data)

  n <- length(y)
  p <- ncol(X)

  if ((!NodeRotateFun %in% ls("package:ODRF")) && (!NodeRotateFun %in% ls(envir = .GlobalEnv))) {
    source(paste0(FunDir, "/", NodeRotateFun, ".R"))
  }

  # get()
  FUN <- match.fun(NodeRotateFun, descend = TRUE)
  #method0 <- strsplit(split, split = "")[[1]][1]


  if (split != "mse") {
    if (!is.integer(y)) {
      y <- as.integer(as.factor(y))
    }
    maxLabel <- length(Levels)

    nodeLabel <- Levels[max.col(nodeNumLabel)] ## "random"
    nodeLabel[which(rowSums(nodeNumLabel) == 0)] <- 0
    sl <- seq(maxLabel)
  } else {
    y <- c(y)
    maxLabel <- 0

    nodeLabel <- as.character(nodeNumLabel[, 1])
  }


  if (all(nodeCutValue == 0)) {
    Nodes <- rep(1, n)
  } else {
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

      X <- cbind(X1, apply(X[, -Xcat], 2, as.numeric))
      p <- ncol(X)
      numCat <- length(unlist(catLabel))
      rm(X1)
      rm(Xj)
    }

    X <- as.matrix(X)
    colnames(X) <- varName

    if (!is.null(subset)) {
      X <- X[subset, ]
    }

    # weights=c(weights,paramList$weights)
    if (!is.null(weights)) {
      X <- X * matrix(weights0, length(y), ncol(X))
    }
    weights <- weights0

    # Variable scaling.
    if (Xscale != "No") {
      indp <- (numCat + 1):p
      X[, indp] <- (X[, indp] - matrix(minCol, n, length(indp), byrow = T)) /
        matrix(maxminCol, n, length(indp), byrow = T)
    }

    if (TreeRandRotate) {
      X[, rotdims] <- X[, rotdims, drop = FALSE] %*% rotmat
    }

    # Nodes = ppCARTNode(X, nodeRotaMat, nodeCutValue, childNode)
    # Nodes = .Call(`_ppRF_ppCARTPredict`, X, nodeRotaMat, nodeCutValue, childNode, nodeLabel)$node
    Nodes <- .Call("_ODRF_predict_ODT", PACKAGE = "ODRF", X, nodeRotaMat, nodeCutValue, childNode, nodeLabel)$node
    # , as.character(ppTree$nodeLabel)
  }

  dimProj <- paramList$dimProj
  numProj <- paramList$numProj
  paramList <- defaults(paramList, split, p, weights, catLabel)


  if (is.infinite(MaxDepth)) {
    numNode <- max(numNode, sum(2^(0:ceiling(log2(n / MinLeaf)))))
  } else {
    MaxDepth <- min(MaxDepth, ceiling(log2(n / MinLeaf)), n / MinLeaf - 1)
    numNode <- min(numNode, sum(2^(0:MaxDepth)))
  }

  numNode0 <- length(nodeCutValue)
  nodeX <- sort(unique(Nodes))
  nodeXIndx <- vector("list", numNode + 1)
  nodeXIndx[seq(numNode0)] <- NA
  for (nx in nodeX) {
    nodeXIndx[[nx]] <- which(Nodes == nx)
  }
  # cutNode=which(is.na(nodeXIndx[seq(numNode0)]))

  rep0 <- rep(0, numNode - numNode0)
  childNode0 <- childNode
  # nodeNumLabel0=c(nodeNumLabel,rep0);
  nodeNumLabel0 <- nodeNumLabel
  # names(nodeNumLabel0)=paste0("v",seq(numNode-numNode0))
  # names(nodeNumLabel0)[seq(numNode0)]=nodeLabel

  nodeDepth <- c(nodeDepth, rep0)
  nodeCutValue <- c(nodeCutValue, rep0)
  nodeCutIndex <- c(nodeCutIndex, rep0)
  childNode <- c(childNode, rep0)
  # nodeLabel = c(nodeLabel,rep0);
  # nodeNumLabel=c(nodeNumLabel,rep0);


  # start create pptree
  ##############################################################################
  currentNode <- nodeX[1]
  freeNode <- ifelse(currentNode == 1, 2, childNode[max(which(childNode[seq(currentNode)] != 0))] + 2)
  while (!is.null(nodeXIndx[[currentNode]])) {
    if (is.na(nodeXIndx[currentNode])) {
      if (nodeCutValue[currentNode] != 0) {
        freeNode <- freeNode + 2
      }
      currentNode <- currentNode + 1
      next
    }

    if ((length(unique(y[nodeXIndx[[currentNode]]])) == 1) ||
      (length(nodeXIndx[[currentNode]]) <= (MinLeaf + 1)) ||
      (nodeDepth[currentNode] >= MaxDepth) ||
      (freeNode >= numNode)) {
      nn <- length(y[nodeXIndx[[currentNode]]])
      # r=ceiling(nn*nodeNumLabel0[currentNode])
      # parentLabel=nn*nodeNumLabel0[currentNode,]
      if (split != "mse") {
        leafLabel <- table(Levels[c(sl, y[nodeXIndx[[currentNode]]])]) - 1 + nn * nodeNumLabel0[currentNode, ]
        # leafLabel = table(c(Levels[y[nodeXIndx[[currentNode]]]],rep(names(nodeNumLabel0[currentNode]),r)))
        # nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        # nodeNumLabel[currentNode]=max(leafLabel)
      } else {
        # nodeLabel[currentNode]=(r*as.numeric(names(nodeNumLabel0[currentNode]))+nn*mean(y[nodeXIndx[[currentNode]]]))/(r+nn);
        # nodeNumLabel[currentNode]=nn
        leafLabel <- c((nn * nodeNumLabel0[currentNode, 2] * nodeNumLabel0[currentNode, 1] +
          nn * mean(y[nodeXIndx[[currentNode]]])) / (nn * nodeNumLabel0[currentNode, 2] + nn))
        leafLabel <- c(leafLabel, (nn * nodeNumLabel0[currentNode, 2] + nn))
      }
      # nodeNumLabel0=rbind(nodeNumLabel0,leafLabel)

      if (currentNode > numNode0) {
        nodeRotaMat <- rbind(nodeRotaMat, c(0, currentNode, 0))
        nodeNumLabel <- rbind(nodeNumLabel, leafLabel)
      } else {
        nodeNumLabel[currentNode, ] <- leafLabel
      }

      nodeXIndx[currentNode] <- NA
      currentNode <- currentNode + 1
      next
    }

    if (!is.null(weights)) {
      Wcd <- weights[nodeXIndx[[currentNode]]]
    } else {
      Wcd <- 1
    }

    ##########################################
    if (NodeRotateFun == "RotMatMake") {
      sparseM <- RotMatMake(
        X[nodeXIndx[[currentNode]], ], y[nodeXIndx[[currentNode]]],
        paramList$RotMatFun, paramList$PPFun, FunDir, paramList
      )
    }

    if (!NodeRotateFun %in% ls("package:ODRF")) {
      paramList$X <- X[nodeXIndx[[currentNode]], ]
      paramList$y <- y[nodeXIndx[[currentNode]]]
      sparseM <- do.call(FUN, paramList)
      paramList$X <- NULL
      paramList$y <- NULL
    }

    if (NodeRotateFun %in% c("RotMatRF", "RotMatRand")) {
      sparseM <- do.call(FUN, paramList)
    }

    if (NodeRotateFun == "RotMatPPO") {
      paramList$dimProj <- dimProj
      paramList$numProj <- numProj
      if (is.null(paramList$dimProj)) {
        paramList$dimProj <- min(ceiling(length(y[nodeXIndx[[currentNode]]])^0.4), ceiling(p * 2 / 3))
      }
      if (is.null(paramList$numProj)) {
        paramList$numProj <- ifelse(paramList$dimProj == "Rand",sample(floor(p / 3), 1),ceiling(p / paramList$dimProj))
      }
      sparseM <- RotMatPPO(
        X = X[nodeXIndx[[currentNode]], ], y = y[nodeXIndx[[currentNode]]], model = paramList$model,
        split = paramList$split, weights = paramList$weights, dimProj = paramList$dimProj,
        numProj = paramList$numProj, catLabel = paramList$catLabel
      )
    }

    numDr <- unique(sparseM[, 2])
    rotaX <- matrix(0, p, length(numDr))
    for (i in seq_along(numDr)) {
      lrows <- which(sparseM[, 2] == numDr[i])
      rotaX[sparseM[lrows, 1], i] <- sparseM[lrows, 3]
    }

    ###################################################################

    rotaX <- X[nodeXIndx[[currentNode]], , drop = FALSE] %*% rotaX

    #bestCut <- best_cut_node(method0, rotaX, y[nodeXIndx[[currentNode]]], Wcd, MinLeaf, maxLabel)
    bestCut <- best.cut.node(rotaX, y[nodeXIndx[[currentNode]]], split, lambda, Wcd, MinLeaf, maxLabel)

    if (bestCut$BestCutVar == -1) {
      TF <- TRUE
    } else {
      Lindex <- which(rotaX[, bestCut$BestCutVar] < bestCut$BestCutVal)
      TF <- min(length(Lindex), length(nodeXIndx[[currentNode]]) - length(Lindex)) <= MinLeaf
    }
    if (TF) {
      nn <- length(y[nodeXIndx[[currentNode]]])
      # r=ceiling(nn*nodeNumLabel0[currentNode])
      # parentLabel=nn*nodeNumLabel0[currentNode,]
      if (split != "mse") {
        leafLabel <- table(Levels[c(sl, y[nodeXIndx[[currentNode]]])]) - 1 + nn * nodeNumLabel0[currentNode, ]
        # leafLabel = table(c(Levels[y[nodeXIndx[[currentNode]]]],rep(names(nodeNumLabel0[currentNode]),r)))
        # nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        # nodeNumLabel[currentNode]=max(leafLabel)
      } else {
        # nodeLabel[currentNode]=(r*as.numeric(names(nodeNumLabel0[currentNode]))+nn*mean(y[nodeXIndx[[currentNode]]]))/(r+nn);
        # nodeNumLabel[currentNode]=nn
        leafLabel <- c((nn * nodeNumLabel0[currentNode, 2] * nodeNumLabel0[currentNode, 1] +
          nn * mean(y[nodeXIndx[[currentNode]]])) / (nn * nodeNumLabel0[currentNode, 2] + nn))
        leafLabel <- c(leafLabel, (nn * nodeNumLabel0[currentNode, 2] + nn))
      }
      # nodeNumLabel0=rbind(nodeNumLabel0,leafLabel)

      if (currentNode > numNode0) {
        nodeRotaMat <- rbind(nodeRotaMat, c(0, currentNode, 0))
        nodeNumLabel <- rbind(nodeNumLabel, leafLabel)
      } else {
        nodeNumLabel[currentNode, ] <- leafLabel
      }

      nodeXIndx[currentNode] <- NA
      currentNode <- currentNode + 1
      next
    }

    sparseM <- sparseM[sparseM[, 2] == numDr[bestCut$BestCutVar], , drop = FALSE]
    sparseM[, 2] <- currentNode
    if (currentNode <= numNode0) {
      # r=rep(nodeNumLabel0[currentNode]/length(y[nodeXIndx[[currentNode]]]),2)
      # names(r)=rep(nodeLabel[currentNode],2)
      LRnode <- rbind(nodeNumLabel0[currentNode, ], nodeNumLabel0[currentNode, ])
      if (!is.na(childNode0[currentNode])) {
        if (split != "mse") {
          # r=rep(nodeNumLabel0[currentNode],2)
          # names(r)=rep(names(nodeNumLabel0[currentNode]),2)
          LRnode <- LRnode / length(y[nodeXIndx[[currentNode]]])
        } else {
          LRnode[, 2] <- LRnode[, 2] / length(y[nodeXIndx[[currentNode]]])
        }
      }

      idx <- which(nodeRotaMat[, 2] == currentNode)
      nodeRotaMat <- rbind(nodeRotaMat[seq(idx - 1), , drop = FALSE], sparseM, nodeRotaMat[-seq(idx), , drop = FALSE])

      idx <- which(childNode0 != 0)
      idx <- idx[idx > currentNode]
      childNode0[idx] <- childNode0[idx] + 2
      childNode0[currentNode] <- freeNode

      if (freeNode <= numNode0) {
        numNode0 <- numNode0 + 2
        idx <- seq(freeNode - 1)
        nodeXIndx <- c(nodeXIndx[idx], list(NA, NA), nodeXIndx[-idx])
        childNode0 <- c(childNode0[idx], NA, NA, childNode0[-idx])

        nodeDepth <- c(nodeDepth[idx], 0, 0, nodeDepth[-idx])
        nodeCutValue <- c(nodeCutValue[idx], 0, 0, nodeCutValue[-idx])
        nodeCutIndex <- c(nodeCutIndex[idx], 0, 0, nodeCutIndex[-idx])

        # nodeLabel=c(nodeLabel[idx],0,0,nodeLabel[-idx])
        # nodeNumLabel=c(nodeNumLabel[idx],0,0,nodeNumLabel[-idx])
        # nodeNumLabel0=c(nodeNumLabel0[idx],r,nodeNumLabel0[-idx])
        nodeNumLabel <- rbind(nodeNumLabel[idx, ], 0, 0, nodeNumLabel[-idx, , drop = FALSE])
        nodeNumLabel0 <- rbind(nodeNumLabel0[idx, ], LRnode, nodeNumLabel0[-idx, , drop = FALSE])

        idx <- 1:max(which(nodeRotaMat[, 2] == (freeNode - 1)))
        nodeRotaMat <- rbind(nodeRotaMat[idx, , drop = FALSE], 0, 0, nodeRotaMat[-idx, , drop = FALSE])
        nodeRotaMat[-idx, 2] <- c(freeNode, freeNode + 1, nodeRotaMat[-c(1, 2, idx + 2), 2] + 2)
      } else {
        # nodeNumLabel0[freeNode+c(0,1)]=r
        # names(nodeNumLabel0)[freeNode+c(0,1)]=names(r)
        nodeNumLabel0 <- rbind(nodeNumLabel0, LRnode)
        # nodeNumLabel=rbind(nodeNumLabel,0,0)
      }
    } else {
      # nodeNumLabel0[freeNode+c(0,1)]=rep(nodeNumLabel0[currentNode],2)
      # names(nodeNumLabel0)[freeNode+c(0,1)]=rep(names(nodeNumLabel0[currentNode]),2)
      nodeNumLabel <- rbind(nodeNumLabel, 0)
      nodeNumLabel0 <- rbind(nodeNumLabel0, nodeNumLabel0[currentNode, ], nodeNumLabel0[currentNode, ])
      nodeRotaMat <- rbind(nodeRotaMat, sparseM)
      # childNode[currentNode]=freeNode;
    }
    rm(rotaX)
    rm(sparseM)

    nodeCutValue[currentNode] <- bestCut$BestCutVal
    nodeCutIndex[currentNode] <- bestCut$BestIndex[bestCut$BestCutVar]
    childNode[currentNode] <- freeNode

    # Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
    nodeXIndx[[freeNode]] <- nodeXIndx[[currentNode]][Lindex]
    nodeXIndx[[freeNode + 1]] <- nodeXIndx[[currentNode]][setdiff(seq_along(nodeXIndx[[currentNode]]), Lindex)]
    nodeDepth[freeNode + c(0, 1)] <- nodeDepth[currentNode] + 1

    nodeXIndx[currentNode] <- NA
    currentNode <- currentNode + 1
    freeNode <- freeNode + 2
  }
  childNode0[is.na(childNode0)] <- 0
  childNode[1:numNode0] <- childNode0
  nodeDepth <- nodeDepth[1:(currentNode - 1)]
  colnames(nodeRotaMat) <- c("var", "node", "coef")
  rownames(nodeRotaMat) <- rep(nodeDepth, table(nodeRotaMat[, 2]))
  rownames(nodeNumLabel) <- nodeDepth

  ppTree <- list(call = Call, terms = Terms, split = split, Levels = Levels, NodeRotateFun = NodeRotateFun, paramList = paramList)
  if ("projections" %in% ls(envir = .GlobalEnv)) {
    ppTree <- c(ppTree[seq(5)], list(projections = projections), ppTree[-seq(5)])
  }

  ppTree$data <- list(
    subset = subset, weights = weights, na.action = na.action, n = n, p = p, varName = varName,
    Xscale = Xscale, minCol = minCol, maxminCol = maxminCol, Xcat = Xcat, catLabel = catLabel,
    TreeRandRotate = TreeRandRotate, rotdims = rotdims, rotmat = rotmat
  )
  ppTree$tree <- list(lambda = lambda, FunDir = FunDir, MaxDepth = MaxDepth, MinLeaf = MinLeaf, numNode = numNode)
  ppTree$structure <- list(
    nodeRotaMat = nodeRotaMat, nodeNumLabel = nodeNumLabel, nodeCutValue = nodeCutValue[1:(currentNode - 1)],
    nodeCutIndex = nodeCutIndex[1:(currentNode - 1)], childNode = childNode[1:(currentNode - 1)],
    nodeDepth = nodeDepth
  )
  # class(ppTree) <- "ODT"
  class(ppTree) <- append(class(ppTree), "ODT")
  return(ppTree)
}
