#' using training data to update an existing \code{ODRF}.
#'
#' Update existing \code{\link{ODRF}} using batches of data to improve the model.
#'
#' @param ppForest An object of class \code{ODRF}.
#' @param X An n by d numeric matrix (preferable) or data frame is used to update the object of class \code{ODRF}.
#' @param y A response vector of length n is used to update the object of class \code{ODRF}.
#' @param weights Vector of non-negative observational weights; fractional weights are allowed (default NULL).
#' @param ... optional parameters to be passed to the low level function.
#'
#' @return The same result as \code{ODRF}.
#'
#' @seealso \code{\link{ODRF}} \code{\link{prune.ODRF}} \code{\link{online.ODT}}
#'
#' @examples
#' # Classification with Oblique Decision Random Forest
#' data(seeds)
#' set.seed(221212)
#' train <- sample(1:209, 100)
#' train_data <- data.frame(seeds[train, ])
#' test_data <- data.frame(seeds[-train, ])
#' index <- seq(floor(nrow(train_data) / 2))
#' forest <- ODRF(varieties_of_wheat ~ ., train_data[index, ],
#'   type = "i-classification",parallel = FALSE)
#' online.forest <- online(forest, train_data[-index, -8], train_data[-index, 8])
#' pred <- predict(online.forest, test_data[, -8])
#' # classification error
#' (mean(pred != test_data[, 8]))
#'
#' # Regression with Oblique Decision Random Forest
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 100)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#' index <- seq(floor(nrow(train_data) / 2))
#' forest <- ODRF(Density ~ ., train_data[index, ],
#'  type = "regression", parallel = FALSE)
#' online.forest <- online(online.forest, train_data[-index, -1],
#'  train_data[-index, 1])
#' pred <- predict(online.forest, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#'
#' @rdname online.ODRF
#' @aliases online.ODRF
#' @method online ODRF
#' @export
online.ODRF <- function(ppForest, X, y, weights = NULL, ...) {
  weights0 <- weights
  ppTrees <- ppForest$ppTrees
  Levels <- ppForest$Levels
  type <- ppForest$type
  NodeRotateFun <- ppForest$NodeRotateFun
  Call <- ppForest$call
  Terms <- ppForest$terms
  paramList <- ppForest$paramList

  ppForest <- ppForest[c(9, 10, 11)]
  ppForestVar <- c(names(ppForest$data), names(ppForest$tree), names(ppForest$forest))
  ppForest <- do.call("c", ppForest)

  for (v in seq_along(ppForestVar)) {
    assign(ppForestVar[v], ppForest[[v]])
  }
  rm(ppForest)

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
    warning("NA values exist in data matrix")
  }
  # y= data[,setdiff(colnames(data),vars[-1])]
  # X= data[,vars[-1]]
  # X=as.matrix(X)

  y <- data[, 1]
  X <- data[, -1]
  rm(data)

  # if(!is.null(subset))
  #  X=X[subset,]
  # if(!is.null(weights))
  #  X <- X * matrix(weights0,length(y),ncol(X))
  weights <- weights0


  ppForest <- list(
    call = Call, terms = Terms, type = type, Levels = NULL,
    NodeRotateFun = NodeRotateFun, paramList = paramList, oobErr = NULL, oobConfusionMat = NULL
  )
  if (type != "regression") {
    # adjust y to go from 1 to numClass if needed
    if (is.factor(y)) {
      ppForest$Levels <- levels(y)
      y <- as.integer(y)
    } else if (is.numeric(y)) {
      ppForest$Levels <- sort(unique(y))
      y <- as.integer(as.factor(y))
    } else {
      ppForest$Levels <- levels(as.factor(y))
      y <- as.integer(as.factor(y))
      # stop("Incompatible X type. y must be of type factor or numeric.")
    }

    numClass <- length(ppForest$Levels)
    classCt <- cumsum(table(y))
    if (stratify) {
      Cindex <- vector("list", numClass)
      for (m in 1L:numClass) {
        Cindex[[m]] <- which(y == m)
      }
    }
  }

  Levels <- ppForest$Levels
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

  # Variable scaling.
  if (Xscale != "No") {
    indp <- (numCat + 1):p
    X[, indp] <- (X[, indp] - matrix(minCol, n, length(indp), byrow = T)) /
      matrix(maxminCol, n, length(indp), byrow = T)
  }

  ppForest$data <- list(
    subset = subset, weights = weights, na.action = na.action, n = n, p = p, varName = varName,
    Xscale = Xscale, minCol = minCol, maxminCol = maxminCol, Xcat = Xcat, catLabel = catLabel
  )
  ppForest$tree <- list(FunDir = FunDir, MaxDepth = MaxDepth, MinLeaf = MinLeaf, numNode = numNode, TreeRandRotate = TreeRandRotate)
  ppForest$forest <- list(
    ntrees = ntrees, numOOB = numOOB, storeOOB = storeOOB, replacement = replacement, stratify = stratify,
    parallel = parallel, numCores = numCores, seed = seed
  )


  PPtree <- function(itree, ...) {
    ppTree <- ppTrees[[itree]] # [1:7]
    if ((numOOB > 0) && storeOOB) {
      ppTree <- ppTree[-(length(ppTree) - c(2, 1, 0))]
    }
    class(ppTree) <- "ODT"
    set.seed(seed + itree)

    TDindx0 <- seq(n)
    TDindx <- TDindx0
    if (replacement) {
      go <- TRUE
      while (go) {
        # make sure each class is represented in proportion to classes in initial dataset
        if (stratify && (type != "regression")) {
          if (classCt[1L] != 0L) {
            TDindx[1:classCt[1L]] <- sample(Cindex[[1L]], classCt[1L], replace = TRUE)
          }
          for (z in 2:numClass) {
            if (classCt[z - 1L] != classCt[z]) {
              TDindx[(classCt[z - 1L] + 1L):classCt[z]] <- sample(Cindex[[z]], classCt[z] - classCt[z - 1L], replace = TRUE)
            }
          }
        } else {
          TDindx <- sample(TDindx0, n, replace = TRUE)
        }
        go <- all(TDindx0 %in% TDindx)
      }
    } else {
      TDindx <- sample.int(TDindx0, n - numOOB, replace = FALSE)
    }

    weights1 <- weights[TDindx]
    ppForestT <- online(ppTree, X[TDindx, ], y[TDindx], weights1)

    if ((numOOB > 0) && storeOOB) {
      oobErr <- 1
      NTD <- setdiff(TDindx0, TDindx)
      pred <- predict(ppForestT, X[NTD, ])

      if (type != "regression") {
        oobErr <- mean(pred != Levels[y[NTD]])
      } else {
        oobErr <- mean((pred - y[NTD])^2)
      }

      ppForestT <- c(ppForestT, list(oobErr = oobErr, oobIndex = NTD, oobPred = pred))
    }

    return(ppForestT)
  }



  if (parallel) {
    # RNGkind("L'Ecuyer-CMRG")
    if (is.infinite(numCores)) {
      # Use all but 1 core if numCores=0.
      numCores <- parallel::detectCores() - 1L # logical = FALSE
    }
    numCores <- min(numCores, ntrees)
    gc()

    # cl <- parallel::makePSOCKcluster(num.cores)
    # library("ODRF1")
    # library(foreach)
    # foreach::registerDoSEQ()
    cl <- parallel::makeCluster(numCores, type = ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK"))
    chunks <- parallel::clusterSplit(cl, seq_len(ntrees))
    doParallel::registerDoParallel(cl, numCores)
    # set.seed(seed)
    ppForestT <- foreach::foreach(
      icore = seq_along(chunks), .combine = list, .multicombine = TRUE,
      .packages = "ODRF", .noexport = "ppForest"
    ) %dopar% {
      lapply(chunks[[icore]], PPtree)
    }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)

    # do.call(rbind.fill,list1)
    ppForest$ppTrees <- do.call("c", ppForestT)
    # ppForest$ppTrees=NULL
    # for (i in 1:numCores) {
    #  ppForest$ppTrees=c(ppForest$ppTrees,ppForestT[[i]])
    # }
  } else {
    # Use just one core.
    ppForest$ppTrees <- lapply(1:ntrees, PPtree)
  }


  ####################################
  if ((numOOB > 0) && storeOOB) {
    oobVotes <- matrix(NA, n, ntrees)
    for (t in seq_len(ntrees)) {
      oobVotes[ppForest$ppTrees[[t]]$oobIndex, t] <- ppForest$ppTrees[[t]]$oobPred
    }
    idx <- which(rowSums(is.na(oobVotes)) < ntrees)
    oobVotes <- oobVotes[idx, , drop = FALSE]
    yy <- y[idx]

    if (type != "regression") {
      ny <- length(yy)
      nC <- numClass
      weights <- rep(1, ny * ntrees)
      Votes <- factor(c(t(oobVotes)), levels = Levels)
      Votes <- as.integer(Votes) + nC * rep(0:(ny - 1), rep(ntrees, ny))
      Votes <- aggregate(c(rep(0, ny * nC), weights), by = list(c(1:(ny * nC), Votes)), sum)[, 2]

      prob <- matrix(Votes, ny, nC, byrow = TRUE)
      pred <- max.col(prob) ## "random"
      oobPred <- Levels[pred]
      ppForest$oobErr <- mean(Levels[yy] != oobPred)

      # oobPred=rep(NA,noob)
      # for (i in 1:noob) {
      #  oobTable = table(oobVotes[i,])
      #  oobPred[i]=names(oobTable)[which.max(oobTable)];
      # }

      # oobErr=mean(oobPred!=Levels[y[idx]]);
      XConfusionMat <- table(oobPred, Levels[yy])
      class_error <- (rowSums(XConfusionMat) - diag(XConfusionMat)) / rowSums(XConfusionMat)
      XConfusionMat <- cbind(XConfusionMat, class_error)
      ppForest$oobConfusionMat <- XConfusionMat
    } else {
      oobPred <- rowMeans(oobVotes, na.rm = TRUE)
      ppForest$oobErr <- mean((oobPred - yy)^2) / mean((yy - mean(y))^2)
    }
  }

  # class(ppTree) <- append(class(ppTree),"ODRF")
  class(ppForest) <- "ODRF"
  return(ppForest)
}
