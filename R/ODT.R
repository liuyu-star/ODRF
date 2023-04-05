#' Classification and Regression with Oblique Decision Tree
#'
#' Classification and regression using an oblique decision tree (ODT) in which each node is split by a linear combination of predictors. Different methods are provided for selecting the linear combinations, while the splitting values are chosen by one of three criteria.
#'
#' @param formula Object of class \code{formula} with a response describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} containing variables named in the formula. If \code{data} is missing it is obtained from the current environment by \code{formula}.
#' @param X An n by d numeric matrix (preferable) or data frame.
#' @param y A response vector of length n.
#' @param split The criterion used for splitting the nodes. "entropy": information gain and "gini": gini impurity index for classification; "mse": mean square error for regression;
#' 'auto' (default): If the response in \code{data} or \code{y} is a factor, "gini" is used, otherwise regression is assumed.
#' @param lambda The argument of \code{split} is used to determine the penalty level of the partition criterion. Three options are provided including, \code{lambda=0}: no penalty; \code{lambda=2}: AIC penalty; \code{lambda='log'} (Default): BIC penalty. In Addition, lambda can be any value from 0 to n (training set size).
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictors in the split node.
#' including \itemize{
#' \item{"RotMatPPO": projection pursuit optimization model (\code{\link{PPO}}), see \code{\link{RotMatPPO}} (default, model="PPR").}
#' \item{"RotMatRF": single feature similar to CART, see \code{\link{RotMatRF}}.}
#' \item{"RotMatRand": random rotation, see \code{\link{RotMatRand}}.}
#' \item{"RotMatMake": users can define this function, for details see \code{\link{RotMatMake}}.}
#' }
#' @param FunDir The path to the \code{function} of the user-defined \code{NodeRotateFun} (default current working directory).
#' @param paramList List of parameters used by the functions \code{NodeRotateFun}. If left unchanged, default values will be used, for details see \code{\link[ODRF]{defaults}}.
#' @param MaxDepth The maximum depth of the tree (default \code{Inf}).
#' @param numNode Number of nodes that can be used by the tree (default \code{Inf}).
#' @param MinLeaf Minimal node size (Default 10).
#' @param Levels The category label of the response variable when \code{split} is not equal to 'mse'.
#' @param subset An index vector indicating which rows should be used. (NOTE: If given, this argument must be named.)
#' @param weights Vector of non-negative observational weights; fractional weights are allowed (default NULL).
#' @param na.action A function to specify the action to be taken if NAs are found. (NOTE: If given, this argument must be named.)
#' @param catLabel A category labels of class \code{list} in predictors. (default NULL, for details see Examples)
#' @param Xcat A class \code{vector} is used to indicate which predictor is the categorical variable. The default Xcat=0 means that no special treatment is given to category variables.
#' When Xcat=NULL, the predictor x that satisfies the condition "\code{(length(table(x))<10) & (length(x)>20)}" is judged to be a category variable.
#' @param Xscale Predictor standardization methods. " Min-max" (default), "Quantile", "No" denote Min-max transformation, Quantile transformation and No transformation respectively.
#' @param TreeRandRotate If or not to randomly rotate the training data before building the tree (default FALSE, see \code{\link[ODRF]{RandRot}}).
#' @param ... Optional parameters to be passed to the low level function.
#'
#' @return An object of class ODT containing a list of components::
#' \itemize{
#' \item{\code{call}: The original call to ODT.}
#' \item{\code{terms}: An object of class \code{c("terms", "formula")} (see \code{\link{terms.object}}) summarizing the formula. Used by various methods, but typically not of direct relevance to users.}
#' \item{\code{split}, \code{Levels} and \code{NodeRotateFun} are important parameters for building the tree.}
#' \item{\code{predicted}: the predicted values of the training data.}
#' \item{\code{projections}: Projection direction for each split node.}
#' \item{\code{paramList}: Parameters in a named list to be used by \code{NodeRotateFun}.}
#' \item{\code{data}: The list of data related parameters used to build the tree.}
#' \item{\code{tree}: The list of tree related parameters used to build the tree.}
#' \item{\code{structure}: A set of tree structure data records.
#' \itemize{
#' \item{\code{nodeRotaMat}: Record the split variables (first column), split node serial number (second column) and rotation direction (third column) for each node. (The first column and the third column are 0 means leaf nodes)}
#' \item{\code{nodeNumLabel}: Record each leaf node's category for classification or predicted value for regression (second column is data size). (Each column is 0 means it is not a leaf node)}
#' \item{\code{nodeCutValue}: Record the split point of each node. (0 means leaf nodes)}
#' \item{\code{nodeCutIndex}: Record the index values of the partitioning variables selected based on the partition criterion \code{split}.}
#' \item{\code{childNode}: Record the number of child nodes after each splitting.}
#' \item{\code{nodeDepth}: Record the depth of the tree where each node is located.}
#' }}
#' }
#'
#' @seealso \code{\link{online.ODT}} \code{\link{prune.ODT}} \code{\link{as.party}} \code{\link{predict.ODT}} \code{\link{print.ODT}} \code{\link{plot.ODT}} \code{\link{plot_ODT_depth}}
#'
#' @author Yu Liu and Yingcun Xia
#' @references Zhan, H., Liu, Y., & Xia, Y. (2022). Consistency of The Oblique Decision Tree and Its Random Forest. arXiv preprint arXiv:2211.12653.
#' @keywords tree
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
#' tree <- ODT(Density ~ ., train_data,
#'   split = "mse",
#'   NodeRotateFun = "RotMatPPO", paramList = list(model = "Log", dimProj = "Rand")
#' )
#' pred <- predict(tree, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#'
#' # Projection analysis of the oblique decision tree.
#' data(iris)
#' tree <- ODT(Species ~ ., data = iris, split="gini",
#'             paramList = list(model = "PPR", numProj = 1))
#' print(round(tree[["projections"]],3))
#'
#' ### Train ODT on one-of-K encoded categorical data ###
#' # Note that the category variable must be placed at the beginning of the predictor X
#' # as in the following example.
#' set.seed(22)
#' Xcol1 <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' Xcol2 <- sample(c("1", "2", "3", "4", "5"), 100, replace = TRUE)
#' Xcon <- matrix(rnorm(100 * 3), 100, 3)
#' X <- data.frame(Xcol1, Xcol2, Xcon)
#' Xcat <- c(1, 2)
#' catLabel <- NULL
#' y <- as.factor(sample(c(0, 1), 100, replace = TRUE))
#' tree <- ODT(X, y, split = "entropy", Xcat = NULL)
#' head(X)
#' #>   Xcol1 Xcol2          X1         X2          X3
#' #> 1     B     5 -0.04178453  2.3962339 -0.01443979
#' #> 2     A     4 -1.66084623 -0.4397486  0.57251733
#' #> 3     B     2 -0.57973333 -0.2878683  1.24475578
#' #> 4     B     1 -0.82075051  1.3702900  0.01716528
#' #> 5     C     5 -0.76337897 -0.9620213  0.25846351
#' #> 6     A     5 -0.37720294 -0.1853976  1.04872159
#'
#' # one-of-K encode each categorical feature and store in X1
#' numCat <- apply(X[, Xcat, drop = FALSE], 2, function(x) length(unique(x)))
#' # initialize training data matrix X
#' X1 <- matrix(0, nrow = nrow(X), ncol = sum(numCat))
#' catLabel <- vector("list", length(Xcat))
#' names(catLabel) <- colnames(X)[Xcat]
#' col.idx <- 0L
#' # convert categorical feature to K dummy variables
#' for (j in seq_along(Xcat)) {
#'   catMap <- (col.idx + 1):(col.idx + numCat[j])
#'   catLabel[[j]] <- levels(as.factor(X[, Xcat[j]]))
#'   X1[, catMap] <- (matrix(X[, Xcat[j]], nrow(X), numCat[j]) ==
#'     matrix(catLabel[[j]], nrow(X), numCat[j], byrow = TRUE)) + 0
#'   col.idx <- col.idx + numCat[j]
#' }
#' X <- cbind(X1, X[, -Xcat])
#' colnames(X) <- c(paste(rep(seq_along(numCat), numCat), unlist(catLabel),
#'   sep = "."
#' ), "X1", "X2", "X3")
#'
#' # Print the result after processing of category variables.
#' head(X)
#' #>   1.A 1.B 1.C 2.1 2.2 2.3 2.4 2.5          X1         X2          X3
#' #> 1   0   1   0   0   0   0   0   1 -0.04178453  2.3962339 -0.01443979
#' #> 2   1   0   0   0   0   0   1   0 -1.66084623 -0.4397486  0.57251733
#' #> 3   0   1   0   0   1   0   0   0 -0.57973333 -0.2878683  1.24475578
#' #> 4   0   1   0   1   0   0   0   0 -0.82075051  1.3702900  0.01716528
#' #> 5   0   0   1   0   0   0   0   1 -0.76337897 -0.9620213  0.25846351
#' #> 6   1   0   0   0   0   0   0   1 -0.37720294 -0.1853976  1.04872159
#' catLabel
#' #> $Xcol1
#' #> [1] "A" "B" "C"
#' #>
#' #> $Xcol2
#' #> [1] "1" "2" "3" "4" "5"
#'
#' tree <- ODT(X, y, split = "gini", Xcat = c(1, 2), catLabel = catLabel,NodeRotateFun = "RotMatRF")
#'
#' @import Rcpp
#' @importFrom stats model.frame model.extract model.matrix na.fail
#' @export
ODT <- function(X, ...) {
  UseMethod("ODT")
  # formula X
}

#' @rdname ODT
#' @method ODT formula
#' @aliases ODT.formula
#' @export
ODT.formula <- function(formula, data = NULL, split = "auto", lambda = "log", NodeRotateFun = "RotMatPPO", FunDir = getwd(), paramList = NULL,
                        MaxDepth = Inf, numNode = Inf, MinLeaf = 10, Levels = NULL, subset = NULL, weights = NULL, na.action = na.fail,
                        catLabel = NULL, Xcat = 0, Xscale = "Min-max", TreeRandRotate = FALSE, ...) {
  Call <- match.call()
  indx <- match(c("formula", "data", "subset", "na.action"), names(Call), nomatch = 0L) # , "weights"
  # formula=X
  if (indx[[1]] == 0) {
    stop("A 'formula' or 'X', 'y' argument is required")
  } else if (indx[[2]] == 0) {
    # stop("a 'data' argument is required")
    # data <- environment(formula)
    X <- eval(formula[[3]]) # ,envir =.BaseNamespaceEnv)
    y <- eval(formula[[2]]) # ,envir =.BaseNamespaceEnv)
    if (sum(match(class(X), c("data.frame", "matrix"), nomatch = 0L)) == 0) {
      stop("argument 'X' can only be the classes 'data.frame' or 'matrix'")
    }
    if (ncol(X) == 1) {
      stop("argument 'X' dimension must exceed 1")
    }

    if (is.null(colnames(X))) {
      colnames(X) <- paste0("X", seq_len(ncol(X)))
    }
    data <- data.frame(y, X)
    # varName <- colnames(X)
    yname <- ls(envir = .GlobalEnv)
    colnames(data) <- c(as.character(formula)[2], colnames(X))
    formula <- as.formula(paste0(as.character(formula)[2], "~."))
    Call$formula <- formula
    Call$data <- quote(data)
  } else {
    if (sum(match(class(data), c("data.frame"), nomatch = 0L)) == 0) {
      stop("argument 'data' can only be the classe 'data.frame'")
    }
    if (ncol(data) == 2) {
      stop("The predictor dimension of argument 'data' must exceed 1.")
    }

    # varName <- setdiff(colnames(data), as.character(formula)[2])
    # X <- data[, varName]
    # y <- data[, as.character(formula)[2]]
    # Call$data <- quote(data)
    yname <- colnames(data)
    data <- model.frame(formula, data, drop.unused.levels = TRUE)
    y <- data[, 1]
    X <- data[, -1]
    Call$data <- quote(data)
  }

  varName <- colnames(X)
  yname <- names(unlist(sapply(yname, function(x) grep(x, as.character(formula)[2]))))
  yname <- yname[which.max(nchar(yname))]
  if (yname != as.character(formula)[2]) {
    varName <- c(yname, varName)
  }

  if(NodeRotateFun=="RotMatRF"&&is.null(paramList$numProj)){
    if (is.null(Xcat)) {
      n=length(y)
      Xcat <- which(apply(X, 2, function(x) {
        (length(table(x)) < 10) & (n > 20)
      }))
    }
    if (sum(Xcat) > 0) {
      numCat <- apply(X[, Xcat, drop = FALSE], 2, function(x) length(table(x)))
      p=ncol(X) - sum(numCat) - length(Xcat)
      paramList$numProj <- ceiling(sqrt(p))
    }
  }


  ppTree <- ODT.compute(
    formula, Call, varName, X, y, split, lambda, NodeRotateFun, FunDir, paramList, MaxDepth, numNode,
    MinLeaf, Levels, subset, weights, na.action, catLabel, Xcat, Xscale, TreeRandRotate
  )

  nodeRotaMat <- ppTree$structure$nodeRotaMat
  cutNode <- unique(nodeRotaMat[, 2][nodeRotaMat[, 1] != 0])
  projections <- NULL
  if (length(cutNode) > 0) {
    projections <- matrix(0, length(cutNode), ppTree$data$p)
    for (cn in seq_along(cutNode)) {
      idx <- which(nodeRotaMat[, 2] == cutNode[cn])
      projections[cn, nodeRotaMat[idx, 1]] <- nodeRotaMat[idx, 3]
    }
    colnames(projections) <- ppTree$data$varName
    rownames(projections) <- paste("proj", seq_along(cutNode), sep = "")
  }

  ppTree <- c(ppTree[-(length(ppTree)-(3:0))], list(projections = projections), ppTree[(length(ppTree)-(3:0))])
  class(ppTree) <- append(class(ppTree), "ODT")

  return(ppTree)
}


#' @rdname ODT
#' @method ODT default
#' @aliases ODT.default
#' @export
ODT.default <- function(X, y, split = "auto", lambda = "log", NodeRotateFun = "RotMatPPO", FunDir = getwd(), paramList = NULL,
                        MaxDepth = Inf, numNode = Inf, MinLeaf = 10, Levels = NULL, subset = NULL, weights = NULL, na.action = na.fail,
                        catLabel = NULL, Xcat = 0, Xscale = "Min-max", TreeRandRotate = FALSE, ...) {
  Call <- match.call()
  indx <- match(c("X", "y", "subset", "na.action"), names(Call), nomatch = 0L) # , "weights"
  if (indx[[1]] == 0 || indx[[2]] == 0) {
    stop("A 'formula' or 'X', 'y' argument is required")
  } else {
    if (sum(match(class(X), c("data.frame", "matrix"), nomatch = 0L)) == 0) {
      stop("argument 'X' can only be the classes 'data.frame' or 'matrix'")
    }
    if (ncol(X) == 1) {
      stop("argument 'X' dimension must exceed 1")
    }

    if (is.null(colnames(X))) {
      colnames(X) <- paste0("X", seq_len(ncol(X)))
    }
    data <- data.frame(y = y, X)
    varName <- colnames(X)

    formula <- y~.
    Call$formula <- formula
    Call$data <- quote(data)
    Call$X <- NULL
    Call$y <- NULL
  }

  if(NodeRotateFun=="RotMatRF"&&is.null(paramList$numProj)){
    if (is.null(Xcat)) {
      n=length(y)
      Xcat <- which(apply(X, 2, function(x) {
        (length(table(x)) < 10) & (n > 20)
      }))
    }
    if (sum(Xcat) > 0) {
      numCat <- apply(X[, Xcat, drop = FALSE], 2, function(x) length(table(x)))
      p=ncol(X) - sum(numCat) - length(Xcat)
      paramList$numProj <- ceiling(sqrt(p))
    }
  }


  ppTree <- ODT.compute(
    formula, Call, varName, X, y, split, lambda, NodeRotateFun, FunDir, paramList, MaxDepth, numNode,
    MinLeaf, Levels, subset, weights, na.action, catLabel, Xcat, Xscale, TreeRandRotate
  )

  nodeRotaMat <- ppTree$structure$nodeRotaMat
  cutNode <- which(ppTree$structure$nodeCutValue!= 0)# unique(nodeRotaMat[nodeRotaMat[, 1] != 0, 2])
  projections <- NULL
  if (length(cutNode) > 0) {
    projections <- matrix(0, length(cutNode), ppTree$data$p)
    for (cn in seq_along(cutNode)) {
      idx <- which(nodeRotaMat[, 2] == cutNode[cn])
      projections[cn, nodeRotaMat[idx, 1]] <- nodeRotaMat[idx, 3]
    }
    colnames(projections) <- ppTree$data$varName
    rownames(projections) <- paste("proj", seq_along(cutNode), sep = "")
  }

  ppTree <- c(ppTree[-(length(ppTree)-(3:0))], list(projections = projections), ppTree[(length(ppTree)-(3:0))])
  class(ppTree) <- append(class(ppTree), "ODT")

  return(ppTree)
}

#' @keywords internal
#' @noRd
ODT.compute <- function(formula, Call, varName, X, y, split, lambda, NodeRotateFun, FunDir, paramList, MaxDepth, numNode,
                        MinLeaf, Levels, subset, weights, na.action, catLabel, Xcat, Xscale, TreeRandRotate) {
  if (is.factor(y) && (split == "auto")) {
    split <- "gini"
    warning("You are creating a tree for classification")
  }
  if (is.numeric(y) && (split == "auto")) {
    split <- "mse"
    warning("You are creating a tree for regression")
  }
  if (is.factor(y) && (split == "mse")) {
    stop(paste0("When ", formula[[2]], " is a factor type, 'split' cannot take 'mse'."))
  }

  MinLeaf <- (MinLeaf == 1) + MinLeaf
  # if (MinLeaf == 5) {
  #  MinLeaf <- ifelse(split == "mse", 10, 5)
  # }

  if ((!NodeRotateFun %in% ls("package:ODRF")) && (!NodeRotateFun %in% ls(envir = .GlobalEnv))) {
    source(paste0(FunDir, "/", NodeRotateFun, ".R"))
  }

  FUN <- match.fun(NodeRotateFun, descend = TRUE)
  # split0 <- strsplit(split, split = "")[[1]][1]

  n <- length(y)
  p <- ncol(X)
  yname <- NULL
  if (length(varName) > p) {
    yname <- varName[1]
    varName <- varName[-1]
  }

  if (split != "mse") {
    if (is.null(Levels)) {
      Levels <- levels(as.factor(y))
      y <- as.integer(as.factor(y))
    }
    maxLabel <- length(Levels)
    if (length(Levels) == 1) {
      stop("the number of factor levels of response variable must be greater than one")
    }
  } else {
    y <- c(y)
    maxLabel <- 0
  }

  if (is.null(Xcat)) {
    Xcat <- which(apply(X, 2, function(x) {
      (length(table(x)) < 10) & (n > 20)
    }))
  }

  numCat <- 0
  if ((sum(Xcat) > 0) && (is.null(catLabel))) {
    warning(paste0("The categorical variable ", paste(Xcat, collapse = ", "), " has been transformed into an one-of-K encode variables!"))
    numCat <- apply(X[, Xcat, drop = FALSE], 2, function(x) length(unique(x)))
    X1 <- matrix(0, nrow = n, ncol = sum(numCat)) # initialize training data matrix X
    catLabel <- vector("list", length(Xcat))
    names(catLabel) <- colnames(X)[Xcat]
    col.idx <- 0L
    # one-of-K encode each categorical feature and store in X
    for (j in seq_along(Xcat)) {
      catMap <- (col.idx + 1L):(col.idx + numCat[j])
      # convert categorical feature to K dummy variables
      catLabel[[j]] <- levels(as.factor(X[, Xcat[j]]))
      X1[, catMap] <- (matrix(X[, Xcat[j]], n, numCat[j]) == matrix(catLabel[[j]], n, numCat[j], byrow = TRUE)) + 0
      col.idx <- col.idx + numCat[j]
    }
    X <- cbind(X1, X[, -Xcat])
    varName <- c(paste(rep(seq_along(numCat), numCat), unlist(catLabel), sep = "."), varName[-Xcat])
    rm(X1)
    p <- ncol(X)
  }
  X <- as.matrix(X)
  colnames(X) <- varName
  if (!is.numeric(X)){
    X=apply(X, 2, as.numeric)
  }


  if (all(apply(X, 2, is.character)) && (sum(Xcat) > 0)) {
    stop("The training data 'data' contains categorical variables, so that 'Xcal=NULL' can be automatically transformed into an one-of-K encode variables.")
  }

  # address na values.
  data <- data.frame(y, X)
  if (any(is.na(data))) {
    warning("NA values exist in data frame")
  }

  Call0 <- Call
  colnames(data) <- c(as.character(formula)[2], varName)
  if (!is.null(yname)) {
    colnames(data)[1] <- yname
    temp <- model.frame(formula, data, drop.unused.levels = TRUE)
    Terms <- attr(temp, "terms")

    colnames(data)[1] <- "y" # as.character(formula)[2]
    Call0$formula[[2]] <- quote(y)
  }

  indx <- match(c("formula", "data", "subset", "na.action"), names(Call0), nomatch = 0L)
  temp <- Call0[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  temp$drop.unused.levels <- TRUE
  temp <- eval(temp) # , parent.frame())
  Terms0 <- attr(temp, "terms")
  if (is.null(yname)) {
    Terms <- Terms0
    Call <- Call0
  }

  # data=model.frame(formula, data, drop.unused.levels = TRUE)
  # y <- data[,1]
  # X <- data[,-1]
  y <- c(model.extract(temp, "response"))
  X <- model.matrix(Terms0, temp)
  int <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  if (int > 0) {
    X <- X[, -int, drop = FALSE]
  }
  n <- length(y)
  p <- ncol(X)
  if (!is.integer(y) && split != "mse") {
    y <- as.integer(as.factor(y))
  }

  rm(data)


  # weights=c(weights,paramList$weights)
  if (!is.null(subset)) {
    weights <- weights[subset]
  }
  if (!is.null(weights)) {
    X <- X * matrix(weights, n, p)
  }

  # Variable scaling.
  minCol <- NULL
  maxminCol <- NULL
  if (Xscale != "No") {
    indp <- (sum(numCat) + 1):p
    if (Xscale == "Min-max") {
      minCol <- apply(X[, indp], 2, min)
      maxminCol <- apply(X[, indp], 2, function(x) {
        max(x) - min(x)
      })
    }
    if (Xscale == "Quantile") {
      minCol <- apply(X[, indp], 2, quantile, 0.05)
      maxminCol <- apply(X[, indp], 2, function(x) {
        quantile(x, 0.95) - quantile(x, 0.05)
      })
    }
    X[, indp] <- (X[, indp] - matrix(minCol, n, length(indp), byrow = T)) / matrix(maxminCol, n, length(indp), byrow = T)
  }

  # rotate the data?
  rotdims <- NULL
  rotmat <- NULL
  if (TreeRandRotate) {
    # if p > 1000 then rotate only a random subset of 1000 of the dimensions
    if (p > 1000L) {
      rotmat <- RandRot(1000L)
      rotdims <- sample.int(p, 1000L)
      X[, rotdims] <- X[, rotdims] %*% rotmat
    } else {
      rotdims <- 1:p
      rotmat <- RandRot(p)
      X <- X %*% rotmat
    }
  }

  dimProj <- paramList$dimProj
  numProj <- paramList$numProj
  paramList <- defaults(paramList, split, p, weights, catLabel)

  if ((split == "mse") && (!paramList$model %in% c("PPR", "Rand", "Log"))) {
    stop(paste0("'model = ", paramList$model, "' can only be used for classification"))
  }

  if (is.infinite(MaxDepth)) {
    numNode <- min(numNode, sum(2^(0:ceiling(log2(n / MinLeaf)))))
  } else {
    MaxDepth <- min(MaxDepth, ceiling(log2(n / MinLeaf)), n / MinLeaf - 1)
    numNode <- min(numNode, sum(2^(0:MaxDepth)))
  }

  nodeXIndx <- vector("list", numNode + 1)
  nodeXIndx[[1]] <- 1:n

  if (split != "mse") {
    nodeNumLabel <- matrix(0, 0, maxLabel)
    colnames(nodeNumLabel) <- Levels
    sl <- seq(maxLabel)
  } else {
    nodeNumLabel <- matrix(0, 0, 2)
    colnames(nodeNumLabel) <- c("prediction", "number")
  }

  nodeRotaMat <- NULL
  nodeDepth <- rep(1, numNode)
  nodeCutValue <- rep(0, numNode)
  nodeCutIndex <- nodeCutValue
  # nodeLabel = nodeCutValue;
  # nodeNumLabel = nodeCutValue;
  childNode <- nodeCutValue
  nodeLR <- nodeCutValue

  # start create Oblique Decision Tree.
  ##############################################################################
  currentNode <- 1
  freeNode <- 2
  while (!is.null(nodeXIndx[[currentNode]])) {
    if ((length(unique(y[nodeXIndx[[currentNode]]])) == 1) ||
      (length(nodeXIndx[[currentNode]]) <= (2 * MinLeaf)) ||
      (nodeDepth[currentNode] >= MaxDepth) ||
      (freeNode >= numNode)) {
      if (split != "mse") {
        leafLabel <- table(Levels[c(sl, y[nodeXIndx[[currentNode]]])]) - 1
        # nodeLabel[currentNode]=names(leafLabel)[which.max(leafLabel)];
        # nodeNumLabel[currentNode]=max(leafLabel)
        # nodeNumLabel[currentNode,]=leafLabel
      } else {
        # nodeLabel[currentNode] = mean(y[nodeXIndx[[currentNode]]]);
        # nodeNumLabel[currentNode]=length(y[nodeXIndx[[currentNode]]])
        leafLabel <- c(mean(y[nodeXIndx[[currentNode]]]), length(y[nodeXIndx[[currentNode]]]))
      }
      nodeNumLabel <- rbind(nodeNumLabel, leafLabel)
      nodeRotaMat <- rbind(nodeRotaMat, c(0, currentNode, 0))
      #nodeXIndx[currentNode] <- NA

      TF <- ifelse(currentNode > 1, (nodeLR[currentNode - 1] == nodeLR[currentNode]) && (nodeCutValue[currentNode - 1] == 0), FALSE)
      if (TF && (split != "mse") && (length(unique(max.col(nodeNumLabel[currentNode - c(1, 0), ]))) == 1)) {
        idx <- which(nodeRotaMat[, 2] == nodeLR[currentNode])
        nodeRotaMat[idx[1], ] <- c(0, nodeLR[currentNode], 0)
        nodeRotaMat <- nodeRotaMat[-c(idx[-1], nrow(nodeRotaMat) - c(1, 0)), , drop = FALSE]
        nodeCutValue[nodeLR[currentNode]] <- 0
        # nodeCutValue = nodeCutValue[-(currentNode-c(1,0))]
        nodeCutIndex[nodeLR[currentNode]] <- 0
        # nodeCutIndex = nodeCutIndex[-(currentNode-c(1,0))]

        childNode[nodeLR[currentNode]] <- 0
        # childNode = childNode[-(currentNode-c(1,0))]
        idx <- which(childNode != 0)
        idx <- idx[idx > nodeLR[currentNode]]
        childNode[idx] <- childNode[idx] - 2

        nodeNumLabel[nodeLR[currentNode], ] <- colSums(nodeNumLabel[currentNode - c(1, 0), ])
        nodeNumLabel <- nodeNumLabel[-(currentNode - c(1, 0)), ]
        nodeDepth <- nodeDepth[-(currentNode - c(1, 0))]
        nodeXIndx[[nodeLR[currentNode]]] <- unlist(nodeXIndx[currentNode - c(1, 0)])
        nodeXIndx <- nodeXIndx[-(currentNode - c(1, 0))]
        nodeLR <- nodeLR[-(currentNode - c(1, 0))]
        freeNode <- freeNode - 2
        currentNode <- currentNode - 1
      } else {
        currentNode <- currentNode + 1
      }

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

    if (!NodeRotateFun %in% ls("package:ODRF", pattern = "RotMat")) {
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
        paramList$numProj <- ifelse(paramList$dimProj == "Rand", sample(floor(p / 3), 1), ceiling(p / paramList$dimProj))
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

    # bestCut <- best_cut_node(split0, rotaX, y[nodeXIndx[[currentNode]]], Wcd, MinLeaf, maxLabel)
    bestCut <- best.cut.node(rotaX, y[nodeXIndx[[currentNode]]], split, lambda, Wcd, MinLeaf, maxLabel)
    if (bestCut$BestCutVar == -1) {
      TF <- TRUE
    } else {
      Lindex <- which(rotaX[, bestCut$BestCutVar] < bestCut$BestCutVal)
      TF <- min(length(Lindex), length(nodeXIndx[[currentNode]]) - length(Lindex)) <= MinLeaf
    }
    if (TF) {
      if (split != "mse") {
        leafLabel <- table(Levels[c(sl, y[nodeXIndx[[currentNode]]])]) - 1
      } else {
        leafLabel <- c(mean(y[nodeXIndx[[currentNode]]]), length(y[nodeXIndx[[currentNode]]]))
      }
      nodeNumLabel <- rbind(nodeNumLabel, leafLabel)
      nodeRotaMat <- rbind(nodeRotaMat, c(0, currentNode, 0))
      #nodeXIndx[currentNode] <- NA

      TF <- ifelse(currentNode > 1, (nodeLR[currentNode - 1] == nodeLR[currentNode]) && (nodeCutValue[currentNode - 1] == 0), FALSE)
      if (TF && (split != "mse") && (length(unique(max.col(nodeNumLabel[currentNode - c(1, 0), ]))) == 1)) {
        idx <- which(nodeRotaMat[, 2] == nodeLR[currentNode])
        nodeRotaMat[idx[1], ] <- c(0, nodeLR[currentNode], 0)
        nodeRotaMat <- nodeRotaMat[-c(idx[-1], nrow(nodeRotaMat) - c(1, 0)), , drop = FALSE]
        nodeCutValue[nodeLR[currentNode]] <- 0
        # nodeCutValue = nodeCutValue[-(currentNode-c(1,0))]
        nodeCutIndex[nodeLR[currentNode]] <- 0
        # nodeCutIndex = nodeCutIndex[-(currentNode-c(1,0))]

        childNode[nodeLR[currentNode]] <- 0
        # childNode = childNode[-(currentNode-c(1,0))]
        idx <- which(childNode != 0)
        idx <- idx[idx > nodeLR[currentNode]]
        childNode[idx] <- childNode[idx] - 2

        nodeNumLabel[nodeLR[currentNode], ] <- colSums(nodeNumLabel[currentNode - c(1, 0), ])
        nodeNumLabel <- nodeNumLabel[-(currentNode - c(1, 0)), ]
        nodeDepth <- nodeDepth[-(currentNode - c(1, 0))]
        nodeXIndx[[nodeLR[currentNode]]] <- unlist(nodeXIndx[currentNode - c(1, 0)])
        nodeXIndx <- nodeXIndx[-(currentNode - c(1, 0))]
        nodeLR <- nodeLR[-(currentNode - c(1, 0))]
        freeNode <- freeNode - 2
        currentNode <- currentNode - 1
      } else {
        currentNode <- currentNode + 1
      }

      next
    }

    sparseM <- sparseM[sparseM[, 2] == numDr[bestCut$BestCutVar], , drop = FALSE]
    sparseM[, 2] <- currentNode
    nodeRotaMat <- rbind(nodeRotaMat, sparseM)
    rm(rotaX)
    rm(sparseM)

    nodeCutValue[currentNode] <- bestCut$BestCutVal
    nodeCutIndex[currentNode] <- bestCut$BestIndex[bestCut$BestCutVar]
    childNode[currentNode] <- freeNode

    # Lindex=which(rotaX[,bestCut$BestCutVar]<bestCut$BestCutVal);
    nodeXIndx[[freeNode]] <- nodeXIndx[[currentNode]][Lindex]
    nodeXIndx[[freeNode + 1]] <- nodeXIndx[[currentNode]][setdiff(seq_along(nodeXIndx[[currentNode]]), Lindex)]
    nodeDepth[freeNode + c(0, 1)] <- nodeDepth[currentNode] + 1
    nodeLR[freeNode + c(0, 1)] <- currentNode

    nodeNumLabel <- rbind(nodeNumLabel, 0)

    nodeXIndx[currentNode] <- NA
    freeNode <- freeNode + 2
    currentNode <- currentNode + 1
  }

  nodeDepth <- nodeDepth[1:(currentNode - 1)]
  colnames(nodeRotaMat) <- c("var", "node", "coef")
  #
  if (length(nodeDepth) > 1) {
    rownames(nodeRotaMat) <- rep(nodeDepth, table(nodeRotaMat[, 2]))
    rownames(nodeNumLabel) <- nodeDepth
  }

  predicted=rep(0,n)
  idx=which(!is.na(nodeXIndx[1:(currentNode - 1)]))
  if (split %in% c("gini","entropy")) {
    if(all(nodeCutValue == 0)){
      nodeNumLabel=matrix(nodeNumLabel,nrow = 1, ncol = length(Levels))
    }
    nodeLabel <- Levels[max.col(nodeNumLabel)]
    nodeLabel[which(rowSums(nodeNumLabel) == 0)] <- "0"
  }
  if(split == "mse"){
    nodeLabel <- nodeNumLabel[, 1]
  }
  for (i in idx) {
    predicted[nodeXIndx[[i]]]=nodeLabel[i]
    names(predicted)[nodeXIndx[[i]]]=i
  }

  ppTree <- list(call = Call, terms = Terms, split = split, Levels = Levels, NodeRotateFun = NodeRotateFun, predicted=predicted, paramList = paramList)
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
