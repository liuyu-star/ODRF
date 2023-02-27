#' Classification and Regression using Oblique Decision Random Forest
#'
#' Classification and regression implemented by the oblique decision random forest. ODRF usually produces more accurate predictions than RF, but needs longer computation time.
#'
#' @param formula Object of class \code{formula} with a response describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} containing variables named in the formula. If \code{data} is missing it is obtained from the current environment by \code{formula}.
#' @param X An n by d numeric matrix (preferable) or data frame.
#' @param y A response vector of length n.
#' @param split The criterion used for splitting the nodes. "entropy": information gain and "gini": gini impurity index for classification; "mse": mean square error for regression;
#' 'auto' (default): If the response in \code{data} or \code{y} is a factor, "gini" is used, otherwise regression is assumed.
#' @param lambda The adjustment parameter of \code{split} is used to determine whether to split or not, with the available values being 0, 1 and 'log' (Default).
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictors in the split node.
#' including \itemize{
#' \item{"RotMatPPO": projection pursuit optimization model (\code{\link{PPO}}), see \code{\link{RotMatPPO}} (default, model="PPR").}
#' \item{"RotMatRF": single feature similar to Random Forest, see \code{\link{RotMatRF}}.}
#' \item{"RotMatRand": random rotation, see \code{\link{RotMatRand}}.}
#' \item{"RotMatMake": users can define this function, for details see \code{\link{RotMatMake}}.}
#' }
#' @param FunDir The path to the \code{function} of the user-defined \code{NodeRotateFun} (default current working directory).
#' @param paramList List of parameters used by the functions \code{NodeRotateFun}. If left unchanged, default values will be used, for details see \code{\link[ODRF]{defaults}}.
#' @param ntrees The number of trees in the forest (default 100).
#' @param storeOOB If TRUE then the samples omitted during the creation of a tree are stored as part of the tree (default TRUE).
#' @param replacement if TRUE then n samples are chosen, with replacement, from training data (default TRUE).
#' @param stratify If TRUE then class sample proportions are maintained during the random sampling. Ignored if replacement = FALSE (default TRUE).
#' @param numOOB  Ratio of 'out-of-bag' (default 1/3).
#' @param parallel Parallel computing or not (default TRUE).
#' @param numCores Number of cores to be used for parallel computing (default \code{Inf}).
#' @param MaxDepth The maximum depth of the tree (default \code{Inf}).
#' @param numNode Number of nodes that can be used by the tree (default \code{Inf}).
#' @param MinLeaf Minimal node size (Default 5).
#' @param subset An index vector indicating which rows should be used. (NOTE: If given, this argument must be named.)
#' @param weights Vector of non-negative observational weights; fractional weights are allowed (default NULL).
#' @param na.action A function to specify the action to be taken if NAs are found. (NOTE: If given, this argument must be named.)
#' @param catLabel A category labels of class \code{list} in predictors. (default NULL, for details see Examples)
#' @param Xcat A class \code{vector} is used to indicate which predictor is the categorical variable, the default \code{Xcat}=0 means that no special treatment is given to category variables.
#' When Xcat=NULL, the predictor x that satisfies the condition (length(unique(x))<10) & (n>20) is judged to be a category variable.
#' @param Xscale Predictor standardization methods. " Min-max" (default), "Quantile", "No" denote Min-max transformation, Quantile transformation and No transformation respectively.
#' @param TreeRandRotate If or not to randomly rotate the Training data before building the tree (default FALSE, see \code{\link[ODRF]{RandRot}}).
#' @param ... Optional parameters to be passed to the low level function.
#'
#' @return An object of class ODRF Containing a list components:
#' \itemize{
#' \item{\code{call}: The original call to ODRF.}
#' \item{\code{terms}: An object of class \code{c("terms", "formula")} (see \code{\link{terms.object}}) summarizing the formula. Used by various methods, but typically not of direct relevance to users.}
#' \item{\code{ppTrees}: Each tree used to build the forest. \itemize{
#' \item{\code{oobErr}: 'out-of-bag' error for tree, misclassification rate (MR) for classification or mean square error (MSE) for regression.}
#' \item{\code{oobIndex}: Which training data to use as 'out-of-bag'.}
#' \item{\code{oobPred}: Predicted value for 'out-of-bag'.}
#' \item{\code{other}: For other tree related values \code{\link{ODT}}.}
#' }}
#' \item{\code{oobErr}: 'out-of-bag' error for forest, misclassification rate (MR) for classification or mean square error (MSE) for regression.}
#' \item{\code{oobConfusionMat}: 'out-of-bag' confusion matrix for forest.}
#' \item{\code{split}, \code{Levels} and \code{NodeRotateFun} are important parameters for building the tree.}
#' \item{\code{paramList}: Parameters in a named list to be used by \code{NodeRotateFun}.}
#' \item{\code{data}: The list of data related parameters used to build the forest.}
#' \item{\code{tree}: The list of tree related parameters used to build the tree.}
#' \item{\code{forest}: The list of forest related parameters used to build the forest.}
#' }
#'
#' @seealso \code{\link{online.ODRF}} \code{\link{prune.ODRF}} \code{\link{predict.ODRF}} \code{\link{print.ODRF}} \code{\link{Accuracy}} \code{\link{VarImp}}
#'
#' @author Yu Liu and Yingcun Xia
#' @references Zhan, H., Liu, Y., & Xia, Y. (2022). Consistency of The Oblique Decision Tree and Its Random Forest. arXiv preprint arXiv:2211.12653.
#' @references Tomita, T. M., Browne, J., Shen, C., Chung, J., Patsolic, J. L., Falk, B., ... & Vogelstein, J. T. (2020). Sparse projection oblique randomer forests. Journal of machine learning research, 21(104).
#' @keywords forest
#'
#' @examples
#' # Classification with Oblique Decision Randome Forest.
#' data(seeds)
#' set.seed(221212)
#' train <- sample(1:209, 80)
#' train_data <- data.frame(seeds[train, ])
#' test_data <- data.frame(seeds[-train, ])
#' forest <- ODRF(varieties_of_wheat ~ ., train_data,
#'   split = "entropy",parallel = FALSE, ntrees = 50
#' )
#' pred <- predict(forest, test_data[, -8])
#' # classification error
#' (mean(pred != test_data[, 8]))
#' \donttest{
#' # Regression with Oblique Decision Randome Forest.
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 80)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#' forest <- ODRF(Density ~ ., train_data,
#'   split = "mse", parallel = FALSE,
#'   NodeRotateFun = "RotMatPPO", paramList = list(model = "Log", dimProj = "Rand")
#' )
#' pred <- predict(forest, test_data[, -1])
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#' }
#'
#' ### Train ODRF on one-of-K encoded categorical data ###
#' set.seed(22)
#' Xcol1 <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' Xcol2 <- sample(c("1", "2", "3", "4", "5"), 100, replace = TRUE)
#' Xcon <- matrix(rnorm(100 * 3), 100, 3)
#' X <- data.frame(Xcol1, Xcol2, Xcon)
#' Xcat <- c(1, 2)
#' catLabel <- NULL
#' y <- as.factor(sample(c(0, 1), 100, replace = TRUE))
#' \donttest{
#' forest <- ODRF(y ~ X, split = "entropy", Xcat = NULL, parallel = FALSE)
#' }
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
#' # initialize training data matrix X1
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
#' \donttest{
#' forest <- ODRF(X, y,
#'   split = "gini", Xcat = c(1, 2),
#'   catLabel = catLabel, parallel = FALSE
#' )
#' }
#'
#' @import rlang
#' @importFrom glue glue
#' @importFrom lifecycle deprecated
#' @export
ODRF <- function(X, ...) {
  UseMethod("ODRF")
}


#' @rdname ODRF
#' @method ODRF formula
#' @aliases ODRF.formula
#' @export
ODRF.formula <- function(formula, data = NULL, split = "auto", lambda = "log", NodeRotateFun = "RotMatPPO", FunDir = getwd(), paramList = NULL,
                         ntrees = 100, storeOOB = TRUE, replacement = TRUE, stratify = TRUE, numOOB = 1 / 3, parallel = TRUE,
                         numCores = Inf, MaxDepth = Inf, numNode = Inf, MinLeaf = 5, subset = NULL, weights = NULL,
                         na.action = na.fail, catLabel = NULL, Xcat = 0, Xscale = "Min-max", TreeRandRotate = FALSE, ...) {
  Call <- match.call()
  indx <- match(c("formula", "data", "subset", "na.action"), names(Call), nomatch = 0L) # , "weights"
  # formula=X
  if (indx[[1]] == 0) {
    stop("A 'formula' or 'X', 'y' argument is required")
  } else if (indx[[2]] == 0) {
    # stop("a 'data' argument is required")
    # data <- environment(formula)
    X <- eval(formula[[3]])
    y <- eval(formula[[2]])
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

  forest <- ODRF.compute(
    formula, Call, varName, X, y, split, lambda, NodeRotateFun, FunDir, paramList,
    ntrees, storeOOB, replacement, stratify, numOOB, parallel,
    numCores, MaxDepth, numNode, MinLeaf, subset, weights,
    na.action, catLabel, Xcat, Xscale, TreeRandRotate
  )

  # class(forest) = append(class(forest),"ODRF.formula")
  return(forest)
}


#' @rdname ODRF
#' @method ODRF default
#' @aliases ODRF.default
#' @export
ODRF.default <- function(X, y, split = "auto", lambda = "log", NodeRotateFun = "RotMatPPO", FunDir = getwd(), paramList = NULL,
                         ntrees = 100, storeOOB = TRUE, replacement = TRUE, stratify = TRUE, numOOB = 1 / 3, parallel = TRUE,
                         numCores = Inf, MaxDepth = Inf, numNode = Inf, MinLeaf = 5, subset = NULL, weights = NULL,
                         na.action = na.fail, catLabel = NULL, Xcat = 0, Xscale = "Min-max", TreeRandRotate = FALSE, ...) {
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

  ODRF.compute(
    formula, Call, varName, X, y, split, lambda, NodeRotateFun, FunDir, paramList,
    ntrees, storeOOB, replacement, stratify, numOOB, parallel,
    numCores, MaxDepth, numNode, MinLeaf, subset, weights,
    na.action, catLabel, Xcat, Xscale, TreeRandRotate
  )
}

#' @useDynLib ODRF, .registration = TRUE
#' @import Rcpp
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @importFrom stats model.frame model.extract model.matrix na.fail
#' @keywords internal
#' @noRd
ODRF.compute <- function(formula, Call, varName, X, y, split, lambda, NodeRotateFun, FunDir, paramList,
                         ntrees, storeOOB, replacement, stratify, numOOB, parallel,
                         numCores, MaxDepth, numNode, MinLeaf, subset, weights,
                         na.action, catLabel, Xcat, Xscale, TreeRandRotate) {
  if (ntrees == 1) {
    stop("argument 'ntrees' must exceed 1")
  }
  if (is.factor(y) && (split == "auto")) {
    split <- "gini"
    warning("You are creating a forest for classification")
  }
  if (is.numeric(y) && (split == "auto")) {
    split <- "mse"
    warning("You are creating a forest for regression")
  }
  if (is.factor(y) && (split == "mse")) {
    stop(paste0("When ", formula[[2]], " is a factor type, 'split' cannot take 'regression'."))
  }
  # if (MinLeaf == 5) {
  #  MinLeaf <- ifelse(split == "mse", 5, 1)
  # }

  n <- length(y)
  p <- ncol(X)
  yname <- NULL
  if (length(varName) > p) {
    yname <- varName[1]
    varName <- varName[-1]
  }

  if (is.null(Xcat)) {
    Xcat <- which(apply(X, 2, function(x) {
      (length(unique(x)) < 10) & (n > 20)
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

  if (any(apply(X, 2, is.character)) && (sum(Xcat) > 0)) {
    stop("The training data 'data' contains categorical variables, so that 'Xcal=NULL' can be automatically transformed into an one-of-K encode variables.")
  }

  # address na values.
  data <- data.frame(y, X)
  if (any(is.na(as.list(data)))) {
    warning("NA values exist in data frame")
  }

  Call0 <- Call
  colnames(data) <- c(as.character(formula)[2], varName)
  if (!is.null(yname)) {
    colnames(data)[1] <- yname
    temp <- model.frame(formula, data, drop.unused.levels = TRUE)
    Terms <- attr(temp, "terms")

    colnames(data)[1] <- "y" # as.character(formula)[2]
    formula[[2]] <- quote(y)
    Call0$formula <- formula
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
  rm(data)

  ppForest <- list(
    call = Call, terms = Terms, split = split, Levels = NULL,
    NodeRotateFun = NodeRotateFun, paramList = paramList, oobErr = NULL, oobConfusionMat = NULL
  )

  if (split != "mse") {
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
    if (length(ppForest$Levels) == 1) {
      stop("the number of factor levels of response variable must be greater than one")
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

  # weights=c(weights,paramList$weights)
  if (!is.null(subset)) {
    weights <- weights[subset]
  }
  # if (!is.null(weights)) {
  #  X <- X * matrix(weights, n, p)
  # }

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

  ppForest$data <- list(
    subset = subset, weights = weights, na.action = na.action, n = n, p = p, varName = varName,
    Xscale = Xscale, minCol = minCol, maxminCol = maxminCol, Xcat = Xcat, catLabel = catLabel
  )
  ppForest$tree <- list(lambda = lambda, FunDir = FunDir, MaxDepth = MaxDepth, MinLeaf = MinLeaf, numNode = numNode, TreeRandRotate = TreeRandRotate)
  ppForest$forest <- list(
    ntrees = ntrees, numOOB = numOOB, storeOOB = storeOOB, replacement = replacement, stratify = stratify,
    parallel = parallel, numCores = numCores
  )

  # Weights=weights
  # vars=all.vars(Terms)
  PPtree <- function(itree, ...) {
    #set.seed(seed + itree)

    TDindx0 <- seq(n)
    TDindx <- TDindx0
    if (replacement) {
      go <- TRUE
      while (go) {
        # make sure each class is represented in proportion to classes in initial dataset
        if (stratify && (split != "mse")) {
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
      TDindx <- sample(TDindx0, n * (1 - numOOB), replace = FALSE)
    }


    # data=data.frame(y[TDindx],X[TDindx,])
    # colnames(data)=vars
    weights1 <- weights[TDindx]
    ppForestT <- ODT.compute(formula, Call0, varName,
      X = X[TDindx, ], y = y[TDindx], split, lambda, NodeRotateFun, FunDir, paramList, MaxDepth, numNode,
      MinLeaf, Levels, subset = NULL, weights = weights1, na.action = NULL, catLabel, Xcat = 0L, Xscale = "No", TreeRandRotate
    )

    if ((numOOB > 0) && storeOOB) {
      oobErr <- 1
      NTD <- setdiff(TDindx0, TDindx)
      pred <- predict(ppForestT, X[NTD, ])

      if (split != "mse") {
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
    icore <- NULL
    ppForestT <- foreach::foreach(
      icore = seq_along(chunks), .combine = list, .multicombine = TRUE, .export = c("ODT.compute"),
      .packages = c("ODRF"), .noexport = "ppForest"
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
    for (t in 1:ntrees) {
      oobVotes[ppForest$ppTrees[[t]]$oobIndex, t] <- ppForest$ppTrees[[t]]$oobPred
    }
    idx <- which(rowSums(is.na(oobVotes)) < ntrees)
    oobVotes <- oobVotes[idx, , drop = FALSE]
    yy <- y[idx]

    if (split != "mse") {
      ny <- length(yy)
      nC <- numClass
      tree_weights <- rep(1, ny * ntrees)
      Votes <- factor(c(t(oobVotes)), levels = Levels)
      Votes <- as.integer(Votes) + nC * rep(0:(ny - 1), rep(ntrees, ny))
      Votes <- aggregate(c(rep(0, ny * nC), tree_weights), by = list(c(1:(ny * nC), Votes)), sum)[, 2]

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
      XConfusionMat <- table(factor(oobPred, levels = Levels), factor(Levels[yy], levels = Levels))
      class_error <- (rowSums(XConfusionMat) - diag(XConfusionMat)) / (rowSums(XConfusionMat) + 1e-4)
      XConfusionMat <- cbind(XConfusionMat, class_error)
      ppForest$oobConfusionMat <- XConfusionMat
    } else {
      oobPred <- rowMeans(oobVotes, na.rm = TRUE)
      ppForest$oobErr <- mean((oobPred - yy)^2) # / mean((yy - mean(y))^2)
    }
  }

  class(ppForest) <- append(class(ppForest), "ODRF")
  # class(ppForest) <- "ODRF"
  return(ppForest)
}
