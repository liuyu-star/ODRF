#' Classification and Regression using the Ensemble of ODT-based Boosting Trees
#'
#' We use ODT as the basic tree model (base learner). To improve the performance of a boosting tree, we apply the feature bagging in this process, in the same
#' way as the random forest. Our final estimator is called the ensemble of ODT-based boosting trees, denoted by \code{ODBT}, is the average of many boosting trees.
#'
#' @param formula Object of class \code{formula} with a response describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} containing variables named in the formula. If \code{data} is missing it is obtained from the current environment by \code{formula}.
#' @param X An n by d numeric matrix (preferable) or data frame.
#' @param y A response vector of length n.
#' @param Xnew An n by d numeric matrix (preferable) or data frame containing predictors for the new data.
#' @param type Use \code{ODBT} for classification ("class") or regression ("reg").'auto' (default): If the response in \code{data} or \code{y} is a factor, "class" is used, otherwise regression is assumed.
#' @param model The basic tree model for boosting. We offer three options: "ODT" (default), "rpart," and "rpart.cpp" (improved "rpart")..
#' @param TreeRotate If or not to rotate the training data with the rotation matrix estimated by logistic regression before building the tree (default TRUE).
#' @param max.terms The maximum number of iterations for boosting trees.
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictors in the split node.
#' including \itemize{
#' \item{"RotMatPPO": projection pursuit optimization model (\code{\link{PPO}}), see \code{\link{RotMatPPO}} (default, model="PPR").}
#' \item{"RotMatRF": single feature similar to Random Forest, see \code{\link{RotMatRF}}.}
#' \item{"RotMatRand": random rotation, see \code{\link{RotMatRand}}.}
#' \item{"RotMatMake": users can define this function, for details see \code{\link{RotMatMake}}.}
#' }
#' @param FunDir The path to the \code{function} of the user-defined \code{NodeRotateFun} (default current working directory).
#' @param paramList List of parameters used by the functions \code{NodeRotateFun}. If left unchanged, default values will be used, for details see \code{\link{defaults}}.
#' @param ntrees The number of trees in the forest (default 100).
#' @param storeOOB If TRUE then the samples omitted during the creation of a tree are stored as part of the tree (default TRUE).
#' @param replacement if TRUE then n samples are chosen, with replacement, from training data (default TRUE).
#' @param stratify If TRUE then class sample proportions are maintained during the random sampling. Ignored if replacement = FALSE (default TRUE).
#' @param ratOOB  Ratio of 'out-of-bag' (default 1/3).
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
#' @param ... Optional parameters to be passed to the low level function.
#'
#' @return An object of class ODBT Containing a list components:
#' \itemize{
#' \item{\code{call}: The original call to ODBT.}
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
#' \item{\code{results}: The prediction results for new data \code{Xnew} using \code{ODBT}.}
#' }
#'
#' @seealso \code{\link{ODT}}
#'
#' @author Yu Liu and Yingcun Xia
#' @references Zhan, H., Liu, Y., & Xia, Y. (2024). Consistency of Oblique Decision Tree and its Boosting and Random Forest. arXiv preprint arXiv:2211.12653.
#' @references Tomita, T. M., Browne, J., Shen, C., Chung, J., Patsolic, J. L., Falk, B., ... & Vogelstein, J. T. (2020). Sparse projection oblique randomer forests. Journal of machine learning research, 21(104).
#' @keywords forest
#'
#' @examples
#' # Classification with Oblique Decision Tree.
#' data(seeds)
#' set.seed(221212)
#' train <- sample(1:209, 100)
#' train_data <- data.frame(seeds[train, ])
#' test_data <- data.frame(seeds[-train, ])
#' forest <- ODBT(varieties_of_wheat ~ ., train_data, test_data[, -8],model="rpart",
#' type = "class", parallel = FALSE, NodeRotateFun = "RotMatRF")
#' pred <- forest$results$prediction
#' # classification error
#' (mean(pred != test_data[, 8]))
#' forest <- ODBT(varieties_of_wheat ~ ., train_data, test_data[, -8],model="rpart.cpp",
#' type = "class", parallel = FALSE, NodeRotateFun = "RotMatRF")
#' pred <- forest$results$prediction
#' # classification error
#' (mean(pred != test_data[, 8]))
#'
#' # Regression with Oblique Decision Randome Forest.
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 80)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#' # To use ODT as the basic tree model for boosting, you need to set
#' the parameters model = "ODT" and NodeRotateFun = "RotMatPPO".
#' forest <- ODBT(Density ~ ., train_data, test_data[, -1],
#'   type = "reg",parallel = FALSE, model="ODT",
#'   NodeRotateFun = "RotMatPPO")
#' pred <- forest$results$prediction
#' # estimation error
#' mean((pred - test_data[, 1])^2)
#' forest <- ODBT(Density ~ ., train_data, test_data[, -1],
#'   type = "reg", parallel = FALSE,model="rpart.cpp",
#'   NodeRotateFun = "RotMatRF")
#' pred <- forest$results$prediction
#' # estimation error
#' mean((pred - test_data[, 1])^2)

#' @export
ODBT <- function(X, ...) {
  UseMethod("ODBT")
}


#' @rdname ODBT
#' @method ODBT formula
#' @aliases ODBT.formula
#' @export
ODBT.formula <- function(formula, data = NULL, Xnew = NULL, type = "auto",model=c("ODT","rpart","rpart.cpp")[1],TreeRotate=TRUE, max.terms=30,NodeRotateFun = "RotMatRF", FunDir = getwd(), paramList=NULL, #= list(numProj=ceiling(ifelse(is.null(data),ncol(eval(formula[[3]])),nrow(data))/2)),
                        ntrees = 100, storeOOB = TRUE, replacement = TRUE, stratify = TRUE, ratOOB = 0.368, parallel = TRUE,
                        numCores = Inf, MaxDepth = Inf, numNode = Inf, MinLeaf = ceiling(sqrt(ifelse(replacement,1,1-ratOOB)*ifelse(is.null(data),length(eval(formula[[2]])),nrow(data)))/3),
                        subset = NULL, weights = NULL,na.action = na.fail, catLabel = NULL, Xcat = 0, Xscale = "No", ...) {
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

  forest <- ODBT.compute(
    formula, Call, varName, X, y, Xnew,type,model,TreeRotate,
    max.terms, NodeRotateFun, FunDir, paramList,
    ntrees, storeOOB, replacement, stratify, ratOOB, parallel,
    numCores, MaxDepth, numNode, MinLeaf, subset, weights,
    na.action, catLabel, Xcat, Xscale
  )

  # class(forest) = append(class(forest),"ODBT.formula")
  return(forest)
}


#' @rdname ODBT
#' @method ODBT default
#' @aliases ODBT.default
#' @export
ODBT.default <- function(X, y, Xnew = NULL,
                         type = "auto",model=c("ODT","rpart","rpart.cpp")[1],TreeRotate=TRUE, max.terms=30,NodeRotateFun = "RotMatRF", FunDir = getwd(), paramList=NULL,
                         #= list(numProj=ceiling(ifelse(is.null(data),ncol(eval(formula[[3]])),nrow(data))/2)),
                         ntrees = 100, storeOOB = TRUE, replacement = TRUE, stratify = TRUE, ratOOB = 0.368, parallel = TRUE,
                         numCores = Inf, MaxDepth = Inf, numNode = Inf, MinLeaf = ceiling(sqrt(ifelse(replacement,1,1-ratOOB)*length(y))/3),
                         subset = NULL, weights = NULL,na.action = na.fail, catLabel = NULL, Xcat = 0, Xscale = "No", ...) {
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

  ODBT.compute(
    formula, Call, varName, X, y, Xnew,type,model,TreeRotate,
    max.terms, NodeRotateFun, FunDir, paramList,
    ntrees, storeOOB, replacement, stratify, ratOOB, parallel,
    numCores, MaxDepth, numNode, MinLeaf, subset, weights,
    na.action, catLabel, Xcat, Xscale
  )
}

# @importFrom RcppArmadillo fastLm
# @import fastmatrix ols.fit
#' @useDynLib ODRF, .registration = TRUE
#' @import Rcpp
#' @import doParallel
#' @import foreach
#' @import nnet
#' @importFrom parallel detectCores makeCluster clusterSplit stopCluster
#' @importFrom stats model.frame model.extract model.matrix na.fail
#' @importFrom rpart rpart rpart.control
#' @keywords internal
#' @noRd
ODBT.compute <- function(formula, Call, varName, X, y, Xnew,type,model,TreeRotate,
                         max.terms, NodeRotateFun, FunDir, paramList,
                         ntrees, storeOOB, replacement, stratify, ratOOB, parallel,
                        numCores, MaxDepth, numNode, MinLeaf, subset, weights,
                        na.action, catLabel, Xcat, Xscale) {
  #if (ntrees == 1) {
  #  stop("argument 'ntrees' must exceed 1")
  #}
  if (is.factor(y) && (type == "auto")) {
    type <- "class"
    warning("You are creating a forest for classification")
  }
  if (is.numeric(y) && (type == "auto")) {
    type <- "reg"
    warning("You are creating a forest for regression")
  }
  if (is.factor(y) && (type == "reg")) {
    stop(paste0("When ", formula[[2]], " is a factor type, 'type' cannot take 'regression'."))
  }
  # if (MinLeaf == 5) {
  #  MinLeaf <- ifelse(type == "mse", 5, 1)
  # }
  if ((ratOOB <= 0) || !storeOOB) {
    stop("out-of-bag indices for each tree are not stored. ODRF must be called with storeOOB = TRUE.")
  }

  n <- length(y)
  p <- ncol(X)
  yname <- NULL
  if (length(varName) > p) {
    yname <- varName[1]
    varName <- varName[-1]
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
  if (!is.numeric(X)){
    X=apply(X, 2, as.numeric)
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


  Levels=NULL; numClass <- 1
  #if (type %in% c("gini","entropy")) {
  if (type == "class") {
    y <- as.factor(y)
    Levels <- levels(y)
    y <- as.integer(y)

    if (length(Levels) == 1) {
      stop("the number of factor levels of response variable must be greater than one")
    }

    numClass <- length(Levels)
    classCt <- cumsum(table(y))
    if (stratify) {
      Cindex <- vector("list", numClass)
      for (m in 1L:numClass) {
        Cindex[[m]] <- which(y == m)
      }
    }
  }

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


  numCat <- 0
  n1=nrow(Xnew)
  if (sum(Xcat) > 0) {
    xj <- 1
    Xnew1 <- matrix(0, nrow = n1, ncol = length(unlist(catLabel))) # initialize training data matrix X
    # one-of-K encode each categorical feature and store in X
    for (j in seq_along(Xcat)) {
      catMap <- which(catLabel[[j]] %in% unique(Xnew[, Xcat[j]]))
      indC <- catLabel[[j]][catMap]
      Xnewj <- (matrix(Xnew[, Xcat[j]], n1, length(indC)) == matrix(indC, n1, length(indC), byrow = TRUE)) + 0

      if (length(indC) > length(catLabel[[j]])) {
        Xnewj <- Xnewj[, seq_along(catLabel[[j]])]
      }

      xj1 <- xj + length(catLabel[[j]])
      Xnew1[, (xj:(xj1 - 1))[catMap]] <- Xnewj
      xj <- xj1
    }

    Xnew <- cbind(Xnew1, Xnew[, -Xcat])
    #p <- ncol(Xnew)
    numCat <- length(unlist(catLabel))
    rm(Xnew1)
    rm(Xnewj)
  }
  if (!is.numeric(Xnew)){
    Xnew=apply(Xnew, 2, as.numeric)
  }

  # Variable scaling.
  if (Xscale != "No") {
    indp <- (sum(numCat) + 1):p
    Xnew[, indp] <- (Xnew[, indp] - matrix(minCol, n1, length(indp), byrow = T)) /
      matrix(maxminCol, n1, length(indp), byrow = T)
  }


  if(is.null(paramList$numProj)){paramList$numProj=ceiling(p/2)}
  #paramList <- defaults(paramList, split="mse", p, weights, catLabel)
  ppForest <- list(
    call = Call, terms = Terms, type = type, Levels = Levels, NodeRotateFun = FALSE,
    predicted=NULL, paramList = paramList, oobErr = NULL, oobConfusionMat = NULL
  )
  ppForest$data <- list(
    subset = subset, weights = weights, na.action = na.action, n = n, p = p, varName = varName,
    Xscale = Xscale, minCol = minCol, maxminCol = maxminCol, Xcat = Xcat, catLabel = catLabel,
    TreeRotate=TreeRotate
  )
  ppForest$tree <- list(lambda = 0, FunDir = FunDir, MaxDepth = MaxDepth, MinLeaf = MinLeaf, numNode = numNode)
  ppForest$forest <- list(
    ntrees = ntrees, ratOOB = ratOOB, storeOOB = storeOOB, replacement = replacement, stratify = stratify,
    parallel = parallel, numCores = numCores
  )

  seqn=seq(n)
  index <- function(...) {
    TDindx <- seqn
    if (replacement) {
      go <- TRUE
      while (go) {
        # make sure each class is represented in proportion to classes in initial dataset
        if (stratify && (type == "class")) {
          if (classCt[1L] != 0L) {
            TDindx[1:classCt[1L]] <- sample(Cindex[[1L]], classCt[1L], replace = TRUE)
          }
          for (z in 2:numClass) {
            if (classCt[z - 1L] != classCt[z]) {
              TDindx[(classCt[z - 1L] + 1L):classCt[z]] <- sample(Cindex[[z]], classCt[z] - classCt[z - 1L], replace = TRUE)
            }
          }
        } else {
          TDindx <- sample(seqn, n, replace = TRUE)
        }
        go <- all(seqn %in% TDindx)
      }
    } else {
      TDindx <- sample(seqn, ceiling(n * (1 - ratOOB)), replace = FALSE)
    }

    return(TDindx)
  }

  if("type"%in%names(Call0)){Call0=Call0[-which("type"==names(Call0))]}
  #if("type"%in%names(Call0)){names(Call0)["type"==names(Call0)]="split"}
  #mtry=ifelse(is.null(paramList$numProj),ceiling(p/2),paramList$numProj)
  mtry=ifelse(ntrees==1,p,paramList$numProj)
  if(TreeRotate)varName=c(varName,"XB")

  runTree <- function(itree, ...) {
    #set.seed(seed + itree)
    #options (warn = -1)
    #ow <- options("warn")
    #options(warn = -1)
    #warnings('off')

    AIC=Inf;nterm=1;nterm.fails=0
    if(type=="class"&&(length(Levels)>2)){
      #nn=length(TDindx)
      #Y=(matrix(y, n, numClass) == matrix(seq(numClass), n, numClass, byrow = TRUE)) + 0
      YRes=Y
      #Ftk=matrix(table(y[TDindx])/nn, nn, numClass, byrow = TRUE)
      #PFtk=matrix(0, length(NTD), numClass)

      COEF=matrix(0,max.terms+1,numClass)
      RES=matrix(0,n,numClass)
      fitted=rep(list(matrix(0, n, max.terms)),numClass)# vector("list", numClass)
      pred=rep(list(matrix(0, n1, max.terms)),numClass)#vector("list", numClass)
      for (t in seq(max.terms)) {
        if (ntrees == 1) {
          TDindx=seqn
          J = 1:p
        }else{
          TDindx=index()
          J = sample(1:p, mtry)
        }
        NTD <- setdiff(seqn, TDindx)


        XB = X
        XBtest=Xnew
        if(model=="rpart"){
          XB = X[,J]
          XBtest=Xnew[,J]
        }
        if(TreeRotate){
          B <- nnet(XB[TDindx, ],  YRes[TDindx,], size=1, trace=FALSE)$wts[2:(1 + ncol(XB))]#linout = TRUE,MaxNWts= p+5
          XB = data.frame(XB, XB=XB%*%B)
          XBtest = data.frame(XBtest, XB=XBtest%*%B)
        }
        XB = data.frame(XB)
        XBtest = data.frame(XBtest)

        #XB = data.frame(X,y=YRes[,k])[TDindx,J]
        #XB = data.frame(cbind(X[, J], X[, I]%*%B))
        #X1B = data.frame(cbind(X1[, I], X1[, I]%*%B))
        #sse=sse.oob=0
        #ppForestT[[t]]=vector("list", numClass)
        #prSum=rowSums(exp(Ftk))
        for (k in seq(numClass)) {
          #pr=exp(Ftk[,k])/prSum
          if(model=="ODT"){
            Tree <- ODT_compute(formula, Call0, varName,
                                X =XB[TDindx, ], y = YRes[TDindx,k], split="mse",lambda=0, NodeRotateFun=NodeRotateFun,FunDir=FunDir, paramList=paramList, MaxDepth=MaxDepth, numNode=numNode,
                                MinLeaf=MinLeaf, Levels=Levels, subset = NULL, weights = weights[TDindx], na.action = NULL, catLabel=catLabel, Xcat = 0L, Xscale = "No",TreeRandRotate=FALSE)
          }
          if(model=="rpart"){
            Tree <-rpart(y~., data.frame(XB,y=YRes[,k])[TDindx,],control = rpart.control(minbucket = MinLeaf))
          }
          fitted[[k]][,nterm]=predict(Tree, XB)
          pred[[k]][,nterm]=predict(Tree, XBtest)

          # Tree[["predicted"]])
          LM = lm(y~.,data.frame(y=Y[,k], fitted[[k]][,seq(nterm)]))
          #LM = fastLm(y~.,data.frame(y=Y[,k], fitted[[k]]))
          #LM <- ols.fit(x = cbind(1,fitted[[k]]), y = Y[,k])
          RES[,k] = LM$residuals
          #sse = sse + sum(LM$residuals^2)
          #sse.oob = sse.oob + sum(LM$residuals[NTD]^2)
          LM$coefficients[is.na(LM$coefficients)] <- 0.0
          COEF[seq(nterm+1),k] = LM$coefficients
        }

        #J = setdiff(1:n, J)
        aic = log((sum(RES^2) + sum(RES[NTD,]^2))/((n+length(NTD))*numClass)) + log(p)*nterm*log(n)/n
        if (aic < AIC)
        {
          YRes = RES
          #pred = cbind(pred,predict(Tree, Xnew))
          #pred[,count.terms] = predict(fit.B, data.frame(cbind(X1[,I], X1[,I]%*%B)))
          #pred[,count.terms] = predict(fit.B, list(Xk = cbind(X1[,I], X1[,I]%*%B)))
          AIC = aic # AIC(A, k = log(n))
          #nterm=t
          COEF0 = COEF[seq(1+nterm),]
          nterm=nterm+1
          #count.terms0 = count.terms
          #ppForestT[[t]]=c(list(rotdims=Tree[["data"]][["rotdims"]],rotmat=Tree[["data"]][["rotmat"]]),Tree$structure)
          #count.terms = min(count.terms + 1,max.terms)
          nterm.fails=0
        }else{
          nterm.fails = nterm.fails + 1
          if (nterm.fails > 5){
            break
          }
        }
      }

      #predictions=matrix(0,n1,numClass)
      #for (k in seq(numClass)) {
      #nterm=min(nterm,max.terms)
      nterm=nterm-1
      #COEF0=COEF0*(1-is.na(COEF0))
      pred=vapply(seq(numClass), function(k){
        cbind(1,pred[[k]][,seq(nterm)]) %*% COEF0[,k]
      }, rep(0, n1))
      predictions=c(pred)

      fitted=vapply(seq(numClass), function(k){
        cbind(1,fitted[[k]][,seq(nterm)]) %*% COEF0[,k]
      }, rep(0, n))
      sse = mean((fitted-Y)^2)
      predictions = c(predictions,sse)
      #predictions=Levels[max.col(predictions)]

    }else{

      #if((type=="regression")||(length(Levels)==2)){
      if(length(Levels)==2){
        yy=y-1
      }else{
        yy=y
      }
      yres=yy

      fitted=matrix(0, n, max.terms)
      pred=matrix(0, n1, max.terms)
      for (t in seq(max.terms)) {
        if (ntrees == 1) {
          TDindx=seqn
          J = 1:p
        }else{
          TDindx=index()
          J = sample(1:p, mtry)
        }
        NTD <- setdiff(seqn, TDindx)

        XB = X
        XBtest=Xnew
        if(model=="rpart"){
          XB = X[,J]
          XBtest=Xnew[,J]
        }
        if(TreeRotate){
          B <- nnet(XB[TDindx, ], yres[TDindx], size=1, trace=FALSE)$wts[2:(1 + ncol(XB))]#linout = TRUE,MaxNWts= p+5
          XB = data.frame(XB, XB=XB%*%B)
          XBtest = data.frame(XBtest, XB=XBtest%*%B)
        }
        XB = data.frame(XB)
        XBtest = data.frame(XBtest)

        if(model=="ODT"){
          Tree <- ODT_compute(formula, Call0, varName,
                              X = XB[TDindx, ], y = yres[TDindx], split="mse",lambda=0, NodeRotateFun=NodeRotateFun,FunDir=FunDir, paramList=paramList, MaxDepth=MaxDepth, numNode=numNode,
                              MinLeaf=MinLeaf, Levels=Levels, subset = NULL, weights = weights[TDindx], na.action = NULL, catLabel=catLabel, Xcat = 0L, Xscale = "No",TreeRandRotate=FALSE)
        }
        if(model=="rpart"){
          Tree <-rpart(y~., data.frame(XB,y=yres)[TDindx,],control = rpart.control(minbucket = MinLeaf))
        }
        fitted[,nterm] =predict(Tree, XB)
        pred[,nterm] = predict(Tree, XBtest)

        LM = lm(y~.,data.frame(y=yy, fitted[,seq(nterm)]))
        #LM = fastLm(y~.,data.frame(y=yy, fitted))#,silent = TRUE)
        #LM <- ols.fit(x = cbind(1,fitted), y = yy)
        res.t = LM$residuals

        #J = setdiff(1:n, J)
        aic = log((sum(res.t^2)+sum(res.t[NTD]^2))/(n+length(NTD)))  + log(p)*nterm*log(n)/n
        if (aic < AIC)
        {
          yres = res.t
          #pred[,count.terms] = predict(fit.B, data.frame(cbind(X1[,I], X1[,I]%*%B)))
          #pred[,count.terms] = predict(fit.B, list(Xk = cbind(X1[,I], X1[,I]%*%B)))
          coef = LM$coefficients
          AIC = aic # AIC(A, k = log(n))
          #nterm=t
          #count.terms0 = count.terms
          #ppForestT[[t]]=c(list(rotdims=Tree[["data"]][["rotdims"]],rotmat=Tree[["data"]][["rotmat"]]),Tree$structure)
          nterm = nterm+1
          nterm.fails = 0
        }else{
          nterm.fails = nterm.fails + 1
          if (nterm.fails > 5){
            break
          }
        }
      }

      nterm=nterm-1
      coef[is.na(coef)] <- 0.0

      predictions = cbind(1,pred[,seq(nterm)]) %*% coef#predict(LM0, data.frame(pred))
      sse = mean((cbind(1,fitted[,seq(nterm)]) %*% coef-yy)^2)
      predictions = c(predictions,sse)

      #if(length(Levels)==2){
      #  predictions = Levels[(predictions>0.5)+1]
      #}
      #coef = LM0$coefficients
    }
    #}

    #options(ow)
    #warnings('on')
    #options (warn = 0)
    return(predictions)
    #return(c(ppForestT, list(oobErr = oobErr0, oobIndex = NTD, oobPred = pred, ts=nterms)))
  }

  #op=options(nwarnings)
  #op=options(nwarnings = 1)
  #VALUE <- rep(ifelse(type == "classification","0",0), n1)
  #VALUE <- rep(ifelse(type=="classification"&&(length(Levels)>2),"0",0), n1)
  VALUE <- rep(0, ifelse(type=="class"&&(length(Levels)>2),n1*numClass,n1)+1)
  if(type=="class"&&(length(Levels)>2)){
    Y=(matrix(y, n, numClass) == matrix(seq(numClass), n, numClass, byrow = TRUE)) + 0
  }else{
    Y=y
  }
  Y=as.matrix(Y)
  #nnet=nnet::nnet.default;rpart=rpart::rpart;
  #control=rpart::rpart.control;predict=rpart:::predict.rpart
  if (parallel&&(ntrees>1)) {
    # RNGkind("L'Ecuyer-CMRG")
    if (is.infinite(numCores)) {
      # Use all but 1 core if numCores=0.
      numCores <- parallel::detectCores() - 1L # logical = FALSE
    }
    numCores <- min(numCores, ntrees)
    gc()

    # cl <- parallel::makePSOCKcluster(num.cores)
    # library("ODBT1")
    # library(foreach)
    # foreach::registerDoSEQ()
    cl <- parallel::makeCluster(numCores, type = ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK"))
    chunks <- parallel::clusterSplit(cl, seq_len(ntrees))
    doParallel::registerDoParallel(cl, numCores)
    # set.seed(seed)
    icore <- NULL
    Votes <- foreach::foreach(
      icore = seq_along(chunks), .combine = 'cbind', .export = c("ODT.compute"),
      .packages = c("ODRF","nnet","rpart"), .noexport = "ppForest"
    ) %dopar% {
      #lapply(chunks[[icore]], runTree)
      vapply(chunks[[icore]], function(t){
        if(model=="rpart.cpp"){
          GBDTCpp(X,y,Xnew,Y,#nnet,rpart,control,predict,
                  numClass, maxTerms=max.terms, ntrees, mtry, MinLeaf,replacement,ratOOB)
        }else{
          runTree()
        }
        }, VALUE)
    }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)

    # do.call(rbind.fill,list1)
    #Votes <- t(do.call("c", ppForestT))
    # ppForest$structure=NULL
    # for (i in 1:numCores) {
    #  ppForest$structure=c(ppForest$structure,ppForestT[[i]])
    # }
  } else {
    # Use just one core.
    #Votes <- vapply(1:ntrees, runTree, VALUE)
    #Votes <- vapply(1:ntrees, function(t){
     # GBDTCpp(X,as.numeric(y),Xnew,Y,nnet,rpart,control,predict,
    #          numClass, maxTerms=max.terms, ntrees, mtry, ratOOB, MinLeaf)}
    #  , VALUE)
      #if(model=="rpart.cpp"){
      #  Votes <- ODBTCpp(X,y,Xnew,numClass, maxTerms=max.terms, ntrees, mtry, MinLeaf,replacement,ratOOB)
      #}else{
      #  Votes <- vapply(1:ntrees, runTree, VALUE)
      #}
    #Votes <- ODBTCpp(X,y,Xnew,numClass, maxTerms=max.terms, ntrees, mtry, ratOOB, MinLeaf)
    Votes <- vapply(1:ntrees, function(t){
      if(model=="rpart.cpp"){
        BODTCpp(X,y,Xnew,Y,#nnet,rpart,control,predict,
                numClass, maxTerms=max.terms, ntrees, mtry, MinLeaf,replacement,ratOOB)
      }else{
        runTree()
      }
    }, VALUE)
  }
  #options(op)

  ##############################################################################
  Votes=as.matrix(Votes)
  #weights <- rep(1, ntrees)
  weights=1/(Votes[length(VALUE),]+1)
  weights <- weights / sum(weights)
  Votes=Votes[-length(VALUE),,drop = FALSE]
  #if(ntrees==1){
  #  prob=1
  #  pred=Votes
  #}else{
  if (type=="class"&&(length(Levels)>2)) {
    if(1==2){
      weights <- rep(weights, n1)
      Votes <- factor(c(Votes), levels = Levels)
      Votes <- as.integer(Votes) + numClass * rep(0:(n1 - 1), rep(ntrees, n1))
      Votes <- aggregate(c(rep(0, n1 * numClass), weights), by = list(c(1:(n1 * numClass), Votes)), sum)[, 2]

      prob <- matrix(Votes, n1, numClass, byrow = TRUE)
    }
    prob <- Votes %*% weights
    prob=matrix(prob, n1, numClass)
    #prob <- prob / matrix(rowSums(prob), n1, numClass)
    prob=exp(prob)/rowSums(exp(prob))
    colnames(prob) <- Levels
    # pred=apply(prob,1,which.max);
    pred <- max.col(prob) ## "random"
    pred <- Levels[pred]
  } else {
    prob <- weights #/ sum(weights)
    pred <- Votes %*% prob
    if(length(Levels)==2){
      pred = Levels[(pred>0.5)+1]
      #prob <- cbind(1-prob,prob)
      #prob <- prob / matrix(rowSums(prob), n1, numClass)
      #colnames(prob) <- Levels
      #pred <- Levels[max.col(prob)]
    }
    # pred=colMeans(Votes);
    # prob=NULL
  }
  #}


  ppForest$results=list(probability=prob, prediction=c(pred))

  class(ppForest) <- append(class(ppForest), "ODBT")
  # class(ppForest) <- "ODBT"
  return(ppForest)
}

