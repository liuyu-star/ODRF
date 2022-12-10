#' Create a Random Matrix: Binary
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param sparsity a real number in \eqn{(0,1)} that specifies the distribution of non-zero elements in the random matrix.
#' @param prob a probability \eqn{\in (0,1)} used for sampling from
#' \eqn{{-1,1}} where \code{prob = 0} will only sample -1 and \code{prob = 1} will only sample 1.
#' @param catMap a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 3
#' sparsity <- 0.25
#' prob <- 0.5
#' set.seed(4)
#' (a <- RandMatBinary(p, d, sparsity, prob))
RandMatBinary <- function(p, d, sparsity, prob, catMap = NULL,...) {
  nnzs <- round(p * d * sparsity)
  ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
  
  ## Determine if categorical variables need to be taken
  ## into consideration
  if (is.null(catMap)) {
    randomMatrix <- cbind(((ind - 1L)%%p) + 1L, floor((ind -
                                                         1L)/p) + 1L, sample(c(1L, -1L), nnzs, replace = TRUE,
                                                                             prob = c(prob, 1 - prob)))
  } else {
    rw <- ((ind - 1L)%%p) + 1L
    for (j in 1:length(catMap)) {
      isj <- rw %in% catMap[[j]]
      rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
    }
    randomMatrix <- cbind(rw, floor((ind - 1L)/p) + 1L, sample(c(1L,
                                                                 -1L), nnzs, replace = TRUE, prob = c(prob, 1 - prob)),
                          deparse.level = 0)
  }
  return(randomMatrix)
}


#' Create a Random Matrix: Continuous
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param sparsity a real number in \eqn{(0,1)} that specifies the distribution of non-zero elements in the random matrix.
#' @param catMap a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom RcppZiggurat zrnorm
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 3
#' sparsity <- 0.25
#' set.seed(4)
#' (a <- RandMatContinuous(p, d, sparsity))
RandMatContinuous <- function(p, d, sparsity, catMap = NULL,
                              ...) {
  nnzs <- round(p * d * sparsity)
  ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
  
  if (is.null(catMap)) {
    randomMatrix <- cbind(((ind - 1L)%%p) + 1L, floor((ind -
                                                         1L)/p) + 1L, zrnorm(nnzs))
  } else {
    rw <- ((ind - 1L)%%p) + 1L
    for (j in 1:length(catMap)) {
      isj <- rw %in% catMap[[j]]
      rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
    }
    randomMatrix <- cbind(rw, floor((ind - 1L)/p) + 1L, zrnorm(nnzs),
                          deparse.level = 0)
  }
  return(randomMatrix)
}


#' Create a Random Matrix: Random Forest (RF)
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param catMap a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 3
#' paramList <- list(p = p, d = d)
#' set.seed(4)
#' (a <- do.call(RandMatRF, paramList))
RandMatRF <- function(p, d, catMap = NULL, ...) {
  if (d > p) {
    stop("ERROR: parameter d is greater than the number of dimensions p.")
  }
  
  if (is.null(catMap)) {
    randomMatrix <- cbind(sample.int(p, d, replace = FALSE),
                          1:d, rep(1L, d))
  } else {
    pnum <- catMap[[1L]][1L] - 1L
    rw <- sample.int(p, d, replace = FALSE)
    isCat <- rw > pnum
    for (j in (pnum + 1L):p) {
      isj <- rw == j
      rw[isj] <- sample(catMap[[j - pnum]], length(rw[isj]),
                        replace = TRUE)
    }
    randomMatrix <- cbind(rw, 1:d, rep(1L, d), deparse.level = 0)
  }
  return(randomMatrix)
}


#' Create a Random Matrix: Poisson
#'
#' Samples a binary projection matrix where sparsity is distributed
#' \eqn{Poisson(\lambda)}.
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param lambda passed to the \code{\link[stats]{rpois}} function for generation of non-zero elements in the random matrix.
#' @param catMap a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom stats rpois
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 8
#' lambda <- 0.5
#' paramList <- list(p = p, d = d, lambda = lambda)
#' set.seed(8)
#' (a <- do.call(RandMatPoisson, paramList))
RandMatPoisson <- function(p, d, lambda, catMap = NULL, ...) {
  if (lambda <= 0) {
    stop("ERROR: Wrong parameter for Poisson, make sure lambda > 0.")
  }
  
  nnzPerCol <- stats::rpois(d, lambda)
  while (!any(nnzPerCol)) {
    nnzPerCol <- stats::rpois(d, lambda)
  }
  
  nnzPerCol[nnzPerCol > p] <- p
  nnz <- sum(nnzPerCol)
  nz.rows <- integer(nnz)
  nz.cols <- integer(nnz)
  start.idx <- 1L
  for (i in seq.int(d)) {
    if (nnzPerCol[i] != 0L) {
      end.idx <- start.idx + nnzPerCol[i] - 1L
      nz.rows[start.idx:end.idx] <- sample.int(p, nnzPerCol[i],
                                               replace = FALSE)
      nz.cols[start.idx:end.idx] <- i
      start.idx <- end.idx + 1L
    }
  }
  
  if (is.null(catMap)) {
    randomMatrix <- cbind(nz.rows, nz.cols, sample(c(-1L,
                                                     1L), nnz, replace = TRUE))
  } else {
    pnum <- catMap[[1L]][1L] - 1L
    isCat <- nz.rows > pnum
    for (j in (pnum + 1L):p) {
      isj <- nz.rows == j
      nz.rows[isj] <- sample(catMap[[j - pnum]], length(nz.rows[isj]),
                             replace = TRUE)
    }
    randomMatrix <- cbind(nz.rows, nz.cols, sample(c(-1L,
                                                     1L), nnz, replace = TRUE), deparse.level = 0)
  }
  return(randomMatrix)
}


#' Create a Random Matrix: FRC
#'
#'
#' @param p integer the number of dimensions.
#' @param d integer the number of desired columns in the projection matrix.
#' @param nmix integer mupliplier to \code{d} to specify the number of non-zeros.
#' @param catMap a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 2
#' nmix <- 5
#' paramList <- list(p = p, d = d, nmix = nmix)
#' set.seed(4)
#' (a <- do.call(RandMatFRC, paramList))
RandMatFRC <- function(p, d, nmix, catMap = NULL, ...) {
  if (nmix > p) {
    stop("ERROR: parameter nmix is greater than the number of dimensions p.")
  }
  
  nnz <- nmix * d
  nz.rows <- integer(nnz)
  nz.cols <- integer(nnz)
  start.idx <- 1L
  for (i in seq.int(d)) {
    end.idx <- start.idx + nmix - 1L
    nz.rows[start.idx:end.idx] <- sample.int(p, nmix, replace = FALSE)
    nz.cols[start.idx:end.idx] <- i
    start.idx <- end.idx + 1L
  }
  
  if (is.null(catMap)) {
    randomMatrix <- cbind(nz.rows, nz.cols, runif(nnz, -1,
                                                  1))
  } else {
    pnum <- catMap[[1L]][1L] - 1L
    isCat <- nz.rows > pnum
    for (j in (pnum + 1L):p) {
      isj <- nz.rows == j
      nz.rows[isj] <- sample(catMap[[j - pnum]], length(nz.rows[isj]),
                             replace = TRUE)
    }
    randomMatrix <- cbind(nz.rows, nz.cols, runif(nnz, -1,
                                                  1), deparse.level = 0)
  }
  return(randomMatrix)
}


#' Create a Random Matrix: FRCN
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param nmix mupliplier to \code{d} to specify the number of non-zeros.
#' @param catMap a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom RcppZiggurat zrnorm
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 8
#' nmix <- 5
#' paramList <- list(p = p, d = d, nmix = nmix)
#' set.seed(8)
#' (a <- do.call(RandMatFRCN, paramList))
RandMatFRCN <- function(p, d, nmix, catMap = NULL, ...) {
  if (d > p) {
    stop("ERROR: parameter d is greater than the number of dimensions p.")
  }
  
  nnz <- nmix * d
  nz.rows <- integer(nnz)
  nz.cols <- integer(nnz)
  start.idx <- 1L
  for (i in seq.int(d)) {
    end.idx <- start.idx + nmix - 1L
    nz.rows[start.idx:end.idx] <- sample.int(p, nmix, replace = FALSE)
    nz.cols[start.idx:end.idx] <- i
    start.idx <- end.idx + 1L
  }
  
  if (is.null(catMap)) {
    randomMatrix <- cbind(nz.rows, nz.cols, zrnorm(nnz))
  } else {
    pnum <- catMap[[1L]][1L] - 1L
    isCat <- nz.rows > pnum
    for (j in (pnum + 1L):p) {
      isj <- nz.rows == j
      nz.rows[isj] <- sample(catMap[[j - pnum]], length(nz.rows[isj]),
                             replace = TRUE)
    }
    randomMatrix <- cbind(nz.rows, nz.cols, zrnorm(nnz),
                          deparse.level = 0)
  }
  return(randomMatrix)
}


#' Create a Random Matrix: ts-patch
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param pwMin the minimum patch size to sample.
#' @param pwMax the maximum patch size to sample.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 8
#' pwMin <- 3
#' pwMax <- 6
#' paramList <- list(p = p, d = d, pwMin = pwMin, pwMax = pwMax)
#' set.seed(8)
#' (a <- do.call(RandMatTSpatch, paramList))
RandMatTSpatch <- function(p, d, pwMin, pwMax, ...) {
  if (pwMin > pwMax) {
    stop("ERROR: parameter pwMin is greater than pwMax.")
  }
  
  # pw holds all sizes of patch to filter on.  There will
  # be d patches of varying sizes
  pw <- sample.int(pwMax - pwMin, d, replace = TRUE) + pwMin
  
  # nnz is sum over how many points the projection will sum
  # over
  nnz <- sum(pw)
  nz.rows <- integer(nnz)  # vector to hold row coordinates of patch points
  nz.cols <- integer(nnz)  # vector to hold column coordinates of patch points
  
  # Here we create the patches and store them
  start.idx <- 1L
  for (i in seq.int(d)) {
    pw.start <- sample.int(p, 1)  # Sample where to start the patch
    end.idx <- start.idx + pw[i] - 1L  # Set the ending point of the patch
    for (j in 1:pw[i]) {
      # Handle boundary cases where patch goes past end of
      # ts
      if (j + pw.start - 1L > p) {
        end.idx <- j + start.idx - 1L
        break
      }
      nz.rows[j + start.idx - 1L] <- pw.start + j - 1L
      nz.cols[j + start.idx - 1L] <- i
    }
    start.idx <- end.idx + 1L
  }
  random.matrix <- cbind(nz.rows, nz.cols, rep(1L, nnz))
  random.matrix <- random.matrix[random.matrix[, 1] > 0, ]  # Trim entries that are 0
}



#' Create a Random Matrix: image-patch
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param ih the height (px) of the image.
#' @param iw the width (px) of the image.
#' @param pwMin the minimum patch size to sample.
#' @param pwMax the maximum patch size to sample.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @export
#'
#' @examples
#'
#' p <- 28^2
#' d <- 8
#' ih <- iw <- 28
#' pwMin <- 3
#' pwMax <- 6
#' paramList <- list(p = p, d = d, ih = ih, iw = iw, pwMin = pwMin, pwMax = pwMax)
#' set.seed(8)
#' (a <- do.call(RandMatImagePatch, paramList))
RandMatImagePatch <- function(p, d, ih, iw, pwMin, pwMax, ...) {
  if (pwMin > pwMax) {
    stop("ERROR: parameter pwMin is greater than pwMax.")
  }
  
  pw <- sample.int(pwMax - pwMin + 1L, 2 * d, replace = TRUE) +
    pwMin - 1L
  sample.height <- ih - pw[1:d] + 1L
  sample.width <- iw - pw[(d + 1L):(2 * d)] + 1L
  nnz <- sum(pw[1:d] * pw[(d + 1L):(2 * d)])
  nz.rows <- integer(nnz)
  nz.cols <- integer(nnz)
  start.idx <- 1L
  for (i in seq.int(d)) {
    top.left <- sample.int(sample.height[i] * sample.width[i],
                           1L)
    top.left <- floor((top.left - 1L)/sample.height[i]) *
      (ih - sample.height[i]) + top.left
    # top.left <- floor((top.left - 1L)/sample.height[i]) +
    # top.left
    end.idx <- start.idx + pw[i] * pw[i + d] - 1L
    nz.rows[start.idx:end.idx] <- sapply((1:pw[i + d]) -
                                           1L, function(x) top.left:(top.left + pw[i] - 1L) +
                                           x * ih)
    nz.cols[start.idx:end.idx] <- i
    start.idx <- end.idx + 1L
  }
  # random.matrix <- cbind(nz.rows, nz.cols,
  # sample(c(-1L,1L), nnz, replace = TRUE))
  random.matrix <- cbind(nz.rows, nz.cols, rep(1L, nnz))
}


#' Create a Random Matrix: image-control
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param ih the height (px) of the image.
#' @param iw the width (px) of the image.
#' @param pwMin the minimum patch size to sample.
#' @param pwMax the maximum patch size to sample.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @export
#'
#' @examples
#'
#' p <- 28^2
#' d <- 8
#' ih <- iw <- 28
#' pwMin <- 3
#' pwMax <- 6
#' paramList <- list(p = p, d = d, ih = ih, iw = iw, pwMin = pwMin, pwMax = pwMax)
#' set.seed(8)
#' (a <- do.call(RandMatImageControl, paramList))
RandMatImageControl <- function(p, d, ih, iw, pwMin, pwMax, ...) {
  if (pwMin > pwMax) {
    stop("ERROR: parameter pwMin is greater than pwMax.")
  }
  
  pw <- sample.int(pwMax - pwMin + 1L, 2 * d, replace = TRUE) +
    pwMin - 1L
  nnzPerCol <- pw[1:d] * pw[(d + 1L):(2 * d)]
  sample.height <- ih - pw[1:d] + 1L
  sample.width <- iw - pw[(d + 1L):(2 * d)] + 1L
  nnz <- sum(nnzPerCol)
  nz.rows <- integer(nnz)
  nz.cols <- integer(nnz)
  start.idx <- 1L
  for (i in seq.int(d)) {
    end.idx <- start.idx + nnzPerCol[i] - 1L
    nz.rows[start.idx:end.idx] <- sample.int(p, nnzPerCol[i],
                                             replace = FALSE)
    nz.cols[start.idx:end.idx] <- i
    start.idx <- end.idx + 1L
  }
  # random.matrix <- cbind(nz.rows, nz.cols,
  # sample(c(-1L,1L), nnz, replace = TRUE))
  random.matrix <- cbind(nz.rows, nz.cols, rep(1L, nnz))
}



#' Create a Random Matrix: custom
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param nnzSample a vector specifying the number of non-zeros to
#' sample at each \code{d}.  Each entry should be less than \code{p}.
#' @param nnzProb a vector specifying probabilities in one-to-one correspondance
#' with \code{nnzSample}.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom RcppZiggurat zrnorm
#'
#' @export
#'
#' @examples
#'
#' p <- 28
#' d <- 8
#' nnzSample <- 1:8
#' nnzProb <- 1 / 36 * 1:8
#' paramList <- list(p = p, d = d, nnzSample, nnzProb)
#' set.seed(8)
#' (a <- do.call(RandMatCustom, paramList))
RandMatCustom <- function(p, d, nnzSample, nnzProb, ...) {
  try({
    if (any(nnzSample > p) | any(nnzSample == 0)) {
      stop("nnzs per projection must be no more than the number of features.")
    }
  })
  nnzPerCol <- sample(nnzSample, d, replace = TRUE, prob = nnzProb)
  nnz <- sum(nnzPerCol)
  nz.rows <- integer(nnz)
  nz.cols <- integer(nnz)
  start.idx <- 1L
  for (i in seq.int(d)) {
    end.idx <- start.idx + nnzPerCol[i] - 1L
    nz.rows[start.idx:end.idx] <- sample.int(p, nnzPerCol[i],
                                             replace = FALSE)
    nz.cols[start.idx:end.idx] <- i
    start.idx <- end.idx + 1L
  }
  random.matrix <- cbind(nz.rows, nz.cols, zrnorm(nnz))
}



#' Default values passed to RandMat*
#'
#' Given the parameter list and the categorical map this function
#' populates the values of the parameter list accoding to our 'best'
#' known general use case parameters.
#'
#' @param ncolX an integer denoting the number of columns in the design
#' matrix X.
#' @param paramList a list (possibly empty), to be populated with a set
#' of default values to be passed to a RandMat* function.
#' @param cat.map a list specifying which columns in X correspond to the
#' same one-of-K encoded feature. Each element of cat.map is a numeric
#' vector specifying the K column indices of X corresponding to the same
#' categorical feature after one-of-K encoding. All one-of-K encoded
#' features in X must come after the numeric features. The K encoded
#' columns corresponding to the same categorical feature must be placed
#' contiguously within X. The reason for specifying cat.map is to adjust
#' for the fact that one-of-K encoding cateogorical features results in
#' a dilution of numeric features, since a single categorical feature is
#' expanded to K binary features. If cat.map = NULL, then RerF assumes
#' all features are numeric (i.e. none of the features have been
#' one-of-K encoded).
#'
#' @return If \code{cat.map} is NULL, then
#' \itemize{
#' \item \code{p} is set to the number of columns of \code{X}
#' \item \code{d} is set to the ceiling of the square root of the number of columns of \code{X}
#' \item \code{sparsity}: if \eqn{\code{ncol(X)} \ge 10}, then sparsity is set
#' to 3 / \code{ncol{X}}, otherwise it is set to 1 / \code{ncol(X)}.
#' \item \code{prob} defaults to 0.5.
#' }
#'
#' @keywords internal
#'
defaults <- function(paramList,p,catLabel) {
  if (is.null(paramList[["p"]])) {
    paramList$p <- p
  }
  if (is.null(paramList[["sparsity"]])) {
    paramList$sparsity <- ifelse(paramList$p >=10, 3/paramList$p, 1/paramList$p)
  }
  if (is.null(paramList[["prob"]])) {
    paramList$prob <- 0.5
  }
  if (is.null(paramList[["catLabel"]])) {
    paramList$catLabel = catLabel
  }
  if (is.null(paramList[["d"]])) {
    #  q<- min(ceiling(length(y)^0.4),ceiling(paramList$p*2/3))
    #  paramList$d <-min(max(5, ceiling(paramList$p/q)),paramList$p)
    paramList$d <- ceiling(sqrt(paramList[["p"]]))
  }

  return(paramList)
}

if(1==2){
  defaults <- function(paramList,x,y,method,catLabel) {
    paramList$x = x
    paramList$y = y
    paramList$p <- ncol(x)
    
    if (is.null(paramList[["sparsity"]])) {
      paramList$sparsity <- ifelse(paramList$p >=10, 3/paramList$p, 1/paramList$p)
    }
    if (is.null(paramList[["prob"]])) {
      paramList$prob <- 0.5
    }
    if (is.null(paramList[["catLabel"]])) {
      paramList$catLabel = catLabel
    }
    if (is.null(paramList[["method"]])) {
      paramList$method = method
    }
    if (is.null(paramList[["d"]])) {
      #  q<- min(ceiling(length(y)^0.4),ceiling(paramList$p*2/3))
      #  paramList$d <-min(max(5, ceiling(paramList$p/q)),paramList$p)
      paramList$d <- ceiling(sqrt(paramList[["p"]]))
    }
    
    #if (is.null(paramList[["numProj"]])) {}
    #if (is.null(paramList[["dimProj"]])) {}
    
    
    #################################################
    #pp optimization 
    
    #if (is.null(paramList[["ppMethod"]])) {
    #  paramList$ppMethod="PPR"# "LDA"
    #}
    
    if(1==2){
      if (is.null(paramList[["q"]])) {
        paramList$q=1
      }
      if (is.null(paramList[["weight"]])) {
        paramList$weight = TRUE
      }
      if (is.null(paramList[["r"]])) {
        paramList$r=1
      }
      if (is.null(paramList[["lambda"]])) {
        paramList$lambda=0.1
      }
      if (is.null(paramList[["energy"]])) {
        paramList$energy=0
      }
      if (is.null(paramList[["cooling"]])) {
        paramList$cooling=0.9
      }
      if (is.null(paramList[["TOL"]])) {
        paramList$TOL=0.0001
      }
      if (is.null(paramList[["maxiter"]])) {
        paramList$maxiter=1000
      }
    }
    
    return(paramList)
  }
  
}

if(1==2){
  defaults <- function(ncolX, paramList, cat.map) {
    if (is.null(paramList[["p"]]) || is.na(paramList[["p"]])) {
      paramList$p <- ifelse(is.null(cat.map), ncolX, length(cat.map) +
                              cat.map[[1L]][1L] - 1L)
    }
    
    if (is.null(paramList[["d"]]) || is.na(paramList[["d"]])) {
      paramList$d <- ifelse(is.null(cat.map), ceiling(sqrt(ncolX)),
                            ceiling(sqrt(length(cat.map) + cat.map[[1L]][1L] -
                                           1L)))
    }
    
    if (is.null(paramList[["sparsity"]]) || is.na(paramList[["sparsity"]])) {
      paramList$sparsity <- ifelse(is.null(cat.map), ifelse(ncolX >=
                                                              10, 3/ncolX, 1/ncolX), ifelse(length(cat.map) + cat.map[[1L]][1L] -
                                                                                              1L >= 10, 3/(length(cat.map) + cat.map[[1L]][1L] -
                                                                                                             1L), 1/(length(cat.map) + cat.map[[1L]][1L] - 1L)))
    }
    
    if (is.null(paramList[["prob"]]) || is.na(paramList[["prob"]])) {
      paramList$prob <- 0.5
    }
    
    if (is.null(paramList$catMap)) {
      paramList$catMap = cat.map
    }
    
    return(paramList)
  }
}


#' Create rotation matrix used to determine linear combination of mtry features.
#'
#' This function is the default option to make the projection matrix for
#' unsupervised random forest. The sparseM matrix is the projection
#' matrix.  The creation of this matrix can be changed, but the nrow of
#' sparseM should remain p.  The ncol of the sparseM matrix is currently
#' set to mtry but this can actually be any integer > 1; can even be
#' greater than p.  The matrix returned by this function creates a
#' sparse matrix with multiple features per column.
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param sparsity a real number in \eqn{(0,1)} that specifies the distribution of non-zero elements in the random matrix.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return rotationMatrix the matrix used to determine which mtry features or combination of features will be used to split a node.
#'
#'


makeAB <- function(p, d, sparsity, ...) {
  nnzs <- round(p * d * sparsity)
  sparseM <- matrix(0L, nrow = p, ncol = d)
  featuresToTry <- sample(1:p, d, replace = FALSE)
  # the line below creates linear combinations of features
  # to try
  sparseM[sample(1L:(p * d), nnzs, replace = FALSE)] <- sample(c(1L,
                                                                 -1L), nnzs, replace = TRUE)
  # The below returns a matrix after removing zero columns
  # in sparseM.
  ind <- which(sparseM != 0, arr.ind = TRUE)
  return(cbind(ind, sparseM[ind]))
}



#' Create rotation matrix used to determine mtry features.
#'
#' This function is the default option to make the projection matrix for
#' unsupervised random forest. The sparseM matrix is the projection
#' matrix.  The creation of this matrix can be changed, but the nrow of
#' sparseM should remain p.  The ncol of the sparseM matrix is currently
#' set to mtry but this can actually be any integer > 1; can even be
#' greater than p.  The matrix returned by this function creates a
#' sparse matrix with one feature per column.
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param sparsity a real number in \eqn{(0,1)} that specifies the distribution of non-zero elements in the random matrix.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return rotationMatrix the matrix used to determine which mtry features or combination of features will be used to split a node.
#'
#'


makeA <- function(p, d, sparsity, ...) {
  nnzs <- round(p * d * sparsity)
  sparseM <- matrix(0L, nrow = p, ncol = d)
  featuresToTry <- sample(1:p, d, replace = FALSE)
  # the for loop below creates a standard random forest set
  # of features to try
  for (j in 1:d) {
    sparseM[featuresToTry[j], j] <- 1
  }
  # The below returns a matrix after removing zero columns
  # in sparseM.
  ind <- which(sparseM != 0, arr.ind = TRUE)
  return(cbind(ind, sparseM[ind]))
}


#################################################################################
#' Create a Random Matrix: RandMatPPR
#'
#' @param x an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @param y an n length vector of class labels.  Class labels must be integer or numeric and be within the range 1 to the number of classes.
#' @param d the number of desired columns in the projection matrix.
#' @param catLabel a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom stats rnorm lm ppr cancor glm
#' @importFrom dummies dummy
#' @importFrom dr dr
#'
#' @export
#'
#' @examples
#'
#' x <- matrix(rnorm(200),20,10)
#' y <- (rnorm(20)>0)+0
#' d <- 3
#' sparsity <- 0.25
#' set.seed(220828)
#' (a <- RandMatPPR(x,y,d,sparsity,catMap=NULL))
#@importFrom dcov dcor2d #' @importFrom energy dcor2d

RandMatPPR <- function(x, y, numProj=NULL,catLabel = NULL,method='g-classification',...) {
  p = ncol(x)
  n=length(y)
  #d=min(100, max(5, ceiling(p/q)))
  
  p0=p-(!is.null(catLabel))*(length(unlist(catLabel))-length(catLabel))
  if(is.null(numProj)){
    q =min(ceiling(n^0.4),ceiling(p0*2/3))
    numProj=min(p0,max(5, ceiling(p0/q)))
  }
  d1=numProj
  
  d2=ceiling(sqrt(p0))
  ind <- sample.int(p0, d2, replace = FALSE)
  sparseM = cbind(ind, 1:d2, rep(1, d2))
  
  #ceiling(sqrt(p))#d#round(p/3)#
  #if((n>Q)&(n>10)&(p>1)){
  sparseM1=NULL
  if(n>10){
    #sparseM=NULL
    #for (k in 1:ceiling(2*q*d/p)) {
    ind <- sample.int(p0, replace = FALSE);
    s = c(0,floor(quantile(1:p0, (1:d1)/d1)))
    sparseM1 <-cbind(ind, d2+rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p0))
  }
  
  sparseM=rbind(sparseM,sparseM1)
  if (!is.null(catLabel)) {
    ind=sparseM[,1]; 
    catVar=which(ind<=length(catLabel))
    sparseM[-catVar, 1L]=sparseM[-catVar, 1L]+length(unlist(catLabel))-length(catLabel)
    
    catMap=1
    for (xj in ind[catVar]) {
      isj=which(ind==xj)
      sparseM[isj, 1L]<- sample(catMap:(catMap+length(catLabel[[xj]])-1),length(isj),replace = FALSE)
      catMap=catMap+length(catLabel[[xj]])
    }
  }
  
  ##########################    
  if(n>10){
    Yi=c(y);indC=0L
    if(method!='regression'){
      y=as.factor(y)
      indC=levels(y)
      if(length(indC)>2){
        Yi=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
      }else{
        Yi=as.integer(y)
      }
    }
    
    sparseM1=sparseM[d2+(1:p0),]
    sparseM=sparseM[-(d2+(1:p0)),]
    indTab=table(sparseM1[,2])
    ind=as.numeric(names(indTab)[which(indTab>1)])
    jx=which(sparseM1[,2]%in%ind)
    d11=min(max(c(indTab,d1)),length(unique(sparseM1[jx,1])))
    #d11=min(c(length(unique(sparseM[jx,1])),ifelse(method=='r',Inf,Inf)))#
    #d11=max(max(indTab),min(length(unique(sparseM1[jx,1])),d1)) 
    
    sparseM1=rbind(sparseM1,matrix(d1+d2+1,d11,3))
    for (ni in unique(sparseM1[,2L])) {
      lrows <- which(sparseM1[,2L] == ni)
      j.I <- sparseM1[lrows, 1L]
      if(ni==(d1+d2+1)){
        #jx=which(sparseM1[, 3L]!=1L)
        ix = order(abs(sparseM1[jx,3L]), decreasing = TRUE)
        j.I =unique(sparseM1[jx[ix],1])[1:d11]
        lrows=p0+(1:d11)
        sparseM1[lrows, 1L]=j.I
      }
      Xi <- x[, j.I, drop = FALSE]
      #q <- length(j.I)
      
      #S = eigen(cov(Xi))                # using PCA to handle the colinearity
      #Ii = 1:min(which((cumsum(S$values)+1e-4)/sum(S$values,1e-4) > 0.99))
      #Xi = Xi%*%S$vectors[,Ii, drop = FALSE]
      
      if (length(j.I) > 1L) {
        PPR <- try(ppr(Xi, Yi, nterms = 1, bass=1), silent = TRUE)# sm.method="spline",sm.method="gcvspline"
        
        if(inherits(PPR, "try-error")){
          LM = lm(Yi~.,data=data.frame(Xi))
          if(length(indC)>2){
            
            theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
            for (j in seq(ncol(theta.i))) {
              theta.i[is.na(theta.i[,j]),j] = 0.0
            }
            theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            
          }else{
            theta.i=LM$coefficients[-1]
            theta.i[is.na(theta.i)] = 0.0
          }
          
        }else{
          theta.i=PPR$alpha
        }
        
        #theta.i=S$vectors[,Ii]%*%theta.i
      }else{
        theta.i <- rep(1,length(j.I))
      }
      
      theta.i=theta.i/sqrt(sum(theta.i^2))
      sparseM1[lrows, 3L] <- c(theta.i)
    }
    
    #sparseM1[,2]=sparseM1[,2]+d2
    #sparseM[d2+(1:p0),]=sparseM1
    sparseM=rbind(sparseM,sparseM1)
  }
  
  return(sparseM)
}


#' Create a Random Matrix: RandMatPPR
#'
#' @param x an n by d numeric matrix (preferable) or data frame. The rows correspond to observations and columns correspond to features.
#' @param y an n length vector of class labels.  Class labels must be integer or numeric and be within the range 1 to the number of classes.
#' @param d the number of desired columns in the projection matrix.
#' @param catLabel a list specifying specifies which one-of-K encoded columns in X correspond to the same categorical feature.
#' @param ... used to handle superfluous arguments passed in using paramList.
#'
#' @return A random matrix to use in running \code{\link{RerF}}.
#'
#' @importFrom stats rnorm lm ppr cancor glm
#' @importFrom dummies dummy
#' @importFrom dr dr
#' @importFrom nnet nnet
#' @importFrom neuralnet neuralnet
#' @export
#'
#' @examples
#'
#' x <- matrix(rnorm(200),20,10)
#' y <- (rnorm(20)>0)+0
#' d <- 3
#' sparsity <- 0.25
#' set.seed(220828)
#' (a <- RandMatPPR(x,y,d,sparsity,catMap=NULL))
#@importFrom dcov dcor2d q=NUll numProj="Rand" #' @importFrom energy dcor2d
RandMatPPO <- function(x, y, numProj=NULL,dimProj=NULL,catLabel = NULL,ppMethod="NNet",method='i-classification',...) {
  p = ncol(x)
  n=length(y)
  #d=min(100, max(5, ceiling(p/q))) d q
  p0=p-(!is.null(catLabel))*(length(unlist(catLabel))-length(catLabel))
  
  d2=ceiling(sqrt(p0))
  ind <- sample(p0, d2, replace = FALSE)
  sparseM = cbind(ind, 1:d2, rep(1, d2))
  
  #ceiling(sqrt(p))#d#round(p/3)#
  #if((n>Q)&(n>10)&(p>1)){
  sparseM1=NULL
  if(n>10){
    if(identical(dimProj,"Rand")){
      if(is.null(numProj)){max(5,sample(floor(p0/3),1))}
      numProj=min(p0,numProj)
      d1=numProj
      spd=sample(p0,d1)
      indp=sum(spd)
      ind=unlist(sapply(spd,function(pd)sample(p0,pd)))
      sparseM1 <-cbind(ind, d2+rep(1:d1,spd) ,rep(1, indp))
    }else{
      if(is.null(dimProj)){dimProj =min(ceiling(n^0.4),ceiling(p0*2/3))}
      if(is.null(numProj)){numProj=max(5, ceiling(p0/dimProj))} 
      numProj=min(p0,numProj)
      d1=numProj
      
      indp=p0
      #sparseM=NULL
      #for (k in 1:ceiling(2*q*d/p)) {
      ind <- sample(1:indp, replace = FALSE);
      s = c(0,floor(quantile(1:indp, (1:d1)/d1)))
      sparseM1 <-cbind(ind, d2+rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, indp))
    }
  }
  
  sparseM=rbind(sparseM,sparseM1)
  if (!is.null(catLabel)) {
    ind=sparseM[,1]; 
    catVar=which(ind<=length(catLabel))
    sparseM[-catVar, 1L]=sparseM[-catVar, 1L]+length(unlist(catLabel))-length(catLabel)
    
    catMap=1
    for (xj in ind[catVar]) {
      isj=which(ind==xj)
      sparseM[isj, 1L]<- sample(catMap:(catMap+length(catLabel[[xj]])-1),length(isj),replace = FALSE)
      catMap=catMap+length(catLabel[[xj]])
    }
  }
  
  ##########################    
  if(n>10){
    Yi=c(y);indC=0L
    if(method!='regression'){
      y=as.factor(y)
      indC=levels(y)
      if(length(indC)>2){
        Yi=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
      }else{
        Yi=as.integer(y)
      }
    }
    
    sparseM1=sparseM[d2+(1:indp),]
    sparseM=sparseM[-(d2+(1:indp)),]
    indTab=table(sparseM1[,2])
    ind=as.numeric(names(indTab)[which(indTab>1)])
    jx=which(sparseM1[,2]%in%ind)
    d11=min(max(c(indTab,d1)),length(unique(sparseM1[jx,1])))
    #d11=min(c(length(unique(sparseM[jx,1])),ifelse(method=='r',Inf,Inf)))#
    #d11=max(max(indTab),min(length(unique(sparseM1[jx,1])),d1)) 
    
    sparseM1=rbind(sparseM1,matrix(d1+d2+1,d11,3))
    for (ni in unique(sparseM1[,2L])) {
      lrows <- which(sparseM1[,2L] == ni)
      j.I <- sparseM1[lrows, 1L]
      if(ni==(d1+d2+1)){
        #jx=which(sparseM1[, 3L]!=1L)
        ix = order(abs(sparseM1[jx,3L]), decreasing = TRUE)
        j.I =unique(sparseM1[jx[ix],1])[1:d11]
        lrows=indp+(1:d11)
        sparseM1[lrows, 1L]=j.I
      }
      Xi <- x[, j.I, drop = FALSE]
      pi <- length(j.I)
      
      #S = eigen(cov(Xi))                # using PCA to handle the colinearity
      #Ii = 1:min(which((cumsum(S$values)+1e-4)/sum(S$values,1e-4) > 0.99))
      #Xi = Xi%*%S$vectors[,Ii, drop = FALSE]
      
      if (pi > 1L) {
        if(ppMethod=="PPR"){
          PP <- try(ppr(Xi, Yi, nterms = 1, bass=1)$alpha, silent = TRUE)# sm.method="spline",sm.method="gcvspline"
        }
        if(ppMethod=="Rand"){
          PP <- sample(c(1L, -1L),ncol(Xi), replace = TRUE,prob = c(0.5, 0.5))
        }
        
        if(ppMethod=="NNet"){
          #PP = nnet(Xi, Yi, size=1, linout=TRUE, trace=FALSE)$wts[2:(1+ncol(Xi))]
          #PP =neuralnet(Yi~.,data=data.frame(Xi))$weights[[1]][[1]][-1]
          #dat=data.frame(Xi,Yi);
          #colnames(dat)=c(paste0('x',seq(length(j.I))),paste0('y',seq(length(indC)-1)))
          #fm=paste0(paste(colnames(dat)[-seq(length(j.I))],collapse ="+"),"~",paste(colnames(dat)[seq(length(j.I))],collapse ="+"))
          #PP =neuralnet(fm,dat)$weights[[1]][[1]][-1]
          #nn <- neuralnet((Species == "setosa") + (Species == "versicolor") + (Species == "virginica")
          #                ~ Petal.Length + Petal.Width, iris_train, linear.output = FALSE)
          #library(NeuralNetTools)
          #data(neuraldat)
          #library(RSNNS)
          #x <- neuraldat[, c('X1', 'X2', 'X3')]
          #y <- neuraldat[, 'Y1']
          #mod <- mlp(x, y, size = 5, linOut = TRUE)
          #neuralweights(mod)
          
          #if((n>5*pi)&(pi<10)){
          #  PP <- try(ppr(Xi, Yi, nterms = 1, bass=1)$alpha, silent = TRUE)
          #}else{
            #PP = try(nnet(Xi, y, size=1,trace=FALSE)$wts[2:(1+pi)], silent = TRUE)
            PP = nnet(Xi, Yi, size=1,linout=TRUE,trace=FALSE)$wts[2:(1+pi)]#  
          #}
        }
        
        if(ppMethod=="NeuralNet"){
          #if((n>5*pi)&(pi<10)){
          #  PP <- try(ppr(Xi, Yi, nterms = 1, bass=1)$alpha, silent = TRUE)
          #}else{
            #PP =neuralnet(Yi~.,data=data.frame(Xi))$weights[[1]][[1]][-1]
            dat=data.frame(Xi,Yi);
            if(length(indC)>2){
              colnames(dat)=c(paste0('x',seq(pi)),paste0('y',seq(length(indC))))
              fm=paste0(paste(colnames(dat)[-seq(pi)],collapse ="+"),"~",paste(colnames(dat)[seq(pi)],collapse ="+"))
            }else{
              colnames(dat)=c(paste0('x',seq(pi)),'y')
              fm=paste0('y',"~",paste(colnames(dat)[seq(pi)],collapse ="+"))
            }
            PP =neuralnet(formula(fm),dat)$weights[[1]][[1]][-1]
          #}
        }
        
        if(!ppMethod%in%c("PPR","Rand","NNet","NeuralNet")){
          PP <- ppRF:::ppOptCpp(y,Xi,q = 1L, PPmethod = ppMethod, weight = TRUE, r = 1L, lambda = 0.1, 
                                    energy = 0, cooling = 0.9, TOL = 0.001, maxiter = 1000L)$projbest
        }

        if(inherits(PP, "try-error")){
          LM = lm(Yi~.,data=data.frame(Xi))
          if(length(indC)>2){
            
            theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
            for (j in seq(ncol(theta.i))) {
              theta.i[is.na(theta.i[,j]),j] = 0.0
            }
            theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            
          }else{
            theta.i=LM$coefficients[-1]
            theta.i[is.na(theta.i)] = 0.0
          }
          
        }else{
          theta.i=PP
        }
        
        #theta.i=S$vectors[,Ii]%*%theta.i
      }else{
        theta.i <- rep(1,pi)
      }
      
      theta.i=theta.i/sqrt(sum(theta.i^2))
      sparseM1[lrows, 3L] <- c(theta.i)
    }
    
    #sparseM1[,2]=sparseM1[,2]+d2
    #sparseM[d2+(1:p0),]=sparseM1
    sparseM=rbind(sparseM,sparseM1)
  }
  
  return(sparseM)
}
###############################################################################

if(1==2){
  RandMatPPR <- function(x, y, d=NULL,catLabel = NULL,method='g-classification',
                         lambda=c(5,log(length(y)))[2],...) {
    p = ncol(x)
    n=length(y)
    #d=min(100, max(5, ceiling(p/q)))
    
    p0=p-(!is.null(catLabel))*(length(unlist(catLabel))-length(catLabel))
    if(is.null(d)){
      q =min(ceiling(n^0.4),ceiling(p0*2/3))
      d=min(p0,max(5, ceiling(p0/q)))
    }
    d1=d
    
    #d2=ceiling(sqrt(p0))
    d2=d
    ind <- sample.int(p0, d2, replace = FALSE)
    sparseM = cbind(ind, 1:d2, rep(1, d2))
    
    #ceiling(sqrt(p))#d#round(p/3)#
    #if((n>Q)&(n>10)&(p>1)){
    sparseM1=NULL
    if(n>10){
      #sparseM=NULL
      #for (k in 1:ceiling(2*q*d/p)) {
      ind <- sample.int(p0, replace = FALSE);
      s = c(0,floor(quantile(1:p0, (1:d1)/d1)))
      sparseM1 <-cbind(ind, d2+rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p0))
    }
    
    sparseM=rbind(sparseM,sparseM1)
    if (!is.null(catLabel)) {
      ind=sparseM[,1]; 
      catVar=which(ind<=length(catLabel))
      sparseM[-catVar, 1L]=sparseM[-catVar, 1L]+length(unlist(catLabel))-length(catLabel)
      
      catMap=1
      for (xj in ind[catVar]) {
        isj=which(ind==xj)
        sparseM[isj, 1L]<- sample(catMap:(catMap+length(catLabel[[xj]])-1),length(isj),replace = FALSE)
        catMap=catMap+length(catLabel[[xj]])
      }
    }
    
    ##########################    
    if(n>10){
      Yi=c(y);indC=0L
      if(method!='regression'){
        y=as.factor(y)
        indC=levels(y)
        if(length(indC)>2){
          Yi=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
        }else{
          Yi=as.integer(y)
        }
      }
      
      sparseM1=sparseM[d2+(1:p0),]
      sparseM=sparseM[-(d2+(1:p0)),]
      indTab=table(sparseM1[,2])
      ind=as.numeric(names(indTab)[which(indTab>1)])
      jx=which(sparseM1[,2]%in%ind)
      d11=min(max(c(indTab,d1)),length(unique(sparseM1[jx,1])))
      #d11=min(c(length(unique(sparseM[jx,1])),ifelse(method=='r',Inf,Inf)))#
      #d11=max(max(indTab),min(length(unique(sparseM1[jx,1])),d1)) 
      
      sparseM1=rbind(sparseM1,matrix(d1+d2+1,d11,3))
      #sparseM1=rbind(sparseM1,matrix(d1+d2+2,p,3))
      for (ni in unique(sparseM1[,2L])) {
        if(ni==(d1+d2+2)){
          #nnzs <- round(p0 * d * sparsity)
          #m = max(20, floor(p^0.5))
          prob=0.5
          B=matrix(sample(c(1L, -1L), p*p, replace = TRUE,prob = c(prob, 1 - prob)),p,p) 
          I=apply(x%*%B,2,dcov::dcor2d,as.integer(y))
          I = order(I, decreasing = TRUE)[1:d11]
          B = B[,I]
          Xi = x%*%B
          j.I <- 1:p
          lrows=p0+d11+j.I
          sparseM1[lrows, 1L]=j.I
          
          S = eigen(cov(Xi))  # using PCA to handle the colinearity
          Ii = 1:min(which((cumsum(S$values)+1e-4)/sum(S$values,1e-4) > 0.99))
          if(length(Ii)>1){S=S$vectors[,Ii];Xi = Xi%*%S}else{S=diag(ncol(Xi))}
        }else{
          lrows <- which(sparseM1[,2L] == ni)
          j.I <- sparseM1[lrows, 1L]
          if(ni==(d1+d2+1)){
            #jx=which(sparseM1[, 3L]!=1L)
            ix = order(abs(sparseM1[jx,3L]), decreasing = TRUE)
            j.I =unique(sparseM1[jx[ix],1])[1:d11]
            lrows=p0+(1:d11)
            sparseM1[lrows, 1L]=j.I
          }
          Xi <- x[, j.I, drop = FALSE]
        }
        
        if (ncol(Xi) > 1L) {
          PPR <- try(ppr(Xi, Yi, nterms = 1, bass=1), silent = TRUE)# sm.method="spline",sm.method="gcvspline"
          
          if(inherits(PPR, "try-error")){
            LM = lm(Yi~.,data=data.frame(Xi))
            if(length(indC)>2){
              
              theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i))) {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
              theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
              
            }else{
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
            
          }else{
            theta.i=PPR$alpha
          }
          
          if(ni==(d1+d2+2)){theta.i=B%*%S%*%theta.i}#theta.i=S%*%theta.i;
          
        }else{
          theta.i <- rep(1,length(j.I))
        }
        
        theta.i=theta.i/sqrt(sum(theta.i^2))
        sparseM1[lrows, 3L] <- c(theta.i)
      }
      
      #sparseM1[,2]=sparseM1[,2]+d2
      #sparseM[d2+(1:p0),]=sparseM1
      sparseM=rbind(sparseM,sparseM1)
    }
    
    
    ######################################
    if(1==1){
      if(nrow(sparseM)>d2){
        #library(ppRF)
        
        numDr=unique(sparseM[,2]);
        rotaX=matrix(0,length(numDr),p);
        for(i in 1:length(numDr)){
          lrows = which(sparseM[,2] == numDr[i]);
          rotaX[i,sparseM[lrows, 1]] = sparseM[lrows, 3];
        }
        rotaX=x %*% t(rotaX);
        
        BestVal=BestCutNode(rotaX,y,method)$BestVal
        BestVal=log(BestVal+1/n^2)
        #MinVal[-seq(d2)]=log(MinVal[-seq(d2)]+1e-4) + 2*table(sparseM[-seq(d2),2])/n
        if(method=="i-classification"){
          BestVal[-seq(d2)]=BestVal[-seq(d2)]-lambda*table(sparseM[-seq(d2),2])/n
          bestVar=which.max(BestVal)
        }else{
          BestVal[-seq(d2)]=BestVal[-seq(d2)]+lambda*table(sparseM[-seq(d2),2])/n#ifelse(n>20,n^2,sqrt(n))  
          #BestVal[-seq(d2)]=BestVal[-seq(d2)]+lambda*(which.min(BestVal)>d2)*table(sparseM[-seq(d2),2])/n*
          #  (min(BestVal[seq(d2)])-min(BestVal[-seq(d2)]))
          bestVar=which.min(BestVal)
        }
        
        
        sparseM=sparseM[sparseM[,2]==bestVar,,drop = FALSE]
        sparseM[,2]=1
      }
    }
    return(sparseM)
  }
}

if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(sqrt(length(y)),ncol(x)/3)),
                         sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)),
                         catMap = NULL,method='c',...) {
    #options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    #d1 = round(d/3)
    #d2 = round(d/3)
    #d3 = d - d1-d2
    
    #d1=round(d/2)
    #d3 = d - d1
    d1=d
    d3=d
    #ceiling(sqrt(p))#d#round(p/3)#
    
    p = ncol(x)
    n=length(y)
    #if(1==2){
    #if(d>2){
    #if(n>(p/d)^2){
    #if((n>p^2)&(p>1))
    #if((n>p)&(n>10)&(p>1)){
    if((n>ceiling(p/d))&(n>10)&(p>1)){
      ind <- sample.int(p, replace = FALSE);
      s = c(0,round(quantile(1:p, (1:d1)/d1)))
      sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p))
      nnzs=p
      
      #prob=0.5;d2=d1
      #nnzs <- round(p * d1 * sparsity)
      #ind <- sort(sample.int((p * d1), nnzs, replace = FALSE))
      #sparseM1 <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, rep(1,nnzs))
      
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM3 = cbind(ind, (1:d3) + d1, rep(1, d3))
      
      sparseM=rbind(sparseM1,sparseM3)#sparseM1#
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      Yi=c(y);indC=0L
      if(method!='r'){
        y=as.factor(y)
        indC=levels(y)
        if(length(indC)>2){
          Yi=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
        }else{
          Yi=as.integer(y)
        }
      }
      
      gcv=Inf; 
      ind=table(sparseM[1:nnzs,2])
      d11=max(ind)
      ind=as.numeric(names(ind)[which(ind>1)])
      jx=which(sparseM[1:nnzs,2]%in%ind)
      #d11=min(c(length(unique(sparseM[jx,1])),ifelse(method=='r',Inf,Inf)))#
      d11=min(c(length(unique(sparseM[jx,1])),d1))
      
      sparseM1=rbind(sparseM[1:nnzs,],matrix(d1+1,d11,3))
      for (ni in unique(sparseM1[,2L])) {
        lrows <- which(sparseM1[,2L] == ni)
        j.I <- sparseM1[lrows, 1L]
        if(ni==(d1+1)){
          #jx=which(sparseM1[1:nnzs, 3L]!=1L)
          ix = order(abs(sparseM1[jx,3L]), decreasing = TRUE)
          j.I =unique(sparseM1[jx[ix],1])[1:d11]
          lrows=nnzs+(1:d11)
          sparseM1[lrows, 1L]=j.I
        }
        Xi <- x[, j.I, drop = FALSE]
        q <- length(j.I)
        
        if (q > 1L) {
          PPR <- try(ppr(Xi, Yi, nterms = 1, bass=1), silent = TRUE)# sm.method="spline",sm.method="gcvspline"
          #tryCatch({
          #  PPR <- ppr(Xi, Yi, nterms = 1)
          #}, error = function(e) {
          #  cat("ERROR:", e$message, "\n")
          #})
          
          if(inherits(PPR, "try-error")){
            LM = lm(Yi~.,data=data.frame(Xi))
            if(length(indC)>2){
              theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i))) {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
              theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            }else{
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
            gcv.i=sum(LM$residuals^2)
          }else{
            theta.i=PPR$alpha
            gcv.i=sum(PPR$residuals^2)
            #fitted=indC[apply(PPR$fitted.values,1,which.max)]
            #gcv.i = mean(fitted!=y)
          }
        }else{
          theta.i <- 1L
          gcv.i=Inf
        }
        #theta.i=theta.i/sqrt(sum(theta.i^2))
        sparseM1[lrows, 3L] <- c(theta.i)
        
        if(gcv.i<gcv){
          #theta=theta.i
          sparseM0=sparseM1[lrows,,drop = FALSE]
          #sparseM0[lrows, 3L] <- c(theta.i)
          gcv=gcv.i
        }
      }
      #sparseM1=sparseM1[-(1:nnzs),];sparseM1[,2]=1
      #sparseM[-(1:nnzs),2]=sparseM[-(1:nnzs),2]+1
      #sparseM=rbind(sparseM1[-(1:nnzs),],sparseM[-(1:nnzs),])
      #sparseM[,2]=c(rep(1,d11),2:(d3+1))
      
      sparseM[-(1:nnzs),2]=sparseM[-(1:nnzs),2]+1
      sparseM=rbind(sparseM1,sparseM[-(1:nnzs),])
      
      #sparseM=rbind(sparseM0,sparseM[-(1:nnzs),])
      #sparseM[,2]=c(rep(1,nrow(sparseM0)),2:(d3+1))
      
      #sparseM=sparseM0
      #sparseM[,2]=1
      #sparseM=sparseM[-(1:p),]
      #sparseM[,2]=1:d3
    }else{
      #prob=0.5
      #nnzs <- round(p * d * sparsity)
      #ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
      #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      d3=ceiling(sqrt(p))
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM = cbind(ind, 1:d3, rep(1, d3))
    }
    
    #options (warn = 1)
    return(sparseM)
  }
  
}

if(1==2){
  RandMatPPR01 <- function(x, y, q = min(ceiling(length(y)^0.4),ceiling(ncol(x)*2/3)),
                           sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)),
                           catMap = NULL,method='c',...) {
    #options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    #d1 = round(d/3)
    #d2 = round(d/3)
    #d3 = d - d1-d2
    
    p = ncol(x)
    n=length(y)
    #d=min(100, max(5, ceiling(p/q)))
    d=max(5, ceiling(p/q))
    #d1=round(d/2)
    #d3 = d - d1
    d1=d
    d2=d
    #ceiling(sqrt(p))#d#round(p/3)#
    
    
    #if(1==2){
    #if(d>2){
    #if(n>(p/d)^2){
    #if((n>p^2)&(p>1))
    #if((n>p)&(n>10)&(p>1)){
    if((n>q)&(n>10)&(p>1)){
      
      Yi=c(y);indC=0L
      if(method!='r'){
        y=as.factor(y)
        indC=levels(y)
        if(length(indC)>2){
          Yi=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
        }else{
          Yi=as.integer(y)
        }
      }
      
      sparseM=NULL
      for (k in 1:ceiling(2*q*d/p)) {
        ind <- sample.int(p, replace = FALSE);
        s = c(0,round(quantile(1:p, (1:d1)/d1)))
        sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p))
        nnzs=p
        
        #prob=0.5;d2=d1
        #nnzs <- round(p * d1 * sparsity)
        #ind <- sort(sample.int((p * d1), nnzs, replace = FALSE))
        #sparseM1 <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, rep(1,nnzs))
        
        #sparseM=sparseM1#rbind(sparseM1,sparseM3)#
        if (!is.null(catMap)) {
          rw <- sparseM1[, 1]
          for (j in 1:length(catMap)) {
            isj <- rw %in% catMap[[j]]
            rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
          }
          sparseM1[, 1L] <- rw
        }
        
        gcv=Inf; 
        ind=table(sparseM1[,2])
        #d11=max(ind)
        ind=as.numeric(names(ind)[which(ind>1)])
        jx=which(sparseM1[,2]%in%ind)
        #d11=min(c(length(unique(sparseM[jx,1])),ifelse(method=='r',Inf,Inf)))#
        d11=min(c(length(unique(sparseM1[jx,1])),d1))
        
        sparseM1=rbind(sparseM1,matrix(d1+1,d11,3))
        for (ni in unique(sparseM1[,2L])) {
          lrows <- which(sparseM1[,2L] == ni)
          j.I <- sparseM1[lrows, 1L]
          if(ni==(d1+1)){
            #jx=which(sparseM1[, 3L]!=1L)
            ix = order(abs(sparseM1[jx,3L]), decreasing = TRUE)
            j.I =unique(sparseM1[jx[ix],1])[1:d11]
            lrows=nnzs+(1:d11)
            sparseM1[lrows, 1L]=j.I
          }
          Xi <- x[, j.I, drop = FALSE]
          q <- length(j.I)
          
          if (q > 1L) {
            PPR <- try(ppr(Xi, Yi, nterms = 1, bass=1), silent = TRUE)# sm.method="spline",sm.method="gcvspline"
            #tryCatch({
            #  PPR <- ppr(Xi, Yi, nterms = 1)
            #}, error = function(e) {
            #  cat("ERROR:", e$message, "\n")
            #})
            
            if(inherits(PPR, "try-error")){
              LM = lm(Yi~.,data=data.frame(Xi))
              if(length(indC)>2){
                theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
                for (j in seq(ncol(theta.i))) {
                  theta.i[is.na(theta.i[,j]),j] = 0.0
                }
                theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
              }else{
                theta.i=LM$coefficients[-1]
                theta.i[is.na(theta.i)] = 0.0
              }
              gcv.i=sum(LM$residuals^2)
            }else{
              theta.i=PPR$alpha
              gcv.i=sum(PPR$residuals^2)
              #fitted=indC[apply(PPR$fitted.values,1,which.max)]
              #gcv.i = mean(fitted!=y)
            }
          }else{
            theta.i <- 1L
            gcv.i=Inf
          }
          #theta.i=theta.i/sqrt(sum(theta.i^2))
          sparseM1[lrows, 3L] <- c(theta.i)
          
          if(gcv.i<gcv){
            #theta=theta.i
            sparseM0=sparseM1[lrows,,drop = FALSE]
            #sparseM0[lrows, 3L] <- c(theta.i)
            gcv=gcv.i
          }
        }
        
        sparseM1[,2]=sparseM1[,2]+(k-1)*(d+1)
        sparseM=rbind(sparseM,sparseM1)
      }
      
      ind <- sample.int(p, d2, replace = FALSE)
      sparseM1 = cbind(ind, (1:d2) + k*(d+1), rep(1, d2))
      sparseM=rbind(sparseM,sparseM1)
      
      #sparseM[-(1:nnzs),2]=sparseM[-(1:nnzs),2]+1
      #sparseM=rbind(sparseM1,sparseM[-(1:nnzs),])
      
      #sparseM1=sparseM1[-(1:nnzs),];sparseM1[,2]=1
      #sparseM[-(1:nnzs),2]=sparseM[-(1:nnzs),2]+1
      #sparseM=rbind(sparseM1[-(1:nnzs),],sparseM[-(1:nnzs),])
      #sparseM[,2]=c(rep(1,d11),2:(d3+1))
      
      #sparseM=rbind(sparseM0,sparseM[-(1:nnzs),])
      #sparseM[,2]=c(rep(1,nrow(sparseM0)),2:(d3+1))
      
      #sparseM=sparseM0
      #sparseM[,2]=1
      #sparseM=sparseM[-(1:p),]
      #sparseM[,2]=1:d3
    }else{
      #prob=0.5
      #nnzs <- round(p * d * sparsity)
      #ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
      #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      d2=ceiling(sqrt(p))
      ind <- sample.int(p, d2, replace = FALSE)
      sparseM = cbind(ind, 1:d2, rep(1, d2))
    }
    
    #options (warn = 1)
    return(sparseM)
  }
}

if(1==2){
  RandMatPPR1 <- function(x, y, method='c',...) {
    p = ncol(x)
    n=length(y)
    
    sparseM = rep(1,p)
    #if((n>p^2)&(p>1))
    if((n>p)&(n>10)&(p>1))
    {
      Y=c(y);indC=0L
      if(method!='r'){
        y=as.factor(y)
        indC=levels(y)
        if(length(indC)>2){
          Y=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+1
          +0.0001*matrix(runif(n*length(indC),-1,1),n,length(indC)) 
        }else{
          Y=as.integer(y)+0.0001*runif(n,-1,1)
        }
      }
      
      sparseM <- try(ppr(x, Y, bass=1,nterms = 1)$alpha, silent = TRUE)# sm.method="spline",sm.method="gcvspline"
      if(inherits(sparseM, "try-error")){
        
        LM = lm(Y~.,data=data.frame(x))
        if(length(indC)>2){
          
          sparseM=as.matrix(LM$coefficients)[-1,,drop = FALSE]
          for (j in seq(ncol(sparseM))) {
            sparseM[is.na(sparseM[,j]),j] = 0.0
          }
          
          sparseM = eigen(sparseM%*%t(sparseM))$vectors[,1]
        }else{
          sparseM=LM$coefficients[-1]
          sparseM[is.na(sparseM)] = 0.0
        }
      }
    }
    
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(sqrt(length(y)),ncol(x)/3)),
                         sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)),
                         catMap = NULL,method='c',...) {
    #options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    #d1 = round(d/3)
    #d2 = round(d/3)
    #d3 = d - d1-d2
    
    #d1=round(d/2)
    #d3 = d - d1
    d1=d
    d3=d#ceiling(sqrt(p))#d#round(p/3)#
    
    p = ncol(x)
    n=length(y)
    #if(1==2){
    if(n>(p/d)^2){
      ind <- sample.int(p, replace = FALSE);
      s = c(0,round(quantile(1:p, (1:d1)/d1)))
      sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p))
      
      #prob=0.5
      #nnzs <- round(p * d2 * sparsity)
      #ind <- sort(sample.int((p * d2), nnzs, replace = FALSE))
      #sparseM2 <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L+d1, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM3 = cbind(ind, (1:d3) + d1, rep(1, d3))
      
      sparseM=rbind(sparseM1,sparseM3)#sparseM1#
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      Yi=c(y);indC=0L
      if(method!='r'){
        y=as.factor(y)
        indC=levels(y)
        if(length(indC)>2){
          Yi=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
        }else{
          Yi=as.integer(y)
        }
      }
      
      gcv=Inf
      sparseM1=rbind(sparseM[1:p,],matrix(d1+1,d1,3))
      sparseM[-(1:p),2]=sparseM[-(1:p),2]+1
      for (ni in 1:(d1+1)) {
        lrows <- which(sparseM1[,2L] == ni)
        j.I <- sparseM1[lrows, 1L]
        if(ni==(d1+1)){
          ix = order(abs(sparseM1[1:p, 3L]), decreasing = TRUE)
          j.I =sparseM1[ix,1][1:d1]
          lrows=p+(1:d1)
          sparseM1[lrows, 1L]=j.I
        }
        Xi <- x[, j.I, drop = FALSE]
        q <- length(j.I)
        
        if (q > 1L) {
          PPR <- try(ppr(Xi, Yi, nterms = 1, bass=1), silent = TRUE)# sm.method="spline",sm.method="gcvspline"
          #tryCatch({
          #  PPR <- ppr(Xi, Yi, nterms = 1)
          #}, error = function(e) {
          #  cat("ERROR:", e$message, "\n")
          #})
          
          if(inherits(PPR, "try-error")){
            LM = lm(Yi~.,data=data.frame(Xi))
            if(length(indC)>2){
              theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i))) {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
              theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            }else{
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
            gcv.i=sum(LM$residuals^2)
          }else{
            theta.i=PPR$alpha
            gcv.i=sum(PPR$residuals^2)
            #fitted=indC[apply(PPR$fitted.values,1,which.max)]
            #gcv.i = mean(fitted!=y)
          }
        }else{
          theta.i <- 1L
          gcv.i=Inf
        }
        #theta.i=theta.i/sqrt(sum(theta.i^2))
        sparseM1[lrows, 3L] <- c(theta.i)
        
        if(gcv.i<gcv){
          #theta=theta.i
          sparseM0=sparseM1[lrows,]
          #sparseM0[lrows, 3L] <- c(theta.i)
          gcv=gcv.i
        }
      }
      sparseM=rbind(sparseM1,sparseM[-(1:p),])
      #sparseM=rbind(sparseM0,sparseM[-(1:p),])
      #sparseM[,2]=c(rep(1,length(sparseM0[,2])),2:(d3+1))
      
      #sparseM=sparseM0
      #sparseM[,2]=1
      #sparseM=sparseM[-(1:p),]
      #sparseM[,2]=1:d3
    }else{
      #prob=0.5
      #nnzs <- round(p * d * sparsity)
      #ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
      #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      d3=ceiling(sqrt(p))
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM = cbind(ind, 1:d3, rep(1, d3))
    }
    
    #options (warn = 1)
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(sqrt(length(y)),ncol(x)/3)),
                         sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)),
                         catMap = NULL,method='c',...) {
    #options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    #d1 = round(d/3)
    #d2 = round(d/3)
    #d3 = d - d1-d2
    
    #d1=round(d/2)
    #d3 = d - d1
    d1=d
    d3=d#ceiling(sqrt(p))#d#round(p/3)#
    
    p = ncol(x)
    n=length(y)
    #if(1==2){
    if(n>(p/d)^2){
      ind <- sample.int(p, replace = FALSE);
      s = c(0,round(quantile(1:p, (1:d1)/d1)))
      sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p))
      
      #prob=0.5
      #nnzs <- round(p * d2 * sparsity)
      #ind <- sort(sample.int((p * d2), nnzs, replace = FALSE))
      #sparseM2 <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L+d1, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM3 = cbind(ind, (1:d3) + d1, rep(1, d3))
      
      sparseM=rbind(sparseM1,sparseM3)#sparseM1#
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      Y=c(y);indC=0L
      if(method!='r'){
        y=as.factor(y)
        indC=levels(y)
        if(length(indC)>2){
          Y=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
        }else{
          Y=as.integer(y)
        }
      }
      Y=as.matrix(Y)
      
      
      sparseMi=NULL
      for (i in seq(ncol(Y))) {
        Yi=Y[,i]
        
        sparseM1=rbind(sparseM[1:p,],matrix(d1+1,d1,3))
        gcv=Inf
        for (ni in 1:(d1+1)) {
          lrows <- which(sparseM1[,2L] == ni)
          j.I <- sparseM1[lrows, 1L]
          if(ni==(d1+1)){
            ix = order(abs(sparseM1[1:p, 3L]), decreasing = TRUE)
            j.I =sparseM1[ix,1][1:d1]
            lrows=p+(1:d1)
            sparseM1[lrows, 1L]=j.I
          }
          Xi <- x[, j.I, drop = FALSE]
          q <- length(j.I)
          
          if (q > 1L) {
            PPR <- try(ppr(Xi, Yi, nterms = 1), silent = TRUE)# sm.method="spline"
            #tryCatch({
            #  PPR <- ppr(Xi, Yi, nterms = 1)
            #}, error = function(e) {
            #  cat("ERROR:", e$message, "\n")
            #})
            
            if(inherits(PPR, "try-error")){
              LM = lm(Yi~.,data=data.frame(Xi))
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
              gcv.i=sum(LM$residuals^2)
            }else{
              theta.i=PPR$alpha
              gcv.i=sum(PPR$residuals^2)
              #fitted=indC[apply(PPR$fitted.values,1,which.max)]
              #gcv.i = mean(fitted!=y)
            }
          }else{
            theta.i <- 1L
            gcv.i=Inf
          }
          #theta.i=theta.i/sqrt(sum(theta.i^2))
          sparseM1[lrows, 3L] <- c(theta.i)
          
          if(gcv.i<gcv){
            #theta=theta.i
            sparseM0=sparseM1[lrows,]
            #sparseM0[lrows, 3L] <- c(theta.i)
            gcv=gcv.i
          }
        }
        #sparseM=rbind(sparseM1,sparseM[-(1:p),])
        #sparseM=rbind(sparseM0,sparseM[-(1:p),])
        #sparseM[,2]=c(rep(1,length(sparseM0[,2])),2:(d3+1))
        
        sparseM1[,2]=sparseM1[,2]+d1*(i-1)
        sparseMi=rbind(sparseMi,sparseM1)
      }
      sparseM[-(1:p),2]=max(sparseMi[,2])+c(1:d3)
      sparseM=rbind(sparseMi,sparseM[-(1:p),])
      #sparseM=sparseM0
      #sparseM[,2]=1
      #sparseM=sparseM[-(1:p),]
      #sparseM[,2]=1:d3
    }else{
      #prob=0.5
      #nnzs <- round(p * d * sparsity)
      #ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
      #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      d3=d#ceiling(sqrt(p))
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM = cbind(ind, 1:d3, rep(1, d3))
    }
    
    #options (warn = 1)
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(sqrt(length(y)),
                                               ncol(x)/3)), sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),
                                                                              1/ncol(x)), catMap = NULL, ...) {
    #options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    
    p = ncol(x)
    n=length(y)
    #if(1==2){
    if(d>=3){
      #d1 = round(d/3)
      #d2 = round(d/3)
      #d3 = d - d1-d2
      
      #d1=round(d/2)
      #d3 = d - d1
      d1=d
      d3=ceiling(sqrt(p))#round(p/3)
      
      ind <- sample.int(p, replace = FALSE);
      s = c(0,round(quantile(1:p, (1:d1)/d1)))
      sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p))
      
      #prob=0.5
      #nnzs <- round(p * d2 * sparsity)
      #ind <- sort(sample.int((p * d2), nnzs, replace = FALSE))
      #sparseM2 <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L+d1, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM3 = cbind(ind, (1:d3) + d1, rep(1, d3))
      
      sparseM=rbind(sparseM1,sparseM3)#sparseM1#
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      indC=levels(as.factor(y))
      if(length(indC)>2){
        Yi=(matrix(y,n,length(indC))==matrix(as.numeric(as.character(indC)),n,length(indC),byrow = TRUE))+0
      }else{
        Yi=as.numeric(as.character(y))
      }
      
      gcv=Inf
      sparseM1=rbind(sparseM[1:p,],matrix(d1+1,d1,3))
      sparseM[-(1:p),2]=sparseM[-(1:p),2]+1
      for (ni in 1:(d1+1)) {
        lrows <- which(sparseM1[,2L] == ni)
        j.I <- sparseM1[lrows, 1L]
        if(ni==(d1+1)){
          ix = order(abs(sparseM1[1:p, 3L]), decreasing = TRUE)
          j.I =sparseM1[ix,1][1:d1]
          lrows=p+(1:d1)
          sparseM1[lrows, 1L]=j.I
        }
        Xi <- x[, j.I, drop = FALSE]
        q <- length(j.I)
        
        if (q > 1L) {
          PPR <- try(ppr(Xi, Yi, nterms = 1), silent = TRUE)# sm.method="spline"
          #tryCatch({
          #  PPR <- ppr(Xi, Yi, nterms = 1)
          #}, error = function(e) {
          #  cat("ERROR:", e$message, "\n")
          #})
          
          if(inherits(PPR, "try-error")){
            LM = lm(Yi~.,data=data.frame(Xi))
            if(length(indC)>2){
              theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i))) {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
              theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            }else{
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
            gcv.i=sum(LM$residuals^2)
          }else{
            theta.i=PPR$alpha
            gcv.i=sum(PPR$residuals^2)
            #fitted=indC[apply(PPR$fitted.values,1,which.max)]
            #gcv.i = mean(fitted!=y)
          }
        }else{
          theta.i <- 1L
          gcv.i=Inf
        }
        #theta.i=theta.i/sqrt(sum(theta.i^2))
        sparseM1[lrows, 3L] <- c(theta.i)
        
        if(gcv.i<gcv){
          #theta=theta.i
          sparseM0=sparseM1[lrows,]
          #sparseM0[lrows, 3L] <- c(theta.i)
          gcv=gcv.i
        }
      }
      #sparseM=rbind(sparseM1,sparseM[-(1:p),])
      sparseM=rbind(sparseM0,sparseM[-(1:p),])
      sparseM[,2]=c(rep(1,length(sparseM0[,2])),2:(d3+1))
      
      #sparseM=sparseM0
      #sparseM[,2]=1
      #sparseM=sparseM[-(1:p),]
      #sparseM[,2]=1:d3
    }else{
      #prob=0.5
      #nnzs <- round(p * d * sparsity)
      #ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
      #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      ind <- sample.int(p, d, replace = FALSE)
      sparseM = cbind(ind, 1:d, rep(1, d))
    }
    
    #options (warn = 1)
    return(sparseM)
  }
}



if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(sqrt(length(y)),
                                               ncol(x)/3)), sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),
                                                                              1/ncol(x)), catMap = NULL, ...) {
    options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    
    p = ncol(x)
    #if(1==2){
    if(d>=3){
      #d1 = round(d/3)
      #d2 = round(d/3)
      #d3 = d - d1-d2
      
      #d1=round(d/2)
      #d3 = d - d1
      d1=d
      d3=ceiling(sqrt(p))#round(p/3)
      
      #d1=d
      
      ind <- sample.int(p, replace = FALSE);
      s = c(0,round(quantile(1:p, (1:d1)/d1)))
      sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p))
      
      #prob=0.5
      #nnzs <- round(p * d2 * sparsity)
      #ind <- sort(sample.int((p * d2), nnzs, replace = FALSE))
      #sparseM2 <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L+d1, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      
      ind <- sample.int(p, d3, replace = FALSE)
      sparseM3 = cbind(ind, (1:d3) + d1, rep(1, d3))
      
      sparseM=rbind(sparseM1,sparseM3)#sparseM1#
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      #y <- as.factor(y)
      library(dummies)
      Y = as.matrix(dummy(as.factor(y)))
      y <- as.numeric(as.character(y))
      
      sparseM1=rbind(sparseM[1:p,],matrix(d1+1,d1,3))
      sparseM[-(1:p),2]=sparseM[-(1:p),2]+1
      for (ni in 1:(d1+1)) {
        lrows <- which(sparseM1[,2L] == ni)
        j.I <- sparseM1[lrows, 1L]
        if(ni==(d1+1)){
          ix = order(abs(sparseM1[1:p, 3L]), decreasing = TRUE)
          j.I =sparseM1[ix,1][1:d1]
          lrows=p+(1:d1)
          sparseM1[lrows, 1L]=j.I
        }
        Xi <- x[, j.I, drop = FALSE]
        q <- length(j.I)
        
        theta.i <- 1L
        if (q > 1L) {
          if(length(unique(y))>2){Yi=Y}else{Yi=y}
          theta.i <- try(ppr(Xi, Yi, nterms = 1)$alpha, silent = TRUE)
          
          if(inherits(theta.i, "try-error")){
            LM = lm(Yi~Xi)
            if(identical(Yi,Y)){
              theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i))) {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
              theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            }else{
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
          }
          
        }
        #theta.i=theta.i/sqrt(sum(theta.i^2))
        sparseM1[lrows, 3L] <- c(theta.i)
      }
      sparseM=rbind(sparseM1,sparseM[-(1:p),])#sparseM1#
    }else{
      #prob=0.5
      #nnzs <- round(p * d * sparsity)
      #ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
      #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind -1L)/p) + 1L, sample(c(1L, -1L),
      #                                                                          nnzs, replace = TRUE,prob = c(prob, 1 - prob)))
      ind <- sample.int(p, d, replace = FALSE)
      sparseM = cbind(ind, 1:d, rep(1, d))
    }
    
    options (warn = 1)
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = 2*ceiling(ncol(x)/min(sqrt(length(y)),
                                                         ncol(x)/3)), sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),
                                                                                        1/ncol(x)), catMap = NULL, ...) {
    options (warn = -1)
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    p = ncol(x)
    #y <- as.factor(y)
    library(dummies)
    Y = as.matrix(dummy(as.factor(y)))
    y <- as.numeric(as.character(y))
    
    #d2 = 1
    #  d1 = d#ceiling(2 * d/3)
    
    d1 = floor(d/2)
    d2 = d - d1
    
    
    #nnzs <- round(p * d1 * sparsity)
    #ind <- sort(sample.int((p * d1), nnzs, replace = FALSE))
    #sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind - 1L)/p) +
    #                  1L, rep(1, nnzs))
    
    nnzs=p
    ind <- sample.int(p, replace = FALSE);
    s = c(0,round(quantile(1:p, (1:d1)/d1)))
    sparseM <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, nnzs))
    
    if(d2>0){
      #ind <- sample.int(p, d2, replace = FALSE)
      #ind1=which(sparseM[,2]==1)
      #if(length(ind1)>0){
      #  ind1=which(ind%in%unique(sparseM[ind1,1]))
      #  if(length(ind1)>0){
      #  ind=ind[-ind1]
      #  }
      #}
      ind <- sample.int(p, d2, replace = FALSE)
      sparseM = rbind(sparseM, cbind(ind, (1:d2) + d1, rep(1, d2)))
      
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      sparseM1=rbind(sparseM[1:p,],matrix(d1+1,d1,3))
      sparseM[-(1:p),2]=sparseM[-(1:p),2]+1
      for (ni in 1:(d1+1)) {
        lrows <- which(sparseM1[,2L] == ni)
        j.I <- sparseM1[lrows, 1L]
        if(ni==(d1+1)){
          ix = order(abs(sparseM1[1:p, 3L]), decreasing = TRUE)
          j.I =sparseM1[ix,1][1:d1]
          lrows=p+(1:d1)
          sparseM1[lrows, 1L]=j.I
        }
        
        Xi <- x[, j.I, drop = FALSE]
        q <- length(j.I)
        
        theta.i <- 1L
        if (q > 1L) {
          if(length(unique(y))>2){Yi=Y}else{Yi=y}
          theta.i <- try(ppr(Xi, Yi, nterms = 1)$alpha, silent = TRUE)
          
          if(inherits(theta.i, "try-error")){
            LM = lm(Yi~Xi)
            if(identical(Yi,Y)){
              theta.i=as.matrix(LM$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i))) {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
              theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
            }else{
              theta.i=LM$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
          }
          
        }
        sparseM1[lrows, 3L] <- c(theta.i)
      }
      
      #sparseM=sparseM1
      sparseM=rbind(sparseM1,sparseM[-(1:p),])
    }
    
    options (warn = 1)
    return(sparseM)
  }
}

if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(length(y)^0.4, ncol(x)*2/3)),
                         sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)), catMap = NULL, ...) {
    p = ncol(x)
    n = nrow(x)
    y <- as.factor(y)
    y.class = unique(y)
    Y = matrix(0, n, length(y.class))
    for (ic in 1:length(y.class))
    {
      Y[y==y.class[ic],ic] = 1
    }
    
    d = ceiling(min(n^0.4, p*2/3))
    K = ceiling(p/d)
    sparseM <- matrix(0, (K+1)*d, 3)
    
    index = sample(1:p, p)
    
    for (k in 1:(K+1))
    {
      if (k <=K)
      {
        Ik = ((k-1)*d+1):min(k*d,p)
        Jk = index[Ik]
        #Jk = sample(1:p, d)
      }else
      {
        Ik = (p+1):(p+d)
        ix = order(abs(sparseM[K*d, 3]), decreasing = TRUE)
        Jk = sparseM[ix[1:min(length(ix), d)],1]
      }
      
      sparseM[Ik,1] = Jk
      sparseM[Ik,2] = k
      Xi = x[,Jk]
      if (length(Ik) > 1)
      {
        j.ppr <- try(ppr(Xi, Y, nterms = 1), silent = TRUE)
        if(inherits(j.ppr, "try-error"))
        {
          j.ppr = lm(Y~Xi)
          theta.i=as.matrix(j.ppr$coefficients)[-1,,drop = FALSE]
          for (j in seq(ncol(theta.i)))
          {
            theta.i[is.na(theta.i[,j]),j] = 0.0
          }
          theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
        }else{
          theta.i <- j.ppr$alpha
        }
      }else
        theta.i = 1
      
      sparseM[Ik, 3] <- c(theta.i)
    }
    
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(length(y)^0.4, ncol(x)*2/3)),
                         sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)),
                         catMap = NULL,method = c('sir', 'ppr','no')[2],ML=c("lm","glm")[1], ...) {
    p = ncol(x)
    n = nrow(x)
    #y <- as.factor(as.character(y))
    #factor(y, labels = c(0,1))
    y.class = unique(y)
    Y = matrix(0, n, length(y.class))
    for (ic in 1:length(y.class))
    {
      Y[y==y.class[ic],ic] = 1
    }
    
    #d = ceiling(min(n^0.4, p*2/3))
    K = ceiling(p/d)
    sparseM <- matrix(0, p, 3)
    
    index = sample(1:p, p)
    for (k in 1:K)
    {
      Ik = ((k-1)*d+1):min(k*d,p)
      Jk = index[Ik]
      sparseM[Ik,1] = Jk
      sparseM[Ik,2] = k
      X = x[,Jk]
      if (length(Ik) > 1)
      {
        if (method == 'sir'){
          theta.i <- try(dr(Y~X)$evectors[,1], silent = TRUE)
        }
        if (method == 'ppr'){
          theta.i <- try(ppr(X, Y, nterms = 1)$alpha, silent = TRUE)
        }
        if (method == "no"){
          theta.i=NULL
        }
        
        if(inherits(theta.i, "try-error")||(is.null(theta.i)))
        {
          if(length(unique(y))>2){
            if(ML=="lm"){
              j.ppr = lm(Y~X)
              theta.i=as.matrix(j.ppr$coefficients)[-1,,drop = FALSE]
              for (j in seq(ncol(theta.i)))
              {
                theta.i[is.na(theta.i[,j]),j] = 0.0
              }
            }
            if(ML=="glm"){
              theta.i=matrix(0,ncol(X),ncol(Y))
              for (jp in seq(ncol(Y))) {
                theta.i[,jp] =glm(Y[,jp]~X, family = "binomial")$coefficients[-1]
              }
            }
            theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
          }else{
            if(ML=="lm"){
              theta.i = lm(Y[,1]~X)$coefficients[-1]
              theta.i[is.na(theta.i)] = 0.0
            }
            if(ML=="glm"){
              theta.i =glm(Y[,1]~X, family = "binomial")$coefficients[-1]
            }
          }
          
        }else{
          if(length(theta.i)<ncol(X)){
            xj=which(!paste0("X",colnames(X))%in%names(theta.i))
            theta.i[xj]=0
            #BB=rbind(BB,matrix(0,p-nrow(BB),5))
          }
        }
      }else{
        theta.i = 1
      }
      
      theta.i=theta.i/sqrt(sum(theta.i^2))
      sparseM[Ik, 3] <- c(theta.i)
    }
    
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(length(y)^0.4, ncol(x)*2/3)),
                         sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),1/ncol(x)), catMap = NULL, ...) {
    p = ncol(x)
    n = nrow(x)
    y <- as.factor(y)
    y.class = unique(y)
    Y = matrix(0, n, length(y.class))
    for (ic in 1:length(y.class))
    {
      Y[y==y.class[ic],ic] = 1
    }
    
    d = ceiling(min(n^0.4, p*2/3))
    K = ceiling(p/d)
    sparseM <- matrix(0, p, 3)
    
    index = sample(1:p, p)
    for (k in 1:K)
    {
      Ik = ((k-1)*d+1):min(k*d,p)
      Jk = index[Ik]
      sparseM[Ik,1] = Jk
      sparseM[Ik,2] = k
      Xi = x[,Jk]
      if (length(Ik) > 1)
      {
        j.ppr <- try(ppr(Xi, Y, nterms = 1), silent = TRUE)
        if(inherits(j.ppr, "try-error"))
        {
          j.ppr = lm(Y~Xi)
          theta.i=as.matrix(j.ppr$coefficients)[-1,,drop = FALSE]
          for (j in seq(ncol(theta.i)))
          {
            theta.i[is.na(theta.i[,j]),j] = 0.0
          }
          theta.i = eigen(theta.i%*%t(theta.i))$vectors[,1]
        }else{
          theta.i <- j.ppr$alpha
        }
      }else
        theta.i = 1
      
      sparseM[Ik, 3] <- c(theta.i)
    }
    
    return(sparseM)
  }
}


if(1==2){
  RandMatPPR <- function(x, y, d = ceiling(min(sqrt(length(y)),
                                               ncol(x)/3)), sparsity = ifelse(ncol(x) >= 10, 3/ncol(x),
                                                                              1/ncol(x)), catMap = NULL, ...) {
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    p = ncol(x)
    # y <- as.numeric(as.character(y))
    y <- as.factor(y)
    
    d1 = ceiling(2 * d/3)
    # d1 = floor(d/2)
    d2 = d - d1
    
    nnzs <- round(p * d1 * sparsity)
    ind <- sort(sample.int((p * d1), nnzs, replace = FALSE))
    sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind - 1L)/p) +
                       1L, rep(1, d1))
    if(d2>0){
      # ind <- sample.int(p, replace = FALSE) s = c(0,
      # round(quantile(1:p, (1:d2)/d2))) sparseM <-
      # rbind(sparseM, cbind(ind, rep(1:d2, s[-1] - s[-d2 - 1])
      # + d1, rnorm(p)))
      
      #ind <- sample.int(p, d2, replace = FALSE)
      #ind1=which(sparseM[,2]==1)
      #if(length(ind1)>0){
      #  ind1=which(ind%in%unique(sparseM[ind1,1]))
      #  if(length(ind1)>0){
      #  ind=ind[-ind1]
      #  }
      #}
      ind <- sample.int(p, d2, replace = FALSE)
      sparseM = rbind(sparseM, cbind(ind, (1:d2) + d1, rep(1, d2)))
      
      if (!is.null(catMap)) {
        rw <- sparseM[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM[, 1L] <- rw
      }
      
      sparseM1=sparseM[-((1:d2) + d1),]
      for (ni in unique(sparseM1[,2L])) {
        lrows <- which(sparseM1[,2L] == ni)
        j.I <- sparseM1[lrows, 1L]
        Xi <- x[, j.I, drop = FALSE]
        q <- length(j.I)
        
        theta.i <- 1L
        if (q > 1L) {
          Y = as.matrix(dummy(y))
          j.ppr <- try(ppr(Xi, Y, nterms = 1))#, silent = TRUE
          if (inherits(j.ppr, "try-error")) {
            #sir = function(X,Y,d=1)
            #{
            # n = dim(X)[1]
            #D=as.matrix(dist(Y))
            #XC=(t(X)-apply(X,2,mean))
            #eig = eigen(cov(X))
            #inv = eig$vectors%*%diag( 1/(eig$values+1/n^2))%*%t(eig$vectors)
            #W=-XC%*%D%*%t(XC)/(n*(n-1))
            #return(svd(inv%*%W)$u[,1:d])
            #}
            colnames(Xi)=paste0("x",seq(q))
            j.ppr <- try(cancor(Xi, Y))
            if(inherits(j.ppr, "try-error")){
              next
            }else{
              theta.i=rep(0,q)
              theta.i[which(colnames(Xi)%in%j.ppr$xcoef[, 1])]=j.ppr$xcoef[, 1]
            }
            # j.ppr <- lm(y ~ Xi) theta.i <-
            # as.numeric(j.ppr$coefficients[-1])
            # theta.i[is.na(theta.i)] <- 0
          } else {
            theta.i <- j.ppr$alpha
          }
        }
        sparseM1[lrows, 3L] <- c(theta.i)
      }
      sparseM[-((1:d2) + d1),]=sparseM1
    }
    
    return(sparseM)
  }
  
}


if (1 == 2) {
  RandMatPPR <- function(x, y, p = ncol(x), d = ceiling(sqrt(p)),
                         sparsity = ifelse(p >= 10, 3/p, 1/p), prob = 0.5) {
    # d=ceiling(sqrt(ncolX))
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))
    dp <- d * prob
    y <- as.numeric(as.character(y))
    
    nnzs <- round(p * d * sparsity)
    ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
    
    ## Determine if categorical variables need to be taken
    ## into consideration sparseM <- cbind(((ind - 1L) %%
    ## p) + 1L, floor((ind - 1L)/p) +1L, sample(c(1L, -1L),
    ## nnzs, replace = TRUE, prob = c(prob,1 - prob)))
    ## sparseM[,3L] <-ppr(x[,sparseM[,1]],y,
    ## nterms=1L)$alpha
    sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind - 1L)/p) +
                       1L, rep(1, nnzs))
    
    for (ni in unique(sparseM[, 2L])) {
      lrows <- which(sparseM[, 2L] == ni)
      j.I <- sparseM[lrows, 1L]
      Xi <- x[, j.I, drop = FALSE]
      q <- length(j.I)
      
      if (q > 1L) {
        S <- eigen(cov(Xi))  # using PCA to handle thecolinearity
        Ii = which(S$values > 0.01)
        Vec <- diag(q)
        Ii <- 1:q
        if (sum(S$values) > 0) {
          Ii <- 1:min(which(cumsum(S$values)/sum(S$values) >
                              0.99))
          Vec <- S$vectors[, Ii, drop = FALSE]
        }
        
        theta.i <- Vec
        Xi <- as.matrix(Xi %*% Vec)
        q <- length(Ii)
        if (q > 1L) {
          j.ppr <- try(ppr(Xi, y, nterms = 1), silent = TRUE)
          if (inherits(j.ppr, "try-error")) {
            j.ppr <- lm(y ~ Xi)
            theta.i <- as.numeric(j.ppr$coefficients[-1])
            theta.i[is.na(theta.i)] <- 0
          } else {
            theta.i <- j.ppr$alpha
          }
          theta.i <- as.matrix(theta.i)
          theta.i <- Vec %*% theta.i
        }
        sparseM[lrows, 3L] <- theta.i
      }
    }
    
    sparsity = max(floor(p/3), 1)
    ind <- sample.int(p, sparsity, replace = FALSE)
    sparseM = rbind(sparseM, cbind(ind, max(sparseM[, 2]) +
                                     (1:sparsity), rep(1, sparsity)))
    
    return(sparseM)
  }
}


if (1 == 2) {
  RandMatPPR <- function(x, y, p = ncol(x), d = ceiling(min(sqrt(length(y)),
                                                            p/3)), sparsity = ifelse(p >= 10, 3/p, 1/p), prob = 0.5) {
    # d=max(2,floor(min(length(y)^0.4,ncol(x)/3)))#ceiling(sqrt(p))
    dp <- d * prob
    y <- as.numeric(as.character(y))
    
    d1 = floor(d/3)
    d2 = floor(d/3)
    d3 = d - d1 - d2
    nnzs <- round(p * d1 * sparsity)
    ind <- sort(sample.int((p * d1), nnzs, replace = FALSE))
    sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind - 1L)/p) +
                       1L, rnorm(nnzs))
    
    ind <- sample.int(p, replace = FALSE)
    s = c(0, round(quantile(1:p, (1:d2)/d2)))
    sparseM <- rbind(sparseM, cbind(ind, rep(1:d2, s[-1] -
                                               s[-d2 - 1]) + d1, rnorm(p)))
    
    for (ni in unique(sparseM[, 2L])) {
      lrows <- which(sparseM[, 2L] == ni)
      j.I <- sparseM[lrows, 1L]
      Xi <- x[, j.I, drop = FALSE]
      q <- length(j.I)
      
      theta.i <- 1L
      if (q > 1) {
        j.ppr <- try(ppr(Xi, dummy(y), nterms = 1), silent = TRUE)
        if (inherits(j.ppr, "try-error")) {
          j.ppr <- lm(y ~ Xi)
          theta.i <- as.numeric(j.ppr$coefficients[-1])
          theta.i[is.na(theta.i)] <- 0
        } else {
          theta.i <- j.ppr$alpha
        }
      }
      sparseM[lrows, 3L] <- theta.i
    }
    
    ind <- sample.int(p, d3, replace = FALSE)
    sparseM = rbind(sparseM, cbind(ind, (1:d3) + d1 + d2,
                                   rep(1, d3)))
    
    return(sparseM)
  }
}


if (1 == 2) {
  RandMatPPR <- function(x, y, p = ncol(x), d = max(2, floor(min(length(y)^0.4,
                                                                 ncol(x)/3))), sparsity = 0.3, prob = 0.5) {
    dp = d * prob
    y = as.numeric(as.character(y))
    nnzs <- round(p * d * sparsity)
    ind <- sort(sample.int((p * d), nnzs, replace = FALSE))
    
    sparseM <- cbind(((ind - 1L)%%p) + 1L, floor((ind - 1L)/p) +
                       1L, rnorm(nnzs))
    
    if (1 == 1) {
      nz.idx <- 1L
      nnz <- nrow(sparseM)
      # Check each projection to determine which splits the
      # best.
      while (nz.idx <= nnz) {
        
        feature.idx <- sparseM[nz.idx, 2L]
        feature.nnz <- 0L
        while (sparseM[nz.idx + feature.nnz, 2L] == feature.idx) {
          feature.nnz <- feature.nnz + 1L
          if (nz.idx + feature.nnz > nnz) {
            break
          }
        }
        # lrows are the elements in sparseM that will be
        # used to rotate the data
        lrows <- nz.idx:(nz.idx + feature.nnz - 1L)
        
        
        j.I = sparseM[lrows, 1]
        Xi = x[, j.I]
        q = length(j.I)
        
        if (q > 1) {
          if (1 == 1) {
            S = eigen(cov(Xi))  # using PCA to handle the colinearity
            # Ii = which(S$values > 0.01)
            Vec = diag(q)
            Ii = 1:q
            if (sum(S$values) > 0) {
              Ii = 1:min(which(cumsum(S$values)/sum(S$values) >
                                 0.99))
              Vec = as.matrix(S$vectors[, Ii])
            }
            
            Xi = Xi %*% Vec
            
            j.ppr = try(ppr(Xi, y, nterms = 1), silent = TRUE)
            if (inherits(j.ppr, "try-error")) {
              j.ppr = lm(y ~ Xi)
              theta.i = as.numeric(j.ppr$coefficients[-1])
              theta.i[is.na(theta.i)] = 0
            } else {
              theta.i = j.ppr$alpha
            }
            theta.i = as.matrix(theta.i)
            
            if (length(Ii) > 1) {
              theta.i = Vec %*% theta.i
            } else {
              theta.i = Vec
            }
          }
          theta.i = ppr(Xi, y, nterms = 1)$alpha
        } else {
          theta.i = 1L
        }
        
        sparseM[lrows, 3] = theta.i
        nz.idx <- nz.idx + feature.nnz
      }
    }
    
    return(sparseM)
  }
}

# roxygen2::roxygenise()
