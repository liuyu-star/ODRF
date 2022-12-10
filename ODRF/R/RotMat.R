#' Create a Random Matrix: Binary
#'
#'
#' @param p the number of dimensions.
#' @param d the number of desired columns in the projection matrix.
#' @param sparsity a real number in \eqn{(0,1)} that specifies the distribution of non-zero elements in the random matrix. sparsity="dist"
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
RotMatRand <- function(dimX,RandDist=c("Binary","Norm","Uniform")[1],numProj=ceiling(sqrt(paramList$dimX)),
                       dimProj = NULL,sparsity=ifelse(paramList$dimX >=10, 3/paramList$dimX, 1/paramList$dimX),
                       prob=0.5, lambda=1, catLabel = NULL, ...) {
  if(!is.null(dimProj)){
    if (dimProj > dimX) {
      stop("ERROR: parameter dimProj is greater than the number of dimensions dimX.")
    }
    
    nnz <- dimProj * numProj
    nz.rows <- integer(nnz)
    nz.cols <- integer(nnz)
    start.idx <- 1L
    for (i in seq.int(numProj)) {
      end.idx <- start.idx + dimProj - 1L
      nz.rows[start.idx:end.idx] <- sample.int(dimProj, dimProj, replace = FALSE)
      nz.cols[start.idx:end.idx] <- i
      start.idx <- end.idx + 1L
    }
  }else if(sparsity=="dist"){
    nnzPerCol <- stats::rpois(numProj, lambda)
    while (!any(nnzPerCol)) {
      nnzPerCol <- stats::rpois(numProj, lambda)
    }
    
    nnzPerCol[nnzPerCol > dimX] <- dimX
    nnz <- sum(nnzPerCol)
    nz.rows <- integer(nnz)
    nz.cols <- integer(nnz)
    start.idx <- 1L
    for (i in seq.int(numProj)) {
      if (nnzPerCol[i] != 0L) {
        end.idx <- start.idx + nnzPerCol[i] - 1L
        nz.rows[start.idx:end.idx] <- sample.int(dimX, nnzPerCol[i],replace = FALSE)
        nz.cols[start.idx:end.idx] <- i
        start.idx <- end.idx + 1L
      }
    }
  }else{
    nnzs <- round(dimX * numProj * sparsity)
    ind <- sort(sample.int((dimX * numProj), nnzs, replace = FALSE))
    nz.rows=((ind - 1L)%%dimX) + 1L
    nz.cols=floor((ind -1L)/dimX) + 1L
  }
  
  
  if(RandDist=="Binary"){
    proj=sample(c(-1L,1L),length(nz.rows), replace = TRUE)
  }
  if(RandDist=="Norm"){
    proj=rnorm(length(nz.rows))
  }
  if(RandDist=="Uniform"){
    proj=runif(length(nz.rows), -1,1)
  }
  
  randomMatrix <- cbind(nz.rows, nz.cols, proj)

  ## Determine if categorical variables need to be taken d
  ## into consideration
  if (!is.null(catLabel)) {
    ind=randomMatrix[,1]; 
    catVar=which(ind<=length(catLabel))
    #catVar=which(ind<=length(catLabel))
    randomMatrix[-catVar, 1L]=randomMatrix[-catVar, 1L]+length(unlist(catLabel))-length(catLabel)
    
    catMap=1
    for (xj in unique(ind[catVar])) {
      isj=which(ind==xj)
      randomMatrix[isj, 1L]<- sample(catMap:(catMap+length(catLabel[[xj]])-1),length(isj),replace = TRUE)
      catMap=catMap+length(catLabel[[xj]])
    }
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
RotMatRF <- function(dimX, numProj, catLabel = NULL, ...) {
  if (numProj > dimX) {
    stop("ERROR: parameter numProj is greater than the number of dimensions dimX.")
  }
  
  randomMatrix <- cbind(sample.int(dimX,numProj, replace = FALSE),1:numProj, rep(1L, numProj))
  if (!is.null(catLabel)) {
    ind=randomMatrix[,1]; 
    catVar=which(ind<=length(catLabel))
    #catVar=which(ind<=length(catLabel))
    randomMatrix[-catVar, 1L]=randomMatrix[-catVar, 1L]+length(unlist(catLabel))-length(catLabel)
    
    catMap=1
    for (xj in unique(ind[catVar])) {
      isj=which(ind==xj)
      randomMatrix[isj, 1L]<- sample(catMap:(catMap+length(catLabel[[xj]])-1),length(isj),replace = TRUE)
      catMap=catMap+length(catLabel[[xj]])
    }
  }
  
  return(randomMatrix)
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
#' @importFrom nnet nnet
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
RotMatPPO <- function(x, y, numProj,dimProj,catLabel = NULL,ppMethod="PPR",
                      method='i-classification',weights=NULL) {
  p = ncol(x)
  n=length(y)
  weights=ifelse(is.null(weights),rep(1, n),weights)
  #d=min(100, max(5, ceiling(p/q))) d q
  p0=p-(!is.null(catLabel))*(length(unlist(catLabel))-length(catLabel))
  
  d2=ceiling(sqrt(p0))
  ind <- sample(p0, d2, replace = FALSE)
  sparseM = cbind(ind, 1:d2, rep(1, d2))
  
  #ceiling(sqrt(p))#d#round(p/3)#
  #if((n>Q)&(n>10)&(p>1)){
  sparseM1=NULL
  if(n>10){
    if(identical(dimProj,"RAND")){
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
    for (xj in unique(ind[catVar])) {
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
          PP <- try(ppr(Xi, Yi,weights,nterms = 1, bass=1)$alpha, silent = TRUE)# sm.method="spline",sm.method="gcvspline"
        }
        if(ppMethod=="RAND"){
          PP <- sample(c(1L, -1L),ncol(Xi), replace = TRUE,prob = c(0.5, 0.5))
        }
        
        if(ppMethod=="NN"){
          #if((n>5*pi)&(pi<10)){
          #  PP <- try(ppr(Xi, Yi, nterms = 1, bass=1)$alpha, silent = TRUE)
          #}else{
            #PP = try(nnet(Xi, y, size=1,trace=FALSE)$wts[2:(1+pi)], silent = TRUE)
            PP = nnet(Xi, Yi,weights, size=1,linout=TRUE,trace=FALSE)$wts[2:(1+pi)]#  
          #}
        }
        
        if(!ppMethod%in%c("PPR","RAND","NN")){
          PP <- ppRF:::ppOptCpp(y,Xi,q = 1L, PPmethod = ppMethod, weight = TRUE, r = 1L, lambda = 0.1, 
                                    energy = 0, cooling = 0.9, TOL = 0.001, maxiter = 1000L)$projbest
        }

        if(inherits(PP, "try-error")){
          LM = lm(Yi~.,data=data.frame(Xi),weights)
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

##################################################################################
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
#' @export
#'
#' @examples
#'
#' p <- 8
#' d <- 3
#' paramList <- list(p = p, d = d)
#' set.seed(4)
#' (a <- do.call(RandMatRF, paramList))
RotMatMake <- function(x, y, RotMatFun, PPFun, FunDir, paramList) {
  if(!RotMatFun%in%ls("package:ODRF")){
    source(paste0(FunDir,"/",RotMatFun,".R"))
  }
  if(!PPFun%in%ls("package:ODRF")){
    source(paste0(FunDir,"/",PPFun,".R"))
  }
  
  RotMatFun <- match.fun(RotMatFun, descend = TRUE)
  PPFun <- match.fun(PPFun, descend = TRUE)
  
  paramList$x=x;paramList$y=y
  sparseM <- do.call(RotMatFun, paramList)
  for (ni in unique(sparseM[,2L])) {
    lrows <- which(sparseM[,2L] == ni)
    paramList$x=x[,sparseM[lrows, 1L],drop = FALSE]
    sparseM[lrows, 3L]=do.call(PPFun, paramList)
  }
  
  paramList$x = NULL;paramList$y = NULL;
  return(sparseM)
}

# roxygen2::roxygenise()
