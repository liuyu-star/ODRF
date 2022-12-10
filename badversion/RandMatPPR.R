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


#################################################
#library(R.matlab)

#load(paste0("Xy",seed,".RData"))
#load("Xy.RData")

#DataIndx=unlist(read.csv('DataIndx.csv',header = FALSE))
#DataIndx=c(read.csv(paste0('DataIndx',seed,'.csv')))
#DataIndx=c(readMat(paste0("DataIndx",seed,".mat"))$currentDataIndx)
#Xy=readMat("Xy.mat")
#X=as.matrix(Xy$X)
#y=c(Xy$y)


#sparseM=RandMatPPR(X[DataIndx,], y[DataIndx]);
#assign(paste0("PPMat",seed),sparseM)
#colnames(sparseM)=c("features","direction_number", "direction");
#q() n
#writeMat(paste0("PPMat",seed,".mat"), eval(parse(text = "PPMat1"))=sparseM)
#writeMat(paste0("PPMat",seed,".mat"), get(paste0("PPMat",seed))=sparseM)
#write.table(sparseM,file =paste0("PPMat",seed,".csv"),sep = ",",
#            append=FALSE,row.names = FALSE,col.names = FALSE)
#write.table(sparseM,file ="PPMat.csv",sep = ",",
#            append=FALSE,row.names = FALSE,col.names = FALSE)


if(1==2){
  # Check each projection to determine which splits the best.
  numDr=unique(sparseM[, 2L])
  Xpp=matrix(0,length(y),length(numDr))
  for (ndr in numDr) {
    lrows <- which(sparseM[, 2L] == ndr)
    
    # Project input into new space
    Xpp[,ndr] <- X[,sparseM[lrows, 1L], drop = FALSE] %*% sparseM[lrows, 3L, drop = FALSE]
  }
}

if(1==2){
  numDr=unique(sparseM[, 2L])
  PPMat=matrix(0,length(numDr),ncol(X))
  for (ndr in numDr) {
    lrows <- which(sparseM[, 2L] == ndr)
    PPMat[ndr,sparseM[lrows, 1L]] <- sparseM[lrows, 3L]
  }
  #Xpp=X%*%t(PPMat)
}


#writeMat("PPMat.mat", PPMat=sparseM)
#write.table(PPMat, file="PPMat.csv",sep = ",",append=FALSE,row.names = FALSE,col.names = TRUE)
#writeMat("Xpp.mat", Xpp=Xpp)
