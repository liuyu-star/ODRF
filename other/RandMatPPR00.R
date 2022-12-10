RandMatPPR01 <- function(x, y, q = min(ceiling(length(y)^0.4),ceiling(ncol(x)*2/3)),
                         catMap = NULL,method='c',...) {
  p = ncol(x)
  n=length(y)
  #d=min(100, max(5, ceiling(p/q)))
  d=max(5, ceiling(p/q))

  d1=d
  d2=d
  #ceiling(sqrt(p))#d#round(p/3)#
  
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
      
      if (!is.null(catMap)) {
        rw <- sparseM1[, 1]
        for (j in 1:length(catMap)) {
          isj <- rw %in% catMap[[j]]
          rw[isj] <- sample(catMap[[j]], length(rw[isj]), replace = TRUE)
        }
        sparseM1[, 1L] <- rw
      }
      
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
        }else{
          theta.i <- 1L
        }
        #theta.i=theta.i/sqrt(sum(theta.i^2))
        sparseM1[lrows, 3L] <- c(theta.i)
      }
      
      sparseM1[,2]=sparseM1[,2]+(k-1)*(d+1)
      sparseM=rbind(sparseM,sparseM1)
    }
    
    ind <- sample.int(p, d2, replace = FALSE)
    sparseM1 = cbind(ind, (1:d2) + k*(d+1), rep(1, d2))
    sparseM=rbind(sparseM,sparseM1)
    
  }else{
    d2=ceiling(sqrt(p))
    ind <- sample.int(p, d2, replace = FALSE)
    sparseM = cbind(ind, 1:d2, rep(1, d2))
  }

  return(sparseM)
}
