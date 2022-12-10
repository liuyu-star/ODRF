RandMatPPR03 <- function(x, y, d=NULL,catMap = NULL,method='c',...) {
  p = ncol(x)
  n=length(y)
  #d=min(100, max(5, ceiling(p/q)))
  
  p0=p-(!is.null(catMap))*(length(unlist(catMap))-length(catMap))
  if(is.null(d)){
    q =min(ceiling(n^0.4),ceiling(p0*2/3))
    d=min(p0,max(5, ceiling(p0/q)))
  }
  d1=d
  
  d2=ceiling(sqrt(p0))
  ind <- sample.int(p0, d2, replace = FALSE)
  sparseM = cbind(ind, 1:d2, rep(1, d2))
  if (!is.null(catMap)) {
    sparseM[which(ind>length(catMap)), 1L]=sparseM[which(ind>length(catMap)), 1L]+
      length(unlist(catMap))-length(catMap)
    for (j in 1:length(catMap)) {
      isj=which(ind==j)
      sparseM[isj, 1L]<- sample(catMap[[j]],length(isj),replace = FALSE)
    }
  }
  
  #ceiling(sqrt(p))#d#round(p/3)#
  #if((n>Q)&(n>10)&(p>1)){
  if(n>10){
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
    
    #sparseM=NULL
    #for (k in 1:ceiling(2*q*d/p)) {
    ind <- sample.int(p0, replace = FALSE);
    s = c(0,round(quantile(1:p0, (1:d1)/d1)))
    sparseM1 <-cbind(ind, rep(1:d1, s[-1] - s[-d1 - 1]), rep(1, p0))
    
    if (!is.null(catMap)) {
      sparseM1[which(ind>length(catMap)), 1L]=sparseM1[which(ind>length(catMap)), 1L]+
        length(unlist(catMap))-length(catMap)
      for (j in 1:length(catMap)) {
        isj=which(ind==j)
        sparseM1[isj, 1L]<- sample(catMap[[j]],length(isj),replace = FALSE)
      }
    }
    
    gcv=Inf
    indTab=table(sparseM1[,2])
    #d11=max(ind)
    ind=as.numeric(names(indTab)[which(indTab>1)])
    jx=which(sparseM1[,2]%in%ind)
    #d11=min(c(length(unique(sparseM[jx,1])),ifelse(method=='r',Inf,Inf)))#
    d11=max(max(indTab),min(length(unique(sparseM1[jx,1])),d1)) 
    
    sparseM1=rbind(sparseM1,matrix(d1+1,d11,3))
    for (ni in unique(sparseM1[,2L])) {
      lrows <- which(sparseM1[,2L] == ni)
      j.I <- sparseM1[lrows, 1L]
      if(ni==(d1+1)){
        #jx=which(sparseM1[, 3L]!=1L)
        ix = order(abs(sparseM1[jx,3L]), decreasing = TRUE)
        j.I =unique(sparseM1[jx[ix],1])[1:d11]
        lrows=p0+(1:d11)
        sparseM1[lrows, 1L]=j.I
      }
      Xi <- x[, j.I, drop = FALSE]
      #q <- length(j.I)
      
      S = eigen(cov(Xi))                # using PCA to handle the colinearity
      Ii = 1:min(which((cumsum(S$values)+1e-4)/sum(S$values,1e-4) > 0.99))
      Xi = Xi%*%S$vectors[,Ii]
      
      if (length(Ii) > 1L) {
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
          
          res=LM$residuals
        }else{
          theta.i=PPR$alpha
          res=PPR$residuals
        }
        
        theta.i=S$vectors[,Ii]%*%theta.i
        gcv.i=sum(res^2)
      }else{
        theta.i <- rep(1,length(j.I))
        gcv.i=Inf
      }
      theta.i=theta.i/sqrt(sum(theta.i^2))
      sparseM1[lrows, 3L] <- c(theta.i)
      
      if(gcv.i<=gcv){
        #theta=theta.i
        sparseM0=sparseM1[lrows,,drop = FALSE]
        #sparseM0[lrows, 3L] <- c(theta.i)
        gcv=gcv.i
      }
    }
    
    
    #sparseM1[,2]=sparseM1[,2]+(k-1)*(d+1)
    #sparseM=rbind(sparseM,sparseM1)
    
    sparseM0[,2]=d2+1
    sparseM=rbind(sparseM,sparseM0)
    #}
    
    #ind <- sample.int(p0, d2, replace = FALSE)
    #sparseM1 = cbind(ind, (1:d2) + k*(d+1), rep(1, d2))
    #sparseM1 = cbind(ind, (1:d2) + max(sparseM[,2]), rep(1, d2))
    #sparseM1 = cbind(ind, (1:d2) +1, rep(1, d2))
    #sparseM=rbind(sparseM,sparseM2)
    #sparseM[,2]=1:(d2+)
  }
  
  return(sparseM)
}
