RandMatPPR01 <- function(x, y, q = min(ceiling(length(y)^0.4),ceiling(ncol(x)*2/3)),
                         catMap = NULL,method='c',...) {
  X = x
  p = ncol(x)
  n = length(y)
  sparsem = c()
  
  d3 = ceiling(sqrt(p))
  ind <- sample.int(p, d3, replace = FALSE)
  sparseM = cbind(ind, 1:d3, rep(1, d3))
  
  
  if(n > 20)  {
    d = min(ceiling(n^0.4), floor(p*2/3))
    n.search = max(5, floor(p/d))

    ind <- sample.int(p, replace = FALSE)
    
    icv = 1.0e10
    for (i in 1:n.search)
    {
      I = ((i-1)*d+1):(i*d)
      j.I = unique(ind)[1:d]
      ind = ind[(d+1):length(ind)]
      if (length(ind) < p)
        ind = c(ind, sample.int(p, replace = FALSE))
      
      Yi = c(y)
      indC = 0L
      if (method != "r") {
        y = as.factor(y)
        indC = levels(y)
        if (length(indC) > 2) {
          Yi = (matrix(y, n, length(indC)) == matrix(indC, 
                            n, length(indC), byrow = TRUE)) + 0
        }
        else {
          Yi = as.integer(y)
        }
      }
      Xi = X[,j.I]
      
      S = eigen(cov(Xi))                # using PCA to handle the colinearity
      Ii = 1:min(which(cumsum(S$values)/sum(S$values) > 0.99))
      Xi = Xi%*%S$vectors[,Ii]
      
      j.ppr = try(ppr(Xi, Yi,  nterms=1))
      if(inherits(j.ppr, "try-error"))
      {
        j.ppr = lm(Yi~Xi)
        theta.i = as.numeric(j.ppr$coefficients[2:(length(Ii)+1)])
        theta.i[is.na(theta.i)] = 0.0
        theta.i = S$vectors[,Ii]%*%theta.i
      }else
      {
        if (length(Ii) > 1)
          theta.i = S$vectors[,Ii]%*%j.ppr$alpha
        else
          theta.i = S$vectors[,Ii]
      }
      
      j.icv = mean(j.ppr$residuals^2) # /(1-2*q/n)^2
      if (j.icv < icv)
      {
        I = j.I
        theta = theta.i
        icv = j.icv
      }
      
      sparsem = cbind(I, d3+1, theta)
    }
  }
    
  
  sparseM = rbind(sparseM, sparsem)
    
  return(sparseM)
}
