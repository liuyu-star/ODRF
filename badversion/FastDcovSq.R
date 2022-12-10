FastDcovSq=function(X,y){
  #fastDcov computes distance correlation between column matrix X and y
  library(Rcpp)
  
  DCovSq=function(y,x,n,index,a_x){
    
    y=y[index];   xy=x*y
    it=matrix(cbind(it1=rep(1,n),it2=y,it3=x,it4=xy),ncol = 4)
    
    # id1 , id2 , id3 , id4 are used to store sum of weights. 
    # On output , for 1 ≤ j ≤ n
    id=matrix(0,nrow = n,ncol =  4); colnames(id)=paste0("id",seq(4))
    id[2:n,]=matrix(sapply(2:n,function(i){colSums(matrix(it[which(y[seq(i-1)]<y[i]),],ncol = 4))})
                    ,ncol = 4,byrow = T)
    
    # b_y is the vector of row sums of distance matrix of y
    idy=order(y); sy=y[idy]
    si = cumsum(sy); s = si[n] ;b_y=rep(0,n)
    b_y[idy] = (2*seq(n)-n)*sy+(s-2*si)
    
    # the Frobenius inner product of the distance matrices 
    covterm = as.numeric(n*((x-mean(x))%*%(y-mean(y))))
    c1 = as.numeric(xy%*%id[,1]) 
    c2 = as.numeric(x%*%id[,2])
    c3 = as.numeric(y%*%id[,3])
    c4 = sum(id[,4])
    
    D = 4*(c1-c2-c3+c4)-2*covterm
    
    #covsq equal is the square of the distance covariance between x and y 
    term1 = D/n^2 
    term2 = as.numeric(2*(a_x%*%b_y))/n^3 
    term3 = (sum(a_x)*sum(b_y))/n^4
    
    covsq = term1-term2+term3 
    #covsq=apply(X,2,Cal_DCor)
    return(as.numeric(covsq))
  }
  
  Y=data.matrix(X); x=c(y); n=length(x); p=ncol(Y)
  #Sort the data by x. This does not change the answer
  #a_x is the vector of row sums of distance matrix of x
  index=order(x); x=x[index]
  si = cumsum(x); s = si[n] 
  a_x = (2*seq(n)-n)*x+(s-2*si)
  
  cppFunction(code = "
  NumericVector fastdcovsq(NumericMatrix Y,NumericVector x,int n,int p,
                NumericVector index,NumericVector a_x,Function dcovsq){
  NumericVector covsq={as<double>(dcovsq(Y(_,0),x,n,index,a_x))};
  for (int j=1; j!=p; ++j) {
  covsq.push_back(as<double>(dcovsq(Y(_,j),x,n,index,a_x)));
  }
  return covsq;
  }")
 
  return(fastdcovsq(Y,x,n,p,index,a_x,DCovSq))
}
