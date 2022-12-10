#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;

NumericVector fastdcov2d(NumericMatrix Y,NumericVector y,int n,int p,
                         NumericVector index,NumericVector a_x,Function dcovsq){
  NumericVector covsq={as<double>(dcovsq(X(_,0),y,n,index,a_x))};
  for (int j=1; j!=p; ++j) {
    covsq.push_back(as<double>(dcovsq(X(_,j),y,n,index,a_x)));
  }
  return covsq;
}

NumericVector FastDcov2d=function(X,y){
//fastDcov computes distance correlation between column matrix X and y
  library(Rcpp)
  

  
  Y=data.matrix(X); x=c(y); n=length(x); p=ncol(Y)
//Sort the data by x. This does not change the answer
//a_x is the vector of row sums of distance matrix of x
    index=order(x); x=x[index]
    si = cumsum(x); s = si[n] 
    a_x = (2*seq(n)-n)*x+(s-2*si)
      
 NumericVector fastdcov2d(NumericMatrix Y,NumericVector x,int n,int p,
                                             NumericVector index,NumericVector a_x,Function dcovsq){
                      NumericVector covsq={as<double>(dcovsq(Y(_,0),x,n,index,a_x))};
   NumericVector::iterator it;
   for(it = x.begin() + 1; it != x.end(); ++it) {
                      for (int j=1; j!=p; ++j) {
                        covsq.push_back(as<double>(dcovsq(Y(_,j),x,n,index,a_x)));
                      }
                      return covsq;
                    }
      
      return(fastdcovsq(Y,x,n,p,index,a_x,DCovSq))
}



  
  double *xy, *y;
  double ah, bh, ch, sum_all, sum_all2, suml, suml2, sumr, sumr2, sr, sl, cl;
  int i, j, mj;
  
  y = new double[n];
  xy =  new double[n];

  
  
    
VectorXd DCov2d(VectorXd z,VectorXd x,VectorXd index, VectorXd a_x){
  int n=x.size(); 
  MatrixXd it(n,4);
  VectorXd y(n), xy(n);
  
  y=z(index);   
  xy=x.dot(y);
  //it=matrix(cbind(it1=rep(1,n),it2=y,it3=x,it4=xy),ncol = 4)
    it.col(1)=MatrixXd::Ones(n,1)
    it.col(2)=y;
    it.col(3)=x;
    it.col(4)=xy;
    
  // id1 , id2 , id3 , id4 are used to store sum of weights. 
  // On output , for 1 ≤ j ≤ n
  MatrixXd id=MatrixXd::Zero(n,4);
  //colnames(id)=paste0("id",seq(4))
    id[2:n,]=matrix(sapply(2:n,function(i){colSums(matrix(it[which(y[seq(i-1)]<y[i]),],ncol = 4))})
                      ,ncol = 4,byrow = T)
    
    // b_y is the vector of row sums of distance matrix of y
    idy=order(y); sy=y[idy]
    si = cumsum(sy); s = si[n] ;b_y=rep(0,n)
      b_y[idy] = (2*seq(n)-n)*sy+(s-2*si)
      
      // the Frobenius inner product of the distance matrices 
      covterm = as.numeric(n*((x-mean(x))%*%(y-mean(y))))
      c1 = as.numeric(xy%*%id[,1]) 
      c2 = as.numeric(x%*%id[,2])
      c3 = as.numeric(y%*%id[,3])
      c4 = sum(id[,4])
      
      D = 4*(c1-c2-c3+c4)-2*covterm
      
      // covsq equal is the square of the distance covariance between x and y 
    term1 = D/n^2 
    term2 = as.numeric(2*(a_x%*%b_y))/n^3 
    term3 = (sum(a_x)*sum(b_y))/n^4
    
    covsq = term1-term2+term3 
//covsq=apply(X,2,Cal_DCor)
             return(as.numeric(covsq))
}

void DCov2d(double* y,double* x, int n, double* index, double* a_x, double* covsq){
  
  double *sorted_data, *sorted_labels;
  double ah, bh, ch, sum_all, sum_all2, suml, suml2, sumr, sumr2, sr, sl, cl;
  int i, j, mj;
  
  sorted_labels = new double[M];
  sorted_data =  new double[M];
  
  sum_all = 0;
  sum_all2 = 0;
  for(j=0;j<M;j++){
    sum_all+=Labels[j];
    sum_all2 += Labels[j]*Labels[j];
  }
  
  //bh = sum_all2 + sum_all*(sum_all/M) - 2*(sum_all/M)*sum_all;
  bh = (sum_all2 - sum_all*(sum_all/M))/M;
  
  for(i = 0;i<N;i++){
    ah=1e+10;
    
    suml = 0;
    suml2 = 0;
    sumr = sum_all;
    sumr2 = sum_all2;
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
    }
    
    quicksort(sorted_data, sorted_labels, 0, M-1);
    
    for(mj = 0 ; mj<minleaf;mj++){
      cl=sorted_labels[mj];
      suml += cl;
      sumr -= cl;
      suml2 += cl*cl;
      sumr2 -= cl*cl;
    }
    
    for(j = minleaf;j<M-minleaf;j++){
      
      cl=sorted_labels[j];
      
      suml += cl;
      sumr -= cl;
      
      suml2 += cl*cl;
      sumr2 -= cl*cl;
      
      //sl = suml2 + suml*(suml/(j+1)) - 2*(suml/(j+1))*suml;
      //sr = sumr2 + sumr*(sumr/(M-j-1)) - 2*(sumr/(M-j-1))*sumr;
      
      sl = suml2 - suml*(suml/(j+1));
      sr = sumr2 - sumr*(sumr/(M-j-1));
      ch = (sl + sr)/M;
      
      //sl = (suml2 - suml*(suml/(j+1)))/(j+1);
      //sr = (sumr2 - sumr*(sumr/(M-j-1)))/(M-j-1);
      //ch = ((j+1)*sl +(M-j-1)*sr)/M;
      //ch = sl + sr;
      // ch = sl + sr + log(M)/M;
      
      if (ch<bh){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          bh=ch;
          *bcvar = i+1;
          *bcval = 0.5*(sorted_data[j+1]+sorted_data[j]);
        }
      }
      
      if (ch<ah){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          ah=ch;
        }
      }
    }
    
    bestval[i]= ah;
  }
  
  delete[] sorted_labels;
  delete[] sorted_data;
}