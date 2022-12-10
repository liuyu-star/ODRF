// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
//using Eigen::Map;
//using Eigen::MatrixXd;
//using Rcpp::as;


// [[Rcpp::export]]
MatrixXd Test(NumericMatrix AA) {
  
  Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  MatrixXd BB;
  NumericMatrix B(as<NumericMatrix >(BB));
  return(B);
}
Rcpp::cumsum()
// [[Rcpp::export]]
List best_cut_node(NumericMatrix Data, NumericVector Labels) {
  
  //bool conv(NumericVector X) {  
   //   Eigen::Map<Eigen::VectorXd> 
   //   XS(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(X)); return true; }
  

  //Rcpp::sapply() Rcpp::colSums() Rcpp::seq_len()
  //NumericMatrix yi=it(y[seq(i-1)]<y[i],_);
  int n=3;
  NumericMatrix ity(n,4);
  ity.column(1)=Rcpp::
    MatrixXd it;
    NumericMatrix ity(as<NumericMatrix >(it));
    yi=ity(y[seq(i-1)]<y[i]);
    ity= sapply(2:n,function(i){colSums(ity(y[seq(i-1)]<y[i]))})
      for (i = 1; i < n; i++)
        id[i,]=colSums(ity(y[seq(i-1)]<y[i]))
      }
    
 id[2:n,]=matrix(,ncol = 4,byrow = T)
  return(bestCut);
}

