#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List best_cut_node(NumericMatrix Data, NumericVector Labels, NumericVector W,int minleaf, int num_labels) {
  int N, M;
  double *Data1[Data.size()];
  double *Labels1[Labels.size()];
  int  *bcval;
  int i=1,j=1;
  int *bcvar;
  

  M=Data.nrow();
  N=Data.ncol();
  
  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      Data1[i+j*M]=&Data(i,j);
    }//best_cut_node(matrix(seq(100),20,5),seq(20),c(0,1),10,4)
    Labels1[i]=&Labels[i];
  }

  
  //bcvar[0] = i+1;
  bcval=&j;
  bcvar=&i;
 
  //if (W.size()==0){
   // bcvar=&M;
   // bcval=&N;
  //}
 
  List bestCut = List::create(
    _["BestCutVar"]= *bcvar,
    _["BestCutVal"]= *bcval,
    Named("data") = Labels  ,
    Named("MN") = M*N);//Labels1[1]+Data1[1]);
  return(bestCut);
}

