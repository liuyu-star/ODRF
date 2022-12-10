#include <Rcpp.h>
#include "node_cuts.h"
using namespace Rcpp;

// [[Rcpp::export]]
List best_cut_node(char method, NumericMatrix Data, NumericVector Labels, NumericVector W,int minleaf, int num_labels) {

  int M=Data.nrow();
  int N=Data.ncol();
  //int *bcvar = new int;
  //double *bcval = new double;
  int bcvar=-1;
  double bcval=0;
  double *minval = new double[N];
  
  //std::vector<double> Data1(Data.size());
  //std::vector<double> Labels1(M);
  //double Data1[Data.size()];
  //double Labels1[M];
  //double* Data1 = (double*) malloc(Data.size()*sizeof(double));
  //double* Labels1 = (double*) malloc(M*sizeof(double));
  double *Data1 = new double[M*N];
  double *Labels1 = new double[M];
  
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++){
      Data1[i+j*M]=Data(i,j);
    }
    Labels1[i]=Labels[i];
  }
  
 
 if (W.size()==1){
   switch (method){
   case 'c':
     GBCC(M, N, Labels1, Data1, minleaf, num_labels, &bcvar, &bcval, minval);
     break;
   case 'g':
     GBCP(M, N, Labels1, Data1, minleaf, num_labels, &bcvar, &bcval, minval);
     break;
   case 'r':
     GBCR(M, N, Labels1, Data1, minleaf, &bcvar, &bcval, minval);
     break;
   }
   
 }else{
   double *W1 = new double[M];
   for(int k=0;k<M;k++){
     W1[k]=W[k];
   }
   
   switch (method){
   case 'c':
     GBCC(M, N, Labels1, Data1, W1, minleaf, num_labels, &bcvar, &bcval, minval);
     break;
   case 'g':
     GBCP(M, N, Labels1, Data1, W1, minleaf, num_labels, &bcvar, &bcval, minval);
     break;
   case 'r':
     GBCR(M, N, Labels1, Data1, W1, minleaf, &bcvar, &bcval, minval);
     break;
   }
   
   delete[] W1;
 }
 
 NumericVector MinVal(N);
 for(int k=0;k<N;k++){
   MinVal[k]=minval[k];
 }

  List bestCut = List::create(
    _["BestCutVar"]= bcvar, 
    _["BestCutVal"]= bcval,
    _["MinVal"]= MinVal);
  
  delete[] Data1;
  delete[] Labels1;
  delete[] minval;
  
  return(bestCut);
}

