#include <Rcpp.h>
#include "node_cuts.h"
using namespace Rcpp;

// [[Rcpp::export]]
List best_cut_node(char method, double lambda, NumericMatrix Data, NumericVector Labels, NumericVector W,int minleaf, int numLabels) {
  
  int M=Data.nrow();
  int N=Data.ncol();
  //int *bcvar = new int;
  //double *bcval = new double;
  int bcvar=-1;
  double bcval=0;
  double *bestval = new double[N];
  
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
    case 'g':
      gini_split(lambda, M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'i':
      infogain_split(M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'r':
      mse_split(lambda, M, N, Labels1, Data1, minleaf, &bcvar, &bcval, bestval);
      break;
    }
  }else{
    double *W1 = new double[M];
    for(int k=0;k<M;k++){
      W1[k]=W[k];
    }
    
    switch (method){
    case 'g':
      gini_split(lambda, M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'i':
      infogain_split(M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'r':
      mse_split(lambda, M, N, Labels1, Data1, minleaf, &bcvar, &bcval, bestval);
      break;
    }
    
    delete[] W1;
  }
  
  NumericVector BestVal(N);
  for(int k=0;k<N;k++){
    BestVal[k]=bestval[k];
  }
  
  List bestCut = List::create(
    _["BestCutVar"]= bcvar, 
    _["BestCutVal"]= bcval,
    _["BestIndex"]= BestVal);
  
  delete[] Data1;
  delete[] Labels1;
  delete[] bestval;
  
  return(bestCut);
}

