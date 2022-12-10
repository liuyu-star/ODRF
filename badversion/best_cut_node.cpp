#include <Rcpp.h>
#include "node_cuts.h"
//#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
//List best_cut_node(char* method, double* Data, double* Labels, double* W, int minleaf, int num_labels, double* bcvar, double* bcval) {
//List best_cut_node(const char method, NumericMatrix Data, NumericVector Labels, NumericVector W,const int minleaf, int const num_labels, double bcvar, double bcval) {
List best_cut_node(char method, NumericMatrix Data, NumericVector Labels, NumericVector W,int minleaf, int num_labels) {
  //CharacterVector
  //double *Data, *Labels, *W, *bcvar, *bcval;
  int N, M;// num_labe
  //method = (char*) malloc(2*sizeof(char));
  double *bcvar, *bcval;
  
  N=Data.ncol();
  M=Data.nrow();
  

  if (W.size()==0){
    switch (method){
    case 'c':
      GBCC(M, N, Labels, Data, minleaf, num_labels, bcvar, bcval);
      break;
    case 'g':
      GBCP(M, N, Labels, Data, minleaf, num_labels, bcvar, bcval);
      break;
    case 'r':
      GBCR(M, N, Labels, Data, minleaf, bcvar, bcval);
      break;
    }
  }
  else{
    switch (method){
    case 'c':
      GBCC(M, N, Labels, Data, W, minleaf, num_labels, bcvar, bcval);
      break;
    case 'g':
      GBCP(M, N, Labels, Data, W, minleaf, num_labels, bcvar, bcval);
      break;
    case 'r':
      GBCR(M, N, Labels, Data, W, minleaf, bcvar, bcval);
      break;
    }
  }
  
  //delete[] method;
  
  List bestCut = List::create(
    _["BestCutVar"]= bcvar, 
    _["BestCutVal"]= bcval) ;
  
  return(bestCut);
}

