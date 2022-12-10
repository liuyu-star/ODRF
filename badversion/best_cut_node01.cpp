#include <Rcpp.h>
//#include "node_cuts.h"
//#include "quickie.h"
//#include <Rmath.h>
#include <cmath>
using namespace Rcpp;
//using namespace std; 

void GBCC(int , int , double* , double* , int , int , int* , double*);

// [[Rcpp::export]]
List best_cut_node(NumericMatrix Data, NumericVector Labels, NumericVector W,int minleaf, int num_labels) {

  int N, M, np=Data.size();
    
  M=Data.nrow();
  N=Data.ncol();
    
  double *bcval = new double;
  int *bcvar = new int;
  double Data1[np];
  double Labels1[M];

  
  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++){
      Data1[i+j*M]=Data(i,j);
    }
    Labels1[i]=Labels[i];
  }
  
 
  if (W.size()!=0){
      GBCC(M, N, Labels1, Data1, minleaf, num_labels, bcvar, bcval);
  }
 

  List bestCut = List::create(
    _["BestCutVar"]= *bcvar, 
    _["BestCutVal"]= *bcval) ;
  
  return(bestCut);
}


template <typename TLabel>
  void quicksort(double *data, TLabel *labels, int left, int right) {
  
  double pivot, tsd, fl=1;
  int pivot_indx, ll=left, rr=right;
  
  if (left < right) {
    pivot = data[left];
    while (fl==1) {
      while (data[ll] < pivot) ll++;
      while (data[rr] > pivot) rr--;
      if (ll < rr) {
        tsd = data[ll];
        data[ll] = data[rr];
        data[rr] = tsd;
        
        tsd = labels[ll];
        labels[ll] = labels[rr];
        labels[rr] = tsd;
        
        rr--;
        
      }
      else {
        pivot_indx = rr;
        fl=0;
      }
    }
    
    quicksort(data, labels, left, pivot_indx);
    quicksort(data, labels, pivot_indx+1, right);
  }
}

void GBCC(int M, int N, double* Labels, double* Data, int minleaf, int num_labels, int* bcvar, double* bcval){
  
  double *saved_logs, *sorted_data;
  double bh, ch;
  int i, j, cl, nl, mj;
  int *diff_labels_l, *diff_labels_r, *diff_labels, *sorted_labels;
  
  diff_labels_l = new int[num_labels];
  diff_labels_r = new int[num_labels];
  diff_labels   = new int[num_labels];
  
  sorted_labels = new int[M];
  saved_logs  =  new double[M];
  sorted_data =  new double[M];
  
  for(nl=0;nl<num_labels;nl++){
    diff_labels[nl]=0;
  }
  
  for(j = 0;j<M;j++) {
    saved_logs[j] = log2(j+1);
    cl = Labels[j];
    diff_labels[cl-1]++;
  }
  
  bh = 0;
  for(nl=0;nl<num_labels;nl++){
    if (diff_labels[nl]>0) bh-=diff_labels[nl]*(saved_logs[diff_labels[nl]-1]-saved_logs[M-1]);
  }
  
  for(i = 0;i<N;i++){
    
    for(nl=0;nl<num_labels;nl++){
      diff_labels_l[nl] = 0;
      diff_labels_r[nl] = diff_labels[nl];
    }
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
    }
    
    quicksort(sorted_data, sorted_labels, 0, M-1);
    
    for(mj = 0 ; mj<minleaf-1;mj++){
      cl=sorted_labels[mj];
      diff_labels_l[--cl]++;
      diff_labels_r[cl]--;
    }
    
    for(j = minleaf-1;j<M-minleaf;j++){
      
      cl=sorted_labels[j];
      diff_labels_l[--cl]++;
      diff_labels_r[cl]--;
      ch = 0;
      
      for(nl=0;nl<num_labels;nl++) {
        if(diff_labels_l[nl]>0) ch-=(diff_labels_l[nl])*(saved_logs[diff_labels_l[nl]-1]-saved_logs[j]);
        if(diff_labels_r[nl]>0) ch-=(diff_labels_r[nl])*(saved_logs[diff_labels_r[nl]-1]-saved_logs[M-j-2]);
      }
      
      if (ch<bh){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          
          bh=ch;
          bcvar[0] = i+1;
          bcval[0] = 0.5*(sorted_data[j+1]+sorted_data[j]);
        }
      }
    }
  }
  delete[] diff_labels_l;
  delete[] diff_labels_r;
  delete[] diff_labels;
  
  delete[] sorted_labels;
  delete[] sorted_data;
  
  delete[] saved_logs;
}


