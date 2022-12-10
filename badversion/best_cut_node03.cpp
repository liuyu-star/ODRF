#include <Rcpp.h>
//#include "node_cuts.h"
//#include "quickie.h"
//#include <Rmath.h>
#include <cmath>
using namespace Rcpp;
//using namespace std; 

template <typename TLabel>
void quicksort(double *data, TLabel *labels, int left, int right);
void gini_split(int , int , double* , double* , int , int , int* , double*, double* );
void infogain_split(int , int , double* , double* , int , int , int* , double*, double* );
void mse_split(int M, int N, double* Labels, double* Data, int minleaf, int* bcvar, double* bcval, double* bestval);
template <typename TLabel>
void quicksort(double *data, TLabel *labels, double *W, int left, int right);
void gini_split(int M, int N, double* Labels, double* Data, double* W, int minleaf, int numLabels, int* bcvar, double* bcval, double* bestval);
void infogain_split(int M, int N, double* Labels, double* Data, double* W, int minleaf, 
                    int numLabels, int* bcvar, double* bcval, double* bestval);
void mse_split(int M, int N, double* Labels, double* Data, double* W, int minleaf, 
               int* bcvar, double* bcval, double* bestval);

// [[Rcpp::export]]
List best_cut_node(char method, NumericMatrix Data, NumericVector Labels, NumericVector W,int minleaf, int numLabels) {
  
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
      gini_split(M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'i':
      infogain_split(M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'r':
      mse_split(M, N, Labels1, Data1, minleaf, &bcvar, &bcval, bestval);
      break;
    }
  }else{
    double *W1 = new double[M];
    for(int k=0;k<M;k++){
      W1[k]=W[k];
    }
    
    switch (method){
    case 'g':
      gini_split(M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'i':
      infogain_split(M, N, Labels1, Data1, minleaf, numLabels, &bcvar, &bcval, bestval);
      break;
    case 'r':
      mse_split(M, N, Labels1, Data1, minleaf, &bcvar, &bcval, bestval);
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
    _["BestVal"]= BestVal);
  
  delete[] Data1;
  delete[] Labels1;
  delete[] bestval;
  
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

void gini_split(int M, int N, double* Labels, double* Data, int minleaf, 
                int numLabels, int* bcvar, double* bcval, double* bestval){
  
  double *sorted_data;
  double ah, bh, ch, gr, gl;
  int i, j, cl, nl, mj;
  int *diff_labels_l, *diff_labels_r, *diff_labels, *sorted_labels;
  
  diff_labels_l = new int[numLabels];
  diff_labels_r = new int[numLabels];
  diff_labels = new int[numLabels];
  sorted_labels = new int[M];
  sorted_data =  new double[M];
  
  for(nl=0;nl<numLabels;nl++){
    diff_labels[nl]=0;
  }
  
  for(j = 0;j<M;j++) {
    cl = Labels[j];
    diff_labels[cl-1]++;
  }
  
  bh=0;
  for(nl=0;nl<numLabels;nl++){
    bh+=diff_labels[nl]*diff_labels[nl];
  }
  bh = 1 - (bh/(M*M));
  
  
  for(i = 0;i<N;i++){
    ah=1e+10;
    
    for(nl=0;nl<numLabels;nl++){
      diff_labels_l[nl] = 0;
      diff_labels_r[nl] = diff_labels[nl];
    }
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
    }
    
    quicksort(sorted_data, sorted_labels, 0, M-1);
    
    for(mj = 0 ; mj<minleaf;mj++){
      //cl=sorted_labels[j];
      //diff_labels_l[--cl]++;
      //diff_labels_r[cl]--;
      cl=sorted_labels[mj]-1;
      diff_labels_l[cl]++;
      diff_labels_r[cl]--;
    }
    
    for(j = minleaf;j<M-minleaf;j++){
      //cl=sorted_labels[j];
      //diff_labels_l[--cl]++;
      //diff_labels_r[cl]--;
      cl=sorted_labels[j]-1;
      diff_labels_l[cl]++;
      diff_labels_r[cl]--;
      
      gr = 0;
      gl = 0;
      for(nl=0;nl<numLabels;nl++) {
        gl+=diff_labels_l[nl]*diff_labels_l[nl];
        gr+=diff_labels_r[nl]*diff_labels_r[nl];
      }
      gl = 1 - gl/((j+1)*(j+1));
      gr = 1 - gr/((M-j-1)*(M-j-1));
      
      ch = ((j+1)*gl/M) + ((M-j-1)*gr/M);
      
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
  
  delete[] diff_labels_l;
  delete[] diff_labels_r;
  delete[] diff_labels;
  
  delete[] sorted_labels;
  delete[] sorted_data;
}

void infogain_split(int M, int N, double* Labels, double* Data, int minleaf, int numLabels, 
                    int* bcvar, double* bcval, double* bestval){
  
  double *saved_logs, *sorted_data;
  double ah, bh, ch, bh0;
  int i, j, cl, nl, mj;
  int *diff_labels_l, *diff_labels_r, *diff_labels, *sorted_labels;
  
  diff_labels_l = new int[numLabels];
  diff_labels_r = new int[numLabels];
  diff_labels   = new int[numLabels];
  
  sorted_labels = new int[M];
  saved_logs  =  new double[M];
  sorted_data =  new double[M];
  
  for(nl=0;nl<numLabels;nl++){
    diff_labels[nl]=0;
  }
  
  for(j = 0;j<M;j++) {
    saved_logs[j] = log2(j+1);
    cl = Labels[j];
    diff_labels[cl-1]++;
  }
  
  bh = 0;
  for(nl=0;nl<numLabels;nl++){
    if (diff_labels[nl]>0) bh-=diff_labels[nl]*(saved_logs[diff_labels[nl]-1]-saved_logs[M-1]);
  }
  bh = bh/M;
  bh0=bh;
  
  for(i = 0;i<N;i++){
    ah=1e-10;
    
    for(nl=0;nl<numLabels;nl++){
      diff_labels_l[nl] = 0;
      diff_labels_r[nl] = diff_labels[nl];
    }
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
    }
    
    quicksort(sorted_data, sorted_labels, 0, M-1);
    
    for(mj = 0 ; mj<minleaf;mj++){
      //cl=sorted_labels[mj]-1;
      //diff_labels_l[--cl]++;
      //diff_labels_r[cl]--;
      cl=sorted_labels[mj]-1;
      diff_labels_l[cl]++;
      diff_labels_r[cl]--;
    }
    
    for(j = minleaf;j<M-minleaf;j++){
      //cl=sorted_labels[j];
      //diff_labels_l[--cl]++;
      //diff_labels_r[cl]--;
      cl=sorted_labels[j]-1;
      diff_labels_l[cl]++;
      diff_labels_r[cl]--;
      
      ch = 0;
      for(nl=0;nl<numLabels;nl++) {
        if(diff_labels_l[nl]>0) ch-=diff_labels_l[nl]*(saved_logs[diff_labels_l[nl]-1]-saved_logs[j]);
        if(diff_labels_r[nl]>0) ch-=diff_labels_r[nl]*(saved_logs[diff_labels_r[nl]-1]-saved_logs[M-j-2]);
      }
      ch = ch/M;
      
      if (ch<bh){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          bh=ch;
          *bcvar = i+1;
          *bcval = 0.5*(sorted_data[j+1]+sorted_data[j]);
        }
      }
      
      if (bh0-ch>ah){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          ah=bh0-ch;
        }
      }
    }
    
    bestval[i]= ah;
  }
  delete[] diff_labels_l;
  delete[] diff_labels_r;
  delete[] diff_labels;
  
  delete[] sorted_labels;
  delete[] sorted_data;
  
  delete[] saved_logs;
}

void mse_split(int M, int N, double* Labels, double* Data, int minleaf, 
               int* bcvar, double* bcval, double* bestval){
  
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

////////////////////////////////////////////////////////////////////////////////
template <typename TLabel>
void quicksort(double *data, TLabel *labels, double *W, int left, int right) {
  
  
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
        
        tsd = W[ll];
        W[ll] = W[rr];
        W[rr] = tsd;
        
        rr--;
      }
      else {
        pivot_indx = rr;
        fl=0;
      }
    }
    
    quicksort(data, labels, W, left, pivot_indx);
    quicksort(data, labels, W, pivot_indx+1, right);
  }
}

void gini_split(int M, int N, double* Labels, double* Data, double* W, int minleaf, 
                int numLabels, int* bcvar, double* bcval, double* bestval){
  
  double *sorted_data, *sorted_w;
  double ah, bh, ch, sum_W, sum_l, gr, gl;
  int i, j, cl, nl, mj;
  double *diff_labels_l, *diff_labels_r, *diff_labels;
  int *sorted_labels;
  
  sorted_labels = new int[M];
  sorted_data   = new double[M];
  sorted_w      = new double[M];
  
  diff_labels_l = new double[numLabels];
  diff_labels_r = new double[numLabels];
  diff_labels   = new double[numLabels];
  
  for(nl=0;nl<numLabels;nl++){
    diff_labels[nl]=0;
  }
  
  sum_W=0;
  for(j = 0;j<M;j++) {
    cl = Labels[j];
    diff_labels[cl-1]+=W[j];
    sum_W+=W[j];
  }
  
  bh=0;
  for(nl=0;nl<numLabels;nl++){
    bh+=diff_labels[nl]*diff_labels[nl];
  }
  bh = 1 - (bh/(sum_W*sum_W));
  
  for(i = 0;i<N;i++){
    ah=1e+10;
    
    for(nl=0;nl<numLabels;nl++){
      diff_labels_l[nl] = 0;
      diff_labels_r[nl] = diff_labels[nl];
    }
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
      sorted_w[j] = W[j];
    }
    
    sum_l=0;
    
    quicksort(sorted_data, sorted_labels, sorted_w, 0, M-1);
    
    for(mj = 0 ; mj<minleaf;mj++){
      //cl=sorted_labels[mj];
      //diff_labels_l[--cl]+=sorted_w[mj];
      //diff_labels_r[cl]-=sorted_w[mj];
      cl=sorted_labels[mj]-1;
      diff_labels_l[cl]+=sorted_w[mj];
      diff_labels_r[cl]-=sorted_w[mj];
      sum_l+=sorted_w[mj];
    }
    
    for(j = minleaf;j<M-minleaf;j++){
      //cl=sorted_labels[j];
      //diff_labels_l[--cl]+=sorted_w[j];
      //diff_labels_r[cl]-=sorted_w[j];
      cl=sorted_labels[j]-1;
      diff_labels_l[cl]+=sorted_w[j];
      diff_labels_r[cl]-=sorted_w[j];
      sum_l += sorted_w[j];
      
      gr = 0;
      gl = 0;
      
      for(nl=0;nl<numLabels;nl++) {
        gl+=diff_labels_l[nl]*diff_labels_l[nl];
        gr+=diff_labels_r[nl]*diff_labels_r[nl];
      }
      gl = 1 - gl/(sum_l*sum_l);
      gr = 1 - gr/((sum_W-sum_l)*(sum_W-sum_l));
      
      ch = ((sum_l)*gl/sum_W) + ((sum_W-sum_l)*gr/sum_W);
      
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
  
  delete[] diff_labels_l;
  delete[] diff_labels_r;
  delete[] diff_labels;
  
  delete[] sorted_labels;
  delete[] sorted_data;
  delete[] sorted_w;
}

void infogain_split(int M, int N, double* Labels, double* Data, double* W, int minleaf, 
                    int numLabels, int* bcvar, double* bcval, double* bestval){
  
  double *sorted_data, *sorted_w;
  double ah, bh, ch, bh0, sum_W, sum_l;
  int i, j, cl, nl, mj;
  double *diff_labels_l, *diff_labels_r, *diff_labels;
  int *sorted_labels;
  
  sorted_labels = new int[M];
  sorted_data   = new double[M];
  sorted_w      = new double[M];
  
  diff_labels_l = new double[numLabels];
  diff_labels_r = new double[numLabels];
  diff_labels   = new double[numLabels];
  
  for(nl=0;nl<numLabels;nl++){
    diff_labels[nl]=0;
  }
  
  sum_W=0;
  for(j = 0;j<M;j++) {
    cl = Labels[j];
    diff_labels[cl-1]+=W[j];
    sum_W+=W[j];
  }
  
  bh = 0;
  for(nl=0;nl<numLabels;nl++){
    if (diff_labels[nl]>0) bh-=diff_labels[nl]*(log2(diff_labels[nl])-log2(sum_W));
  }
  bh = bh/sum_W;
  bh0=bh;
  
  for(i = 0;i<N;i++){
    ah=1e-10;
    
    for(nl=0;nl<numLabels;nl++){
      diff_labels_l[nl] = 0;
      diff_labels_r[nl] = diff_labels[nl];
    }
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
      sorted_w[j] = W[j];
    }
    
    sum_l=0;
    
    quicksort(sorted_data, sorted_labels, sorted_w, 0, M-1);
    
    for(mj = 0 ; mj<minleaf;mj++){
      //cl=sorted_labels[mj];
      //diff_labels_l[--cl]+=sorted_w[mj];
      //diff_labels_r[cl]-=sorted_w[mj];
      cl=sorted_labels[mj]-1;
      diff_labels_l[cl]+=sorted_w[mj];
      diff_labels_r[cl]-=sorted_w[mj];
      sum_l+=sorted_w[mj];
    }
    
    for(j = minleaf;j<M-minleaf;j++){
      
      //cl=sorted_labels[j];
      //diff_labels_l[--cl]+=sorted_w[j];
      //diff_labels_r[cl]-=sorted_w[j];
      cl=sorted_labels[j]-1;
      diff_labels_l[cl]+=sorted_w[j];
      diff_labels_r[cl]-=sorted_w[j];
      sum_l += sorted_w[j];
      
      ch = 0;
      for(nl=0;nl<numLabels;nl++) {
        if(diff_labels_l[nl]>0) ch-=(sum_l)*(diff_labels_l[nl])*(log2(diff_labels_l[nl])-log2(sum_l));
        if(diff_labels_r[nl]>0) ch-=(sum_W-sum_l)*(diff_labels_r[nl])*(log2(diff_labels_r[nl])-log2(sum_W-sum_l));
      }
      ch = ch/sum_W/sum_W;
      
      if (ch<bh){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          bh=ch;
          *bcvar = i+1;
          *bcval = 0.5*(sorted_data[j+1]+sorted_data[j]);
        }
      }
      
      if (bh0-ch>ah){
        if (fabs(sorted_data[j+1]-sorted_data[j])>1e-15){
          ah=bh0-ch;
        }
      }
    }
    
    bestval[i]= ah;
  }
  delete[] diff_labels_l;
  delete[] diff_labels_r;
  delete[] diff_labels;
  
  delete[] sorted_labels;
  delete[] sorted_data;
  delete[] sorted_w;
}

void mse_split(int M, int N, double* Labels, double* Data, double* W, int minleaf, 
               int* bcvar, double* bcval, double* bestval){
  
  double *sorted_data, *sorted_labels, *sorted_w;
  double ah, bh, ch, sum_W, sum_L, sum_wall, sum_wr, sum_wl, sum_all, sum_all2, suml, suml2, sumr, sumr2, sr, sl, cl;
  int i, j, mj;
  
  sorted_labels = new double[M];
  sorted_data   = new double[M];
  sorted_w      = new double[M];
  
  sum_all  = 0;
  sum_wall = 0;
  sum_all2 = 0;
  sum_W = 0;
  
  for(j=0;j<M;j++){
    sum_all+=Labels[j];
    sum_wall += W[j]*Labels[j];
    sum_all2 += W[j]*Labels[j]*Labels[j];
    sum_W += W[j];
  }
  
  //bh = sum_all2 + sum_W*(sum_all/M)*(sum_all/M) - 2*(sum_all/M)*sum_wall;
  bh = (sum_all2 + sum_W*(sum_all/M)*(sum_all/M) - 2*(sum_all/M)*sum_wall)/sum_W;
  
  for(i = 0;i<N;i++){
    ah=1e+10;
    
    suml = 0;
    suml2 = 0;
    sum_wl = 0;
    
    sumr = sum_all;
    sumr2 = sum_all2;
    sum_wr = sum_wall;
    
    sum_L = 0;
    
    for(j = 0;j<M;j++){
      sorted_data[j] = Data[i*M+j];
      sorted_labels[j] = Labels[j];
      sorted_w[j] = W[j];
    }
    
    quicksort(sorted_data, sorted_labels, sorted_w, 0, M-1);
    
    for(mj = 0 ; mj<minleaf;mj++){
      cl=sorted_labels[mj];
      suml += cl;
      sumr -= cl;
      
      sum_wl +=sorted_w[mj]*cl;
      sum_wr -=sorted_w[mj]*cl;
      suml2 += sorted_w[mj]*cl*cl;
      sumr2 -= sorted_w[mj]*cl*cl;
      sum_L +=sorted_w[mj];
    }
    for(j = minleaf;j<M-minleaf;j++){
      
      cl=sorted_labels[j];
      
      suml += cl;
      sumr -= cl;
      
      sum_wl +=sorted_w[j]*cl;
      sum_wr -=sorted_w[j]*cl;
      suml2 += sorted_w[j]*cl*cl;
      sumr2 -= sorted_w[j]*cl*cl;
      sum_L +=sorted_w[j];
      
      sl = suml2 + sum_L*(suml/(j+1))*(suml/(j+1)) - 2*(suml/(j+1))*sum_wl;
      sr = sumr2 + (sum_W-sum_L)*(sumr/(M-j-1))*(sumr/(M-j-1)) - 2*(sumr/(M-j-1))*sum_wr;
      
      ch = (sl + sr)/sum_W;
      //ch = sl + sr + log(M)/M;
      
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
  delete[] sorted_w;
}