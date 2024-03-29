#include <cmath>
#include "quickie.h"

void infogain_split(double lambda, int M, int N, double* Labels, double* Data, int minleaf, int numLabels, 
          int* bcvar, double* bcval, double* bestval){
    
    double *saved_logs, *sorted_data;
    double ah, bh, ch, bh0, gr, gl, t, tl, tr;
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
        //if (diff_labels[nl]>0) bh-=(diff_labels[nl]/M)*(saved_logs[diff_labels[nl]-1]-saved_logs[M-1]);
        if (diff_labels[nl]>0) bh-=(diff_labels[nl])*(saved_logs[diff_labels[nl]-1]-saved_logs[M-1]);
    }
    //bh = bh/M;
    //bh = bh*(M/(M-1))^2;
    if (lambda==M){
      t=log(M);
    }else{
      t=lambda;
    }
    bh = bh*pow(M,3)/pow(M-t,2);
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
            
            gl = 0;
            gr = 0;
            for(nl=0;nl<numLabels;nl++) {
                //if(diff_labels_l[nl]>0) gl-=(diff_labels_l[nl]/(j+1))*(saved_logs[diff_labels_l[nl]-1]-saved_logs[j]);
                //if(diff_labels_r[nl]>0) gr-=(diff_labels_r[nl]/(M-j-1))*(saved_logs[diff_labels_r[nl]-1]-saved_logs[M-j-2]);
                if(diff_labels_l[nl]>0) gl-=(diff_labels_l[nl])*(saved_logs[diff_labels_l[nl]-1]-saved_logs[j]);
                if(diff_labels_r[nl]>0) gr-=(diff_labels_r[nl])*(saved_logs[diff_labels_r[nl]-1]-saved_logs[M-j-2]);
            }
            //ch = ch/M;
            if (lambda==M){
            tl=log(j+1);
            tr=log(M-j-1);
            }else{
            tl=lambda;
            tr=lambda;
            }
            //ch = gl+gr;
            ch = gl*pow(j+1,3)/pow(j+1-tl,2) + gr*pow(M-j-1,3)/pow(M-j-1-tr,2);
            //ch = ((j+1)*gl/M) + ((M-j-1)*gr/M);
            
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
        
        bestval[i] = bh0-ah;
    }
    delete[] diff_labels_l;
    delete[] diff_labels_r;
    delete[] diff_labels;
    
    delete[] sorted_labels;
    delete[] sorted_data;
    
    delete[] saved_logs;
}


void infogain_split(double lambda, int M, int N, double* Labels, double* Data, double* W, int minleaf, 
          int numLabels, int* bcvar, double* bcval, double* bestval){
    
    double *sorted_data, *sorted_w;
    double ah, bh, ch, bh0, sum_W, sum_l, gr, gl, t, tl, tr;
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
        if (diff_labels[nl]>0) bh-=(diff_labels[nl]/sum_W)*(log2(diff_labels[nl])-log2(sum_W));
    }
    //bh = bh/sum_W;
    if (lambda==M){
      t=log(M);
    }else{
      t=lambda;
    }
    bh = bh*pow(M,3)/pow(M-t,2);
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
            
            gl = 0;
            gr = 0;
            for(nl=0;nl<numLabels;nl++) {
                if(diff_labels_l[nl]>0) gl-=(diff_labels_l[nl]/sum_W)*(log2(diff_labels_l[nl])-log2(sum_l));//(sum_l/sum_W)*
                if(diff_labels_r[nl]>0) gr-=(diff_labels_r[nl]/sum_W)*(log2(diff_labels_r[nl])-log2(sum_W-sum_l));//(sum_W-sum_l)/sum_W/sum_W*
            }
            //ch = ch/sum_W/sum_W;
            if (lambda==M){
            tl=log(j+1);
            tr=log(M-j-1);
            }else{
            tl=lambda;
            tr=lambda;
            }
            ch = gl*pow(j+1,3)/pow(j+1-tl,2) + gr*pow(M-j-1,3)/pow(M-j-1-tr,2);
            
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
        
        bestval[i] = bh0-ah;
    }
    delete[] diff_labels_l;
    delete[] diff_labels_r;
    delete[] diff_labels;
    
    delete[] sorted_labels;
    delete[] sorted_data;
    delete[] sorted_w;
}
