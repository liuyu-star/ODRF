#include "quickie.h"
#include <cmath>
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
