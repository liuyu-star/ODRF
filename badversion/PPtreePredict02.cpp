#include <Rcpp.h>
#include <stdio.h>
#include <cmath>

using namespace Rcpp;

//[[Rcpp::export("PPtreePredictCpp")]]
CharacterVector PPtreePredict(NumericMatrix Data,NumericMatrix nodePPMat, NumericVector cut_var, 
                 NumericVector cut_val, NumericVector nodechilds, CharacterVector nodelabel) {

    int i,current_node,cvar;
    int M=Data.nrow();
    int tn=nodePPMat.nrow();
    int N=cut_var.size();
    int k0,cut_var_n0,ppVar,j=0,t0[N+1];
    double xpi,ppDr,cval;
    
    CharacterVector Prediction(M);
    

    t0[j]=0;
    for (int t = 0; t < tn; t++)
     {
        if (nodePPMat(t,0) ==0)
        {   
            j++;
            t0[j] =t+1;
        }
    }
      
    for(i = 0;i<M;i++){
        current_node = 0;
        k0=0;
        cut_var_n0=0;

        while (nodechilds[current_node]!=0){
            cvar = cut_var[current_node];
            
            for (int k = k0; k <= current_node; k++)
            {
                if (cut_var(k) ==0)
                {   
                    cut_var_n0 ++;
                }
            }
            k0=current_node+1;

            xpi=0;
            for (int t = t0[current_node-cut_var_n0]; t < (t0[current_node-cut_var_n0+1]-1); t++)
             {
                if (nodePPMat(t,1)==cvar)
                {
                    ppVar=nodePPMat(t,0);
                    ppDr=nodePPMat(t,2);
                    xpi += Data(i,ppVar-1) * ppDr;
                }
            }
            
            cval=cut_val[current_node];
            if (xpi < cval){
                current_node = nodechilds[current_node]-1;
            }else{ 
                current_node = nodechilds[current_node];
            }
        }

        Prediction[i] = nodelabel[current_node];
    }
    
    return(Prediction);
}

