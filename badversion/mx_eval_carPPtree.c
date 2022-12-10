#include <stdio.h>
#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *Data,*tree_output,*nodePPMat, *cut_var,*cut_val,*nodechilds,*nodelabel;
    int i,current_node,cvar;
    
    int M=mxGetM(prhs[0]);
    int tn=mxGetM(prhs[1]);
    int N=mxGetM(prhs[4]);
    int k0,cut_var_n0,ppVar,j=0,t0[N+1];
    double xpi,ppDr,cval;

    Data= mxGetPr(prhs[0]);
    nodePPMat= mxGetPr(prhs[1]);
    cut_var = mxGetPr(prhs[2]);
    cut_val = mxGetPr(prhs[3]);
    nodechilds = mxGetPr(prhs[4]);
    nodelabel = mxGetPr(prhs[5]);

    t0[j]=0;
    for (int t = 0; t < tn; t++)
     {
        if (*(nodePPMat + t) ==0)
        {   
            j++;
            t0[j] =t+1;
        }
    }

    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
    tree_output = mxGetPr(plhs[0]);

    for(i = 0;i<M;i++){
        current_node = 0;
        k0=0;
        cut_var_n0=0;

        while (nodechilds[current_node]!=0){
            cvar = cut_var[current_node];
            
            for (int k = k0; k <= current_node; k++)
            {
                if (*(cut_var + k) ==0)
                {   
                    cut_var_n0 ++;
                }
            }
            k0=current_node+1;

            xpi=0;
            for (int t = t0[current_node-cut_var_n0]; t < (t0[current_node-cut_var_n0+1]-1); t++)
             {
                if (*(nodePPMat + t + tn)==cvar)
                {
                    ppVar=*(nodePPMat+ t);
                    ppDr=*(nodePPMat+ t+2*tn);
                    xpi += *(Data + i + (ppVar-1)*M) * ppDr;
                }
            }
            
            cval=cut_val[current_node];
            if (xpi < cval){
                current_node = nodechilds[current_node]-1;
            }else{ 
                current_node = nodechilds[current_node];
            }
        }

        tree_output[i] = nodelabel[current_node];
    }
}
