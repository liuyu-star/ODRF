#include <Rcpp.h>
#include <stdio.h>
#include <cmath>

using namespace Rcpp;

//[[Rcpp::export("predict_ODT")]]
Rcpp::List predict_ODT(NumericMatrix Data,NumericMatrix nodeRotaMat, NumericVector nodeCutValue,
                              NumericVector childNode, CharacterVector nodeLabel) {
  
  int i,t,ppVar,current_node;
  int N=Data.nrow();
  int numNode=nodeCutValue.size();
  //int *LenRotaMat;
  //LenRotaMat = new int[numNode+1];
  double rotaX,ppDr,cutVal;
  
  IntegerVector LenRotaMat(numNode+1);
  CharacterVector prediction(N);
  NumericVector node(N);
  
  LenRotaMat[0]=0;
  for (current_node = 1; current_node <= numNode; current_node++)
  {
    LenRotaMat[current_node]=LenRotaMat[current_node-1]+sum(nodeRotaMat.column(1)==current_node);
  }
  
  for(i = 0;i<N;i++){
    current_node = 0;
    
    while (childNode[current_node]!=0){
      
      rotaX=0; 
      for (t = LenRotaMat[current_node]; t < LenRotaMat[current_node+1]; t++)
      {
        ppVar=nodeRotaMat(t,0);
        ppDr=nodeRotaMat(t,2);
        rotaX += Data(i,ppVar-1) * ppDr;
      }
      
      cutVal=nodeCutValue[current_node];
      if (rotaX < cutVal){
        current_node = childNode[current_node]-1;
      }else{ 
        current_node = childNode[current_node];
      }
      
    }
    
    node[i] = current_node+1;
    prediction[i] = nodeLabel[current_node];
  }
  
  List Pred = List::create(
    _["prediction"]= prediction, 
    _["node"]= node);
  
  return(Pred);
}

