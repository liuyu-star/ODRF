#include <RcppArmadillo.h>  
// [[Rcpp::depends(RcppArmadillo)]] 
using namespace arma; 
using namespace Rcpp;

List VecSort(NumericVector IDdata,IntegerVector Auxdata) {
  NumericVector sortX=clone(IDdata);
  IntegerVector sortY=clone(Auxdata);  
  int n=sortX.size();  
  for(int i=0;i<(n-1);i++){
    for(int j=(i+1);j<n;j++){
      if(sortX(j)<sortX(i)){
        double temp=sortX(i);
        sortX(i)=sortX(j);
        sortX(j)=temp;
        int tempA=sortY(i);
        sortY(i)=sortY(j);
        sortY(j)=tempA;
      }
    }
  }   
  return List::create(_["sortID"]=sortX,_["sortAux"]=sortY);
}

NumericMatrix NormalizeProj(NumericMatrix proj){
  int p=proj.nrow(),q=proj.ncol();
  NumericMatrix normproj(p,q);
  double ss=0;
  for(int i=0;i<p;i++){
    ss+=proj(i,0)*proj(i,0);
  }
  for(int i=0;i<p;i++){
    normproj(i,0)=proj(i,0)/sqrt(ss);
  }    
  if(q>1){
    for(int j=1;j<q;j++){
      double temp1=0,temp2=0;
      for(int i=0;i<p;i++){
        temp1+=proj(i,j)*normproj(i,j-1);
        temp2+=normproj(i,j-1)*normproj(i,j-1);
      }
      for(int i=0;i<p;i++){
        normproj(i,j)=proj(i,j)-normproj(i,j-1)*temp1/temp2; 
      }
    }
  }
  return normproj; 
}


double LDAindex(IntegerVector y, NumericMatrix X, 
                NumericMatrix proj=NumericMatrix(0),bool weight=true){
  double index;
  int n=X.nrow(),p=X.ncol(),q=proj.ncol(),p1=proj.nrow();
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(y);
  int g=gn.size();  
  if(p1!=p)
    q=p;
  NumericVector allmean(q);
  NumericMatrix W(q,q),WB(q,q),gsum(q,g);        
  NumericMatrix projdata(n,q);
  if(p1!=p||p1==1){
    projdata=X; 
  } else{
    for(int i=0;i<n;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          projdata(i,j)+=X(i,k)*proj(k,j);
        }
      }
    }
  }
  for(int i=0;i<n;i++){
    for(int k=0;k<q;k++){
      allmean(k)+=projdata(i,k)/n;
      gsum(k,(y(i)-1))+=projdata(i,k);
    }
  }
  for(int i=0;i<n;i++){
    int l=y[i]-1;
    double gn1;
    if(weight){
      gn1=gn(l);
    } else{
      gn1=n/g; 
    }
    for(int j1=0;j1<q;j1++){
      for(int j2=0;j2<=j1;j2++){
        W(j1,j2)+=((projdata(i,j1)-gsum(j1,l)/gn(l))*
          (projdata(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
        W(j2,j1)=W(j1,j2);
        double temp=((projdata(i,j1)-gsum(j1,l)/gn(l))*
                     (projdata(i,j2)-gsum(j2,l)/gn(l))+
                     (gsum(j1,l)/gn(l)-allmean(j1))*
                     (gsum(j2,l)/gn(l)-allmean(j2)))/gn(l)*gn1;
        WB(j1,j2)+=temp;
        WB(j2,j1)=WB(j1,j2);
      }
    }
  }
  Function det=base["det"];
  index=1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
  return index;
}


double PDAindex(IntegerVector y, NumericMatrix X,
                NumericMatrix proj=NumericMatrix(0),bool weight=true,
                double lambda=0.1){
  double index;
  int n=X.nrow(),p=X.ncol(),q=proj.ncol(),p1=proj.nrow();
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(y);
  int g=gn.size();  
  NumericMatrix W(p,p),WB(p,p),gsum(p,g);
  NumericVector allmean(p);
  if(p1!=p)
    q=p;
  for(int i=0;i<n;i++){
    for(int k=0;k<p;k++){
      allmean(k)+=X(i,k)/n;
      gsum(k,(y(i)-1))+=X(i,k);
    }
  }
  for(int i=0;i<n;i++){
    int l=y[i]-1;
    double gn1;
    if(weight){
      gn1=gn(l);
    } else{
      gn1=n/g; 
    }
    for(int j1=0;j1<p;j1++){
      for(int j2=0; j2<=j1; j2++){
        double temp1,temp2;
        if(j1!=j2){
          temp1=(1-lambda)*((X(i,j1)-gsum(j1,l)/gn(l))*
            (X(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
          
          temp2=(1-lambda)*((X(i,j1)-gsum(j1,l)/gn(l))*
            (X(i,j2)-gsum(j2,l)/gn(l)))+
            (gsum(j1,l)/gn(l)-allmean(j1))*
            (gsum(j2,l)/gn(l)-allmean(j2))/gn(l)*gn1;
        } else{
          temp1=((X(i,j1)-gsum(j1,l)/gn(l))*
            (X(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
          
          temp2=((X(i,j1)-gsum(j1,l)/gn(l))*
            (X(i,j2)-gsum(j2,l)/gn(l))+
            (gsum(j1,l)/gn(l)-allmean(j1))*
            (gsum(j2,l)/gn(l)-allmean(j2)))/gn(l)*gn1;              
        }  
        W(j1,j2)+=temp1;  
        WB(j1,j2)+=temp2;
        W(j2,j1)=W(j1,j2);            
        WB(j2,j1)=WB(j1,j2);
      }
    }      
  }
  NumericMatrix Wt(q,p),WBt(q,p);    
  NumericMatrix Wtt(q,q),WBtt(q,q); 
  if(p1!=p||p1==1){
    Wtt=W;
    WBtt=WB;
  } else{
    for(int i=0;i<p;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          Wt(j,i)+=W(k,i)*proj(k,j);
          WBt(j,i)+=WB(k,i)*proj(k,j);               
        }
      }
    }    
    for(int i=0;i<q;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          Wtt(i,j)+=Wt(i,k)*proj(k,j);
          WBtt(i,j)+=WBt(i,k)*proj(k,j);               
        }
      }
    }      
  }   
  Function det=base["det"];
  index=1.0-as<double>(det(wrap(Wtt)))/as<double>(det(wrap(WBtt)));   
  return index;
}


double Lrindex(IntegerVector y, NumericMatrix X,
               NumericMatrix proj=NumericMatrix(0),bool weight=true,int r=1){
  double index,B=0,W=0;
  int n=X.nrow(),p=X.ncol(),q=proj.ncol(),p1=proj.nrow();
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(y);
  int g=gn.size();   
  if(p1!=p)
    q=p;
  NumericVector allmean(q);
  NumericMatrix gsum(q,g);
  NumericMatrix projdata(n,q);
  if(p1!=p||p1==1){
    projdata=X; 
  } else{
    for(int i=0;i<n;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          projdata(i,j)+=X(i,k)*proj(k,j);
        }
      }
    }
  }
  for(int i=0;i<n;i++){
    for(int k=0;k<q;k++){
      allmean(k)=allmean(k)+projdata(i,k)/n;
      gsum(k,(y(i)-1))+=projdata(i,k);
    }
  }
  for(int i=0;i<n;i++){
    int l=y(i)-1;
    double gn1;
    if(weight){
      gn1=gn(l);
    } else{
      gn1=n/g; 
    }         
    for(int j=0;j<q;j++){
      W=W+pow(pow((projdata(i,j)-gsum(j,l)/gn(l)),2*r),0.5)/gn(l)*gn1;
      B=B+pow(pow((gsum(j,l)/gn(l)-allmean(j)),2*r),0.5)/gn(l)*gn1;
    }
  } 
  W=pow(W,1.0/r); 
  double WB=pow(W+B,1.0/r);
  index=1-W/WB;
  return index;
}


double GINIindex1D(IntegerVector y,NumericMatrix X,
                   NumericVector proj=NumericVector(0)){ 
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(y);
  int g=gn.size();
  int n=X.nrow(),p=X.ncol(),p1=proj.size();
  double n1,n2;
  double index=0,tempindex;
  NumericVector projdata(n);
  if(p1!=p||p1 ==1){
    projdata=X(_,0);
  } else{
    for(int i=0;i<n;i++){
      for(int k=0;k<p;k++){
        projdata(i)+=X(i,k)*proj(k);
      }
    }
  } 
  List VecSortdata=VecSort(projdata,y);
  NumericVector sortdata=as<NumericVector>(VecSortdata["sortID"]);
  IntegerVector sortclass=as<IntegerVector>(VecSortdata["sortAux"]);
  IntegerVector part1,part2,temptable1,temptable2;
  for(int i=1;i<(n-1);i++){  
    part1=sortclass[sortdata<=sortdata[i]];
    part2=sortclass[sortdata>sortdata[i]];
    n1=part1.size();
    n2=part2.size();
    temptable1=table(part1); int g1=temptable1.size();  tempindex=0;
    temptable2=table(part2); int g2=temptable2.size(); 
    for(int j=0;j<g1;j++){
      tempindex+=((n1)/n)*(temptable1(j)/n1)*(1.0-temptable1(j)/n1);
    }
    for(int j=0;j<g2;j++){         
      tempindex+=((n2)/n)*(temptable2(j)/n2)*(1.0-temptable2(j)/n2);        
    } 
    tempindex=g-1.0-g*tempindex;
    if(tempindex>index) 
      index=tempindex;             
  }  
  return index;
}


double ENTROPYindex1D(IntegerVector y,NumericMatrix X,
                      NumericVector proj=NumericVector(0)){
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(y);
  int g=gn.size();
  int n=X.nrow(),p=X.ncol(),p1=proj.size();
  double n1,n2;
  double index=0,tempindex;
  NumericVector projdata(n);
  if(p1!=p||p1==1){
    projdata=X(_,0);
  } else{
    for(int i=0;i<n;i++){
      for(int k=0;k<p;k++){
        projdata(i)+=X(i,k)*proj(k);
      }
    }
  } 
  List VecSortdata=VecSort(projdata,y);
  NumericVector sortdata=as<NumericVector>(VecSortdata["sortID"]);
  IntegerVector sortclass=as<IntegerVector>(VecSortdata["sortAux"]);
  IntegerVector part1,part2,temptable1,temptable2;
  for(int i=1;i<(n-1);i++){  
    part1=sortclass[sortdata<=sortdata[i]];
    part2=sortclass[sortdata>sortdata[i]];
    n1=part1.size();
    n2=part2.size();
    temptable1=table(part1); int g1=temptable1.size();  tempindex=0;
    temptable2=table(part2); int g2=temptable2.size(); 
    for(int j=0; j<g1; j++){
      if(temptable1(j)!=0)
      {  tempindex+=((n1)/n)*(temptable1(j)/n1)*log(temptable1(j)/n1);
      }
    }
    for(int j=0; j<g2; j++){         
      if(temptable2(j)!=0)
      {  tempindex+=((n2)/n)*(temptable2(j)/n2)*log(temptable2(j)/n2);        
      }
    } 
    double maxI=log(2.0)-log(g/1.0);
    if((g/2)*2!=g){
      maxI= -0.5*log((g*g-1.0)/4.0)+1.0/(2.0*g)*log((g-1.0)/(g+1.0));
    }      
    tempindex=(1+tempindex/log(g/1.0))/(1+maxI/log(g/1.0));
    if(tempindex>index) index=tempindex;            
  }   
  return index;
}

// [[Rcpp::export(name="ppOptCpp")]]
List ppOpt(IntegerVector y,NumericMatrix X,int q=1, 
           std::string PPmethod="LDA",bool weight=true,int r=1,
           double lambda=0.1,double energy=0,double cooling=0.999, 
           double TOL=0.0001,int maxiter=50000){
  int n=X.nrow(),p=X.ncol();
  Environment base("package:base");
  Function table=base["table"];
  GetRNGstate();
  NumericMatrix projbest(p,q);
  double indexbest=0,newindex=0;
  if((PPmethod=="GINI"||PPmethod=="ENTROPY") && q>1){
    return Rcpp::List::create(Rcpp::Named("indexbest") = 0,
                              Rcpp::Named("projbest") = 0);
  } else {  
    if(PPmethod=="LDA"){
      indexbest=LDAindex(y,X,NumericMatrix(0),weight);
    } else if(PPmethod=="Lr"){
      indexbest=Lrindex(y,X,NumericMatrix(0),weight,r);
    } else if(PPmethod=="PDA"){
      indexbest=PDAindex(y,X,NumericMatrix(0),weight,lambda);  
    } else if(PPmethod=="GINI"){
      double tempindex=0;
      for(int i=0; i<p; i++){
        NumericVector tempproj(p);
        tempproj(i)=1;
        tempindex=GINIindex1D(y,X,tempproj);  
        if(indexbest<tempindex)
          indexbest=tempindex;
      }        
    } else if(PPmethod=="ENTROPY"){
      double tempindex=0;
      for(int i=0; i<p; i++){
        NumericVector tempproj(p);
        tempproj(i)=1;
        tempindex=ENTROPYindex1D(y,X,tempproj);  
        if(indexbest<tempindex)
          indexbest=tempindex;
      } 
    }    
    if(energy==0)
      energy=1-indexbest;
    
    for(int k=0;k<q;k++){
      projbest(_,k)=rnorm(p);
    }
    projbest=NormalizeProj(projbest);  
    NumericMatrix projdata(n,q);    
    for(int i=0;i<n;i++){
      for(int k=0;k<q;k++){
        projdata(i,k)=0;
        for(int j=0; j<p; j++){
          projdata(i,k)+=X(i,j)*projbest(j,k);
        }
      }
    }  
    if(PPmethod=="LDA"){
      indexbest=LDAindex(y,X,projbest,weight);
    } else if(PPmethod=="Lr"){
      indexbest=Lrindex(y,X,projbest,weight,r);
    } else if(PPmethod=="PDA"){
      indexbest=PDAindex(y,X,projbest,weight,lambda);
    } else if(PPmethod=="GINI"){
      indexbest=GINIindex1D(y,X,projbest);        
    } else if(PPmethod=="ENTROPY"){
      indexbest=ENTROPYindex1D(y,X,projbest);        
    }   
    double temp=1;
    int kk=0;  
    double diff=100;
    while(fabs(diff)>TOL&&kk<maxiter){
      double tempp=energy/log(kk+2.0);
      if(kk>1000) {
        temp=temp*cooling;
      } else {
        temp=temp*cooling*cooling;
      }   
      
      NumericMatrix projnew(p,q);
      for(int k=0;k<q;k++){
        projnew(_,k)=temp*rnorm(p)+projbest(_,k);
      }
      projnew=NormalizeProj(projnew);  
      
      for(int i=0;i<n;i++){
        for(int k=0;k<q;k++){
          projdata(i,k)=0;
          for(int j=0; j<p; j++){
            projdata(i,k)+=X(i,j)*projnew(j,k);
          }
        }
      }      
      if(PPmethod=="LDA"){
        newindex=LDAindex(y,X,projnew,weight);
      } else if(PPmethod=="Lr"){
        newindex=Lrindex(y,X,projnew,weight,r);
      } else if(PPmethod=="PDA"){
        newindex=PDAindex(y,X,projnew,weight,lambda);
      } else if(PPmethod=="GINI"){
        NumericVector tempproj= projnew(_,0);
        newindex=GINIindex1D(y,X,tempproj);        
      } else if(PPmethod=="ENTROPY"){
        NumericVector tempproj= projnew(_,0);
        newindex=ENTROPYindex1D(y,X,tempproj);        
      }      
      NumericVector prob=runif(1);
      double difft=newindex-indexbest;
      double e=exp(difft/tempp);
      if(e>1){
        for(int i=0;i<p;i++){
          for(int j=0;j<q;j++){
            projbest(i,j)=projnew(i,j);
          }
        }
        indexbest=newindex;    
        diff=difft;
      } else if(prob[0]<e&&difft>energy){
        for(int i=0;i<p;i++){
          for(int j=0; j<q; j++){
            projbest(i,j)=projnew(i,j);
          }
        }
        indexbest=newindex;    
        diff=difft; 
      }
      kk++;
    }
    PutRNGstate();
    return Rcpp::List::create(Rcpp::Named("indexbest")=indexbest,
                              Rcpp::Named("projbest")=projbest);
  }                           
}

