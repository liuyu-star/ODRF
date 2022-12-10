// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace std;
using namespace Eigen;
using namespace Rcpp;


void sort_vec(const VectorXd& vec, VectorXd& sorted_vec, VectorXi& ind) {
  ind = VectorXi::LinSpaced(vec.size(), 0, vec.size() - 1);//[0 1 2 3 ... N-1]
  auto rule = [vec](int i, int j)->bool {
    return vec(i) < vec(j);
  };
  std::sort(ind.data(), ind.data() + ind.size(), rule);
  sorted_vec.resize(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    sorted_vec(i) = vec(ind(i));
  }
}


double DCov2d(const VectorXd z,const VectorXd x,const VectorXd a_x,const VectorXi index)
{
  int i,j, n=x.size(); 
  MatrixXd it=MatrixXd::Zero(n+1,4);
  VectorXd y(n), xy(n); 
  
  //y=z(index);
  for (i = 0; i < n; i++){y(i)=z(index(i));}
  xy=x.cwiseProduct(y);
  //VectorXi seq1n= VectorXi::LinSpaced(n, 1, n); 
  //std::vector<int> seq1n= VectorXi::LinSpaced(n, 1, n); 
  
  //(Eigen::placeholders::all,ind) 
  it.block(1,0,n,1)=VectorXd::Ones(n);
  it.block(1,1,n,1)=y;
  it.block(1,2,n,1)=x;
  it.block(1,3,n,1)=xy;
  
  // id1 , id2 , id3 , id4 are used to store sum of weights. 
  // On output , for 1 ?? j ?? n.
  // b_y is the vector of row sums of distance matrix of y
  MatrixXd id=MatrixXd::Zero(n,4);
  VectorXi idy(n);
  VectorXd sy(n),si(n);
  sort_vec(y, sy, idy);  
  
  VectorXi keep_cols = VectorXi::LinSpaced(it.cols(), 0, it.cols()-1);
  si(0)=sy(0);
  VectorXi seqi,syi;
  for (i = 1; i < n; i++){
    seqi= VectorXi::LinSpaced(i, 1, i ); 
    si(i)=sy.head(i+1).sum();
    
    syi=(y.head(i).array()<y(i)).cast<int>();
    syi=syi.cwiseProduct(seqi);
    for (j = 0; j < syi.size(); j++){
      id.row(i)+=it.block(syi(j),0,1,n);//.colwise().sum();
    }
  }
  
  VectorXd b_y(n),by(n);
  by=2*VectorXd::LinSpaced(n, 1, n).array()-n;
  by =sy.cwiseProduct(by);
  //b_y(idy) =by.array()+(si(n-1)-(2*si).array()).array();
  by=by.array()+(si(n-1)-(2*si).array()).array();
  for (i = 0; i < n; i++){b_y(idy(i)) =by(i);}
  
  // the Frobenius inner product of the distance matrices 
  VectorXd covtermx=x.array()-x.mean();
  VectorXd covtermy=y.array()-y.mean();
  
  double D = 4*(xy.dot(id.col(0))-x.dot(id.col(1))-y.dot(id.col(2))+id.col(3).sum())
    -2*n*covtermx.dot(covtermy);
  
  //cov2d equal is the square of the distance covariance between x and y 
  double cov2d = D/pow(n,2)-2*a_x.dot(b_y)/pow(n,3)+a_x.sum()*b_y.sum()/pow(n,4);
  
  return(cov2d);
}


//VectorXd FastDcov2d(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> y)
VectorXd FastDcov2d(const MatrixXd X,const VectorXd y)
{
  //Sort the data by x. This does not change the answer
  //a_x is the vector of row sums of distance matrix of x
  int n=y.size(),p=X.cols();
  VectorXi index(n);
  VectorXd s_y(n),si(n),a_x(n);
  
  sort_vec(y, s_y, index);
  
  for (int i = 0; i < n; i++){
    si(i)=s_y.head(i+1).sum();
  }
  
  a_x=2*VectorXd::LinSpaced(n, 1, n).array()-n;
  a_x =s_y.cwiseProduct(a_x);
  a_x=a_x.array()+(si(n-1)-2*si.array()).array();
  
  VectorXd dcov(p);
  for (int j=0; j<p; j++) {
    dcov(j)=DCov2d(X.col(j),s_y,a_x,index);
  }
  
  return(dcov);
}


// [[Rcpp::export]]
SEXP FastDcor2d(const Map<MatrixXd> X,const Map<VectorXd> y)
{
  //Sort the data by x. This does not change the answer
  //a_x is the vector of row sums of distance matrix of x
  int p=X.cols();// n=y.size(),
  VectorXd dcor(p);
  
  for (int j = 0; j < p; j++)
  {
    dcor(j)=*FastDcov2d(X.col(j),X.col(j)).data();
  }
  
  dcor=FastDcov2d(X,y).cwiseQuotient((dcor*FastDcov2d(y,y)).cwiseSqrt());
  return(wrap(dcor));
}

