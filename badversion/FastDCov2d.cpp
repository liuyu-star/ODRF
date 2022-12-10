#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace std;
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::wrap;
using Rcpp::Named;


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
 int i, n=x.size(); 
 MatrixXd it=MatrixXd::Zero(n+1,4);
 VectorXd y(n), xy(n); 
  
  y=z(index);
  xy=x.cwiseProduct(y);
  //VectorXi seq1n= VectorXi::LinSpaced(n, 1, n); 
  std::vector<int> seq1n= VectorXi::LinSpaced(n, 1, n); 

  //(Eigen::placeholders::all,ind) 
  it(seq1n,0)=VectorXd::Ones(n);
  it(seq1n,1)=y;
  it(seq1n,2)=x;
  it(seq1n,3)=xy;
  
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
    id.row(i)=it(syi.cwiseProduct(seqi),keep_cols).colwise().sum();
      }

   VectorXd b_y(n),by(n);
   by=2*VectorXd::LinSpaced(n, 1, n).array()-n;
   by =sy.cwiseProduct(by);
   b_y(idy) =by.array()+(si(n-1)-(2*si).array()).array();

    // the Frobenius inner product of the distance matrices 
    VectorXd covtermx=x.array()-x.mean();
    VectorXd covtermy=y.array()-y.mean();

    double D = 4*(xy.dot(id.col(0))-x.dot(id.col(1))-y.dot(id.col(2))+id.col(3).sum())
                 -2*n*covtermx.dot(covtermy);
      
    //cov2d equal is the square of the distance covariance between x and y 
    double cov2d = D/pow(n,2)-2*a_x.dot(b_y)/pow(n,3)+a_x.sum()*b_y.sum()/pow(n,4);

    return(cov2d);
}

// [[Rcpp::export]]
VectorXd FastDcov2d(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> y,int threads=10)
{
  #ifdef _OPENMP
  if ( threads > 0 )
    omp_set_num_threads( threads );
  //REprintf("Number of threads=%i\\n", omp_get_max_threads());
  #endif

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
#pragma omp parallel for schedule(dynamic)  
for (int j=0; j<p; j++) {
    dcov(j)=DCov2d(X.col(j),s_y,a_x,index);
}

return(dcov);
}


// [[Rcpp::export]]
SEXP FastDcov2d(const Map<MatrixXd> X,const Map<VectorXd> y,int threads=10)
{
#ifdef _OPENMP
  if ( threads > 0 )
    omp_set_num_threads( threads );
  REprintf("Number of threads=%i\\n", omp_get_max_threads());
#endif
  
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
#pragma omp parallel for schedule(dynamic)  
  for (int j=0; j<p; j++) {
    dcov(j)=DCov2d(X.col(j),s_y,a_x,index);
  }
  
  return(wrap(dcov));//dcov.cast<double>()
}


VectorXd FastDcor2d(const MatrixXd X,VectorXd y)
{
//Sort the data by x. This does not change the answer
//a_x is the vector of row sums of distance matrix of x
int n=y.size(),p=X.cols();
VectorXd dcor(p);

for (int j = 0; j < p; j++)
{
  dcor(j)=*FastDcov2d(X.col(j),X.col(j)).data();
}

dcor=FastDcov2d(X,y).cwiseQuotient((dcor*FastDcov2d(y,y)).cwiseSqrt());
return(dcor);
}


int main()
{
/*
VectorXd z(5), x(5), a_x(5);
VectorXi index(5);
//z << 5 , 2.0 , 7;
//x << 3.0 , 2.0 , 5.0;
//index << 1 , 0 , 2;
//a_x << 5.0 , 3.0 , 4.0;

z<<1,5,4,5,9;
x<<1,3,5,7,9;
index<<2,1,0,3,4;
a_x<<20,14,12,14,20;

double DCov= DCov2d(z,x,a_x,index);
cout << "return:\n" << DCov << endl;
*/

MatrixXd X(5,3);
VectorXd y(5);
X<<1,3,4,
   5,7,8,
   4,2,6,
   5,2,7,
   9,3,2;
y<<5,3,1,7,9;

cout << "return:\n" << FastDcov2d(X,y) << endl;
VectorXd dcor;
VectorXd dcor1(3);
dcor=FastDcov2d(X.col(0),X.col(0));//.cast<double>()
//dcor.addTo(dcor1);dcor1(0)=
dcor1(0)= *FastDcov2d(X.col(0),X.col(0)).data();
int data[15] = {1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15};
Matrix<int, 5, 3> mat2x2(data); 
//Matrix2i mat2x2(data);  
cout << "return:\n" << mat2x2 << endl;
cout << "return:\n" << FastDcor2d(X,y) << endl;

int threads=10;
 #ifdef _OPENMP
  if ( threads > 0 )
    omp_set_num_threads( threads );
   REprintf("Number of threads=%i\\n", omp_get_max_threads());
  #endif

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
#pragma omp parallel for schedule(dynamic)  
for (int j=0; j<p; j++) {
    dcov(j)=DCov2d(X.col(j),s_y,a_x,index);
}

cout << "return:\n" << dcov << endl;
}