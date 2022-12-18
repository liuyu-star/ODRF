#' Projection Pursuit Optimization
#' 
#' Find the optimal projection using various projectin pursuit indices with class information.
#' 
#' @param X An n by d numeric matrix (preferable) or data frame.
#' @param y a n vector.
#' @param model \itemize{model for projection pursuit.
#' \item{"PPR"(default): projection projection regression from \code{\link{ppr}}. When y is a category label, it is
#' expanded to K binary features.}
#' \item{"Log": logistic regression from \code{\link{nnet}}.}
#' \item{"Rand": The random projection generated from {-1, 1}. 
#' The following models can only be used for classification, i.e. the type must be 'i-classification' or 'g-classification'.}
#' \item{"LDA", "PDA", "Lr", "GINI", and "ENTROPY" from library \code{PPtreeViz}.}
#' \item{\itemize{The following models from library \code{\link{Pursuit}}.
#' \item{"holes": Holes index}
#' \item{"cm": Central Mass index}
#' \item{"holes": Holes index}
#' \item{"friedmantukey": Friedman Tukey index}
#' \item{"legendre": Legendre index}
#' \item{"laguerrefourier": Laguerre Fourier index}
#' \item{"hermite": Hermite index,}
#' \item{"naturalhermite": Natural Hermite index}
#' \item{"kurtosismax": Maximum kurtosis index,}
#' \item{"kurtosismin": Minimum kurtosis index,}
#' \item{"moment": Moment index}
#' \item{"mf": MF index}
#' \item{"chi": Chi-square index}
#' }}
#' }
#' @param type The criterion used for splitting the variable. 'g-classification': gini impurity index (classification, default), 
#'        'i-classification': information gain (classification) or 'regression': mean square error (regression).
#' @param weights A vector of length same as \code{data} that are positive weights. (default NULL)
#' 
#' @return optimal projection direction.
#' 
#' @keywords projection pursuit
#' 
#' @references Friedman, J. H. and Stuetzle, W. (1981). Projection pursuit regression. Journal of the American Statistical Association, 76, 817–823. doi: 10.2307/2287576.
#' @references Ripley, B. D. (1996) Pattern Recognition and Neural Networks. Cambridge.
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) PPtree: Projection Pursuit Classification Tree, Electronic Journal of Statistics, 7:1369-1386.
#' @references COOK, D., LEE, E. K., BUJA, A., WICKHAM, H.. Grand tours, projection pursuit guided tours and manual controls. In Chen, Chunhouh, Hardle, Wolfgang, Unwin, e Antony (Eds.), Handbook of data Visualization, Springer Handbooks of Computational Statistics, chapter III.2, p. 295-314. Springer, 2008.
#' 
#' @seealso \code{\link{RotMatMake}} \code{\link{RotMatRand}} \code{\link{RotMatRF}} \code{\link{RotMatMake}}
#' 
#' @examples
#' #classification
#' data(iris)
#' (PP <- PPO(iris[,1:4],iris[,5],model = "Log",type='i-classification'))
#' (PP <- PPO(iris[,1:4],iris[,5],model = "PPR",type='i-classification'))
#' (PP <- PPO(iris[,1:4],iris[,5],model = "LDA",type='i-classification'))
#' 
#' #regression
#' data(body_fat)
#' (PP <- PPO(body_fat[,2:15],body_fat[,1],model = "Log",type='regression'))
#' (PP <- PPO(body_fat[,2:15],body_fat[,1],model = "Rand",type='regression'))
#' (PP <- PPO(body_fat[,2:15],body_fat[,1],model = "PPR",type='regression'))
#' 
#' @import Pursuit Rcpp 
#' @importFrom stats ppr
#' @importFrom nnet nnet
#' @export
PPO <- function(X,y,model = "PPR",type='i-classification',weights=NULL,...) {
  X=as.matrix(X)
  p=ncol(X)
  n=length(y)
  weights=if(is.null(weights))rep(1, n)
  
  Y=c(y);indC=0L
  if(type!='regression'){
    y=as.factor(y)
    indC=levels(y)
    if(length(indC)>2){
      Y=(matrix(y,n,length(indC))==matrix(indC,n,length(indC),byrow = TRUE))+0
    }else{
      Y=as.integer(y)
    }
  }
 
  if(model=="PPR"){
    PP <- ppr(X, Y,weights=weights,nterms = 1, bass=1)$alpha# sm.method="spline",sm.method="gcvspline"
  }else if(model=="Rand"){
    PP <- sample(c(1L, -1L),p, replace = TRUE,prob = c(0.5, 0.5))
  }else if(model=="Log"){
    #if((n>5*p)&(p<10)){
     # PP <- ppr(X, Y,weights,nterms = 1, bass=1)$alpha
    #}else{
      #PP = try(nnet(Xi, y, size=1,trace=FALSE)$wts[2:(1+pi)], silent = TRUE)
      PP = nnet(y~X,weights=weights, size=1,linout=TRUE,trace=FALSE)$wts[2:(1+p)]#  
    #}
  }else if(model%in%c("LDA","PDA","Lr","GINI","ENTROPY")){
    PP <- ODRF:::ppOptCpp(y, X, q=1, PPmethod=model, weight=TRUE, r=1, lambda=0.1, energy=0, cooling=0.9, TOL=0.0001, maxiter=1000L)$projbest
  }else{
    PP <- PP_Optimizer(data = X, class = y, findex = model,
                        optmethod = "GTSA", dimproj = 1, sphere = TRUE, 
                        weight = TRUE, lambda = 0.1, r = 1, cooling = 0.9, 
                        eps = 0.0001, maxiter = 1000L, half = 30)$vector.opt#invisible()
    #projbest=res$vector.opt
    #indexbest=res$index[length(res$index)]
  }
  
  #projMat=matrix(0,p,3)
  #projMat[,1]=1:p
  #projMat[,2]=rep(1,p)
  #projMat[,3]=projbest
  #colnames(projMat)=c("Varible","Number of projections","Projection coefficient")
  #return(list(projbest=projbest,indexbest=indexbest))
  return(as.vector(PP))
}
