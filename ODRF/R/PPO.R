#' PP optimization using various projection pursuit indices
#' 
#' Find the q-dim optimal projection using various projectin pursuit indices 
#' with class information
#' @title Projection pursuit optimization
#' 
#' @param X data matrix without class information
#' @param y class information vector
#' @param q dimension of projection matrix
#' @param ppMethod method for projection pursuit; "LDA", "PDA", "Lr", "GINI", 
#'                 and "ENTROPY"
#' @param weight weight flag in LDA, PDA and Lr index
#' @param r r in Lr index
#' @param lambda lambda in PDA index
#' @param energy energy parameter
#' @param cooling cooling parameter
#' @param TOL tolerance
#' @param maxiter number of maximum iteration
#' 
#' @return bestIndex maximum LDA index value
#' @return projMat optimal q-dim projection matrix

#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @import Pursuit 
#' @export
#' @keywords projection pursuit
#' 
#' @examples
#' data(iris)
#' PP.proj.result <- ppOpt(as.matrix(iris[,1:4]),iris[,5])
#' PP.proj.result
## usage ppOpt(y,X,q=1,ppMethod="LDA",weight=TRUE,r=1,
#              lambda=0.1,energy=0,cooling=0.999,TOL=0.0001,maxiter = 50000)
#@import PPtreeViz
PPO <- function(x,y, q = 1L, ppMethod = "LDA", weight = TRUE, r = 1L, lambda = 0.1, 
                  energy = 0, cooling = 0.9, TOL = 0.0001, maxiter = 1000L,...) {
  x=as.matrix(x)
  p=ncol(x)
 
  if(!ppMethod%in%c("LDA","PDA","Lr","GINI","ENTROPY")){
    res <- PP_Optimizer(data = x, class = y, findex = ppMethod,
                        optmethod = "GTSA", dimproj = q, sphere = TRUE, 
                        weight = weight, lambda = lambda, r = r, cooling = cooling, 
                        eps = TOL, maxiter = maxiter, half = 30)#invisible()
    
    projbest=res$vector.opt
    indexbest=res$index[length(res$index)]
  }else{
    result <- ODRF:::ppOptCpp(as.factor(y),x, q, ppMethod, weight, r, lambda, energy, cooling, TOL, maxiter)
    indexbest=result$indexbest
    projbest=result$projbest
    #class(result)<-append(class(result),"ppOptim")
    #.Call(`_ppRF_ppOpt`, y, X, q, PPmethod, weight, r, lambda, energy, cooling, TOL, maxiter)
  }
  
  projMat=matrix(0,p,3)
  projMat[,1]=1:p
  projMat[,2]=rep(1,p)
  projMat[,3]=projbest
  colnames(projMat)=c("Varible","Number of projections","Projection coefficient")
  return(list(bestIndex=indexbest,projMat=projMat))
}
