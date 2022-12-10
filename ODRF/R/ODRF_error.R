#' carPPRF Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param Xnew an n by d numeric matrix (preferable) or Xnew frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction v_subs PPtree
#'
#' @import Rcpp
#' @export
#' 
#' @examples
#' ### Train RerF on numeric Xnew ###
#' library(rerf)
#' forest <- RerF1(as.matrix(iris[, 1:4]), iris[[5L]], num.cores = 1L)
#'
#' ### Train RerF on one-of-K encoded categorical Xnew ###
#' df1 <- as.Xnew.frame(Titanic)
#' nc <- ncol(df1)
#' df2 <- df1[NULL, -nc]
#' for (i in which(df1$Freq != 0L)) {
#'   df2 <- rbind(df2, df1[rep(i, df1$Freq[i]), -nc])
#' }
## @aliases predict #' @method predict ppRF #' @aliases predict_ppRF #' @method predict ppRF
ODRF.error <- function(forest,y, Xnew=NULL, ynew=NULL, ...) {
  if(forest$method!="regression"){
    y <- factor(y,levels = forest$Levels)
  }
  n=length(y)
  nt=ntrees=forest$forest$ntrees
  nC=length(forest$Levels)
  ny=length(ynew)
  
  treeVotes=predict(forest,Xnew,weight=FALSE)$TreePrediction
  err.test=rep(0,ntrees)
  if(forest$method=="regression"){
    e.0=mean((ynew-mean(y))^2)
    pred=rowSums(treeVotes);
    err.test[nt]=mean((ynew-pred/nt)^2)/e.0;
    for (t in seq(nt-1,1)) {
      pred=pred-treeVotes[,t+1]
      err.test[t]=mean((ynew-pred/t)^2)/e.0;
    }
  }else{
    weights=rep(1,ny*nt)
    Votes=factor(c(t(treeVotes)),levels =forest$Levels)
    treeVotes=matrix(as.integer(Votes),nt,ny)
    
    Votes=c(treeVotes)+nC*rep(0:(ny-1),rep(nt,ny));
    Votes=aggregate(c(rep(0,ny*nC),weights), by=list(c(1:(ny*nC),Votes)),sum)[,2];
    #Votes=aggregate(c(rep(0,ny*nC),weights), by=list(c(1:(ny*nC),Votes)),cumsum)[,2];
    
    #prob=matrix(Votes,ny,nC,byrow = TRUE);
    Votes=matrix(Votes,ny,nC,byrow = TRUE);
    pred=forest$Levels[max.col(Votes)]## "random"
    err.test[nt]=mean(ynew!=pred)
    treeC=matrix(seq(nC),ny,nC,byrow = TRUE)
    for (t in seq(nt-1,1)) {
      Votes=Votes-(treeC==treeVotes[t+1,])*1
      #pred=apply(prob,1,which.max);
      pred=forest$Levels[max.col(Votes)]## "random"
      err.test[t]=mean(ynew!=pred)
    }
  }
    
  
  err.oob=rep(0,ntrees)
  for (tt in 1:ntrees) {
    oobVotes=matrix(NA,n,tt)
    for (t in 1:tt) {
      oobVotes[forest$ppTrees[[t]]$oobIndex,t]=forest$ppTrees[[t]]$oobPred
    }
    idx=which(rowSums(is.na(oobVotes))<tt) 
    oobVotes=oobVotes[idx,,drop = FALSE]
    
    if(forest$method=="regression"){
      pred=rowMeans(oobVotes,na.rm = TRUE);
      err=mean((y[idx]-pred)^2)/mean((y[idx]-mean(y))^2);
    }else{
      ny=length(y[idx])
      nt=ncol(oobVotes)
      weights=rep(1,ny*nt)
      Votes=factor(c(t(oobVotes)),levels =forest$Levels)
      Votes=as.integer(Votes)+nC*rep(0:(ny-1),rep(nt,ny));
      Votes=aggregate(c(rep(0,ny*nC),weights), by=list(c(1:(ny*nC),Votes)),sum)[,2];
      
      prob=matrix(Votes,ny,nC,byrow = TRUE);
      #pred=apply(prob,1,which.max);
      pred=max.col(prob) ## "random"
      pred=forest$Levels[pred]
      err=mean(y[idx]!=pred)
    }
    err.oob[tt]=err
  }
  
  error=list(err.oob=err.oob,err.test=err.test,method=forest$method)
  
  class(error)="ODRF.error"
  return(error)
}


