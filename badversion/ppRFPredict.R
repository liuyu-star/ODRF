#' carPPRF Prediction
#'
#' Prediction a decision tree based on an input matrix and class vector.  This is the main function in the ppRF1 package.
#'
#' @param Xnew an n by d numeric matrix (preferable) or Xnew frame. The rows correspond to observations and columns correspond to features.
#' @return carPPtree Prediction v_subs
#'
#' @useDynLib ppRF
#' @import Rcpp
#' 
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

ppRFPredict= function(Xnew,ppForest,weight=FALSE){
#Returns the output of the ensemble (f_output) as well
#as a [num_treesXnum_samples] matrix (f_votes) containing
#the outputs of the individual trees. 
#
#The 'oobe' flag allows the out-of-bag error to be used to 
#weight the final response (only for classification).

if (!is.matrix(Xnew)) {
  Xnew <- as.matrix(Xnew)
}

n=nrow(Xnew)
nC=length(ppForest$Levels)
ntrees=length(ppForest$trees);

#Votes = matrix(0,ntrees,n);
#oobErr = rep(0,ntrees);
#for(i in 1:ntrees){
  #Votes[i,] = PPtreePredict(Xnew,ppForest$trees[[i]]);
#  oobErr[i] = ppForest$trees[[i]]$oobErr;
#}
#Votes=t(sapply(seq(ntrees), function(i)PPtreePredict(Xnew,ppForest$trees[[i]])))

oobErr=sapply(ppForest$trees,function(trees)trees$oobErr)
VALUE=rep(ifelse(ppForest$params$method=='r',0,as.character(0)),n)
Votes=t(vapply(ppForest$trees, function(trees){ppCARTPredict(Xnew,trees)},VALUE))

if ((!ppForest$params$storeOOB)&weight) {
  stop("out-of-bag indices for each tree are not stored. ppRF must be called with storeOOB = TRUE.")
}
weights=weight*oobErr+(!weight);
if(tolower(ppForest$params$method)%in%c('c','g')){
    #WEIGHTED_HIST Computes weighted histograms
    #   f_hist = weighted_hist(f_votes,weights,nlabels);
    #   f_votes is a 2-D matrix MxN where M is the number of votes and N the 
    #   number of samples.
    #   weights is a Mx1 weight vector
    #   nlabels is the number of labels
   
    # prob=matrix(0,n,nC)
     # for (i in 1:n) {
    #    prob[i,]=aggregate(c(rep(0,nC),weights[,i]), by=list(c(1:nC, f_votes[,i])),sum)[,2];
    #  }
    #f_votes[i,] =sapply(f_votes[i,],function(pred)which(ppForest$Levels%in%pred));
    #f_votes=as.numeric(as.factor(c(f_votes)))
    
    #Levels=levels(as.factor(c(ppForest$Levels[1:nC],Votes[,1])))
    #Levels=max.col(matrix(Levels,nC,nC)==matrix(ppForest$Levels,nC,nC,byrow = T))  
  #prob=matrix(0,n,nC);
  #for (i in seq(n)) {
  #  prob[i,]=aggregate(c(rep(0,nC),weights), by=list(c(1:nC,as.integer(as.factor(Votes[,i])))),sum)[,2];
  #}
  #prob[,Levels]=prob
  
    weights=rep(weights,n)
    Votes=factor(c(Votes),levels =ppForest$Levels)
    Votes=as.integer(Votes)+nC*rep(0:(n-1),rep(ntrees,n));
    Votes=aggregate(c(rep(0,n*nC),weights), by=list(c(1:(n*nC),Votes)),sum)[,2];
    
    prob=matrix(Votes,n,nC,byrow = TRUE);
    prob=prob/matrix(rowSums(prob),n,nC)
    colnames(prob)=ppForest$Levels
    
    #pred=apply(prob,1,which.max);
    pred=max.col(prob) ## "random"
    pred=ppForest$Levels[pred]
}else{
    prob=weights/sum(weights)
    pred=t(Votes)%*%prob
  #pred=colMeans(Votes);
  #prob=NULL
}

return(list(pred=pred,prob=prob))
}
