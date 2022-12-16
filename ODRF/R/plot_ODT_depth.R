#' oblique decision tree depth plot
#' 
#' Draw the error graph of class \code{ODT} at different depths.
#' 
#' @param formula Object of class \code{formula} with a response but no interaction terms describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} in which to interpret the variables named in the formula.If data is missing it is obtained from the current environment by \code{formula}.
#' @param newdata A data frame or matrix containing new data.
#' @param type The criterion used for splitting the nodes. g-classification': gini impurity index(default) and i-classification': information gain for classification; 'regression': mean square error for regression.
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictor variables in the split node. Default is "RotMatPPO" with model="PPR". (see \code{\link{RotMatPPO}}) 
#' Users can define this function, for details see \code{\link{RotMatMake}}.
#' @param paramList Parameters in a named list to be used by \code{ODT}.
#' @param main main title
#' @param ... arguments to be passed to methods.
#' 
#' @return OOB error and test error, classification error rate for classification or RPE(MSE/mean((ytest-mean(y))^2)) for regression.
#' 
#' @keywords tree
#' 
#' @seealso \code{ODT} \code{plot.ODT}
#' 
#' @examples
#' library(ODRF)
#' 
#' data(body_fat)
#' set.seed(221212)
#' train = sample(1:252,100)
#' train_data = data.frame(body_fat[train,])
#' test_data = data.frame(body_fat[-train,])
#' plot_ODT_depth(Density~.,train_data,test_data,type='regression')
#' 
#' @rdname plot.ODT.depth
#' @export
plot_ODT_depth=function(formula,data=NULL,newdata=NULL,type='i-classification',NodeRotateFun="RotMatPPO",
                        paramList=NULL,main=paste0("Oblique ",
                                                   ifelse(type=="regression","Regression","Classification")," Tree"),...){
  set.seed(221109)
  if(is.null(data)){
    data <- data.frame(y=eval(formula[[2]]),eval(formula[[3]]))
    formula=y~.
  }
  
  paramList$formula=formula;paramList$data=data;paramList$MaxDepth=Inf
  paramList$type=type;paramList$NodeRotateFun=NodeRotateFun
  tree <- do.call(ODT, paramList)
  Depth=max(tree$structure$nodeDepth)
  #Depth=Depth+ceiling(Depth/2)
  
  vars=all.vars(tree$terms)
  y= data[,setdiff(colnames(data),vars[-1])]
  if(is.null(newdata)){
    newdata=data
  }
  ynew= newdata[,setdiff(colnames(newdata),vars[-1])]
  Xnew= newdata[,vars[-1]]
  Xnew=as.matrix(Xnew)
  
  err=rep(0,Depth)
  for (d in 1:Depth) {
    set.seed(221109)
    paramList$MaxDepth=d
    tree <- do.call(ODT, paramList)
    pred <- predict(tree,Xnew)
    
    if(type!="regression"){
      err[d]<- mean(pred != ynew)
    }else{
      err[d]<- mean((pred-ynew)^2)/mean((ynew-mean(y))^2);
    }
  }
  
  #strerr=strsplit(as.character(min(err)),split="")[[1]]
  #errid=which(strerr[-seq(which(strerr=="."))]!="0")[2]
  #err=round(err,errid)
  
  matplot(1:Depth, err, type = "l",lty=1, xlab="Depth", ylab="Error", col = c("black"),main=main,xaxt="n",yaxt="n")
  axis(1, seq(1,Depth,length.out = min(6,Depth)),round(seq(1,Depth,length.out = min(6,Depth))),cex.lab = 1.5,cex.axis = 1.25)
  axis(2, seq(err[1],err[Depth],length.out = min(6,Depth)),round(seq(err[1],err[Depth],length.out = min(6,Depth)),2),cex.lab = 1.5,cex.axis = 1.25)
  
  Error=cbind(Depth=1:Depth,Error=err)
  return(invisible(Error))
}
