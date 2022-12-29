#' plot oblique decision tree depth
#' 
#' Draw the error graph of class \code{ODT} at different depths.
#' 
#' @param formula Object of class \code{formula} with a response but no interaction terms describing the model to fit. If this is a data frame, it is taken as the model frame. (see \code{\link{model.frame}})
#' @param data Training data of class \code{data.frame} in \code{\link{ODRF}} is used to calculate the OOB error.
#' @param newdata A data frame or matrix containing new data is used to calculate the test error. If it is missing, let it be \code{data}.
#' @param type The criterion used for splitting the nodes. 'i-classification': information gain and 'g-classification': gini impurity index for classification; 'regression': mean square error for regression. 
#' 'auto' (default): If the response in \code{data} is a factor, 'g-classification' is used, otherwise regression is assumed.
#' @param NodeRotateFun Name of the function of class \code{character} that implements a linear combination of predictors in the split node. including "RotMatPPO": see \code{\link{RotMatPPO}} (default, model="PPR"),
#' "RotMatRand": see \code{\link{RotMatRand}}, and "RotMatRF": see \code{\link{RotMatRF}}. Users can define this function, for details see \code{\link{RotMatMake}}.
#' @param paramList List of parameters used by the functions \code{NodeRotateFun}. If left unchanged, default values will be used, for details see \code{\link[ODRF]{defaults}}.
#' @param digits integer indicating the number of decimal places (round) or significant digits (signif) to be used. 
#' @param main main title
#' @param ... arguments to be passed to methods.
#' 
#' @return test error of \code{newdata}, misclassification rate (MR) for classification or mean square error (MSE) for regression.
#' 
#' @keywords tree
#' 
#' @seealso \code{\link{ODT}} \code{\link{plot.ODT}}
#' 
#' @examples
#' data(body_fat)
#' set.seed(221212)
#' train = sample(1:252,100)
#' train_data = data.frame(body_fat[train,])
#' test_data = data.frame(body_fat[-train,])
#' plot_ODT_depth(Density~.,train_data,test_data,type='regression')
#' 
#' @export
plot_ODT_depth=function(formula,data=NULL,newdata=NULL,type='i-classification',NodeRotateFun="RotMatPPO",
                        paramList=NULL,digits=NULL,main=paste0("Oblique ",
                                                   ifelse(type=="regression","Regression","Classification")," Tree"),...){
  set.seed(221109)
  if(is.null(data)){
    data <- data.frame(y=eval(formula[[2]]),eval(formula[[3]]))
    formula=y~.
  }
  
  paramList$formula=formula;paramList$data=data;paramList$MaxDepth=Inf
  paramList$type=type;paramList$NodeRotateFun=NodeRotateFun
  tree <- do.call(ODT.formula, paramList)
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
    tree <- do.call(ODT.formula, paramList)
    pred <- predict(tree,Xnew)
    
    if(type!="regression"){
      err[d]<- mean(pred != ynew)
    }else{
      err[d]<- mean((pred-ynew)^2)#/mean((ynew-mean(y))^2);
    }
  }
  
  #strerr=strsplit(as.character(min(err)),split="")[[1]]
  #errid=which(strerr[-seq(which(strerr=="."))]!="0")[2]
  #err=round(err,errid)
  minErr=strsplit(as.character(min(err)),"")[[1]]
  id=which(minErr=="e")
  if(type!="regression"){
    digits=0
  }else if(is.null(digits)){
    if(length(id)>0){
      digits=sum(as.numeric(paste0(minErr[c(id+2,length(minErr))]))*c(10,1))
    }else{
      digits=which(minErr[-seq(which(minErr=="."))]!=0)[2]
    }
  }
  
  #type = "l",lty=1,
  matplot(1:Depth, err, pch = 21, bg = "skyblue", type = "b",lty=1, xlab="Depth", ylab=paste0("Error (*",10^-digits,")"), col = c("black"),main=main,xaxt="n",yaxt="n")
  axis(1, seq(1,Depth,length.out = min(6,Depth)),round(seq(1,Depth,length.out = min(6,Depth))),cex.lab = 1.5,cex.axis = 1.25)
  axis(2, seq(err[1],err[Depth],length.out = min(6,Depth)),round(seq(err[1],err[Depth],length.out = min(6,Depth))*10^digits,2),cex.lab = 1.5,cex.axis = 1.25)
  
  Error=cbind(Depth=1:Depth,Error=err)
  return(invisible(Error))
}
