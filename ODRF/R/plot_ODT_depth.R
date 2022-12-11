#' projection pursuit classification tree plot
#' 
#' Draw projection pursuit classification tree with tree structure. It is 
#' modified from a function in party library.
#' @title PPtree plot
#' @param x PPtreeclass object
#' @param font.size font size of plot
#' @param width.size size of eclipse in each node.
#' @param main main title
#' @param sub sub title
#' @param ... arguments to be passed to methods
#' @references Lee, EK(2017) 
#' PPtreeViz: An R Package for Visualizing Projection Pursuit Classification 
#' Trees, Journal of Statistical Software <doi:10.18637/jss.v083.i08>
#' @keywords tree
#' 
#' @export
#' 
#' @examples
#' data(iris)
#' Tree.result <- PPTreeclass(Species~., data = iris,"LDA")
#' Tree.result
#' plot(Tree.result,xjust=3)
plot_ODT_depth=function(formula,data,newdata,method='i-classification',NodeRotateFun="RotMatPPO",
                           paramList=NULL,main=paste0("Oblique ",
                           ifelse(method=="regression","Regression","Classification")," Tree"),...){
  set.seed(221109)
  
  paramList$formula=formula;paramList$data=data;paramList$MaxDepth=Inf
  paramList$method=method;paramList$NodeRotateFun=NodeRotateFun
  tree <- do.call(ODT, paramList)
  Depth=max(tree$structure$nodeDepth)
  #Depth=Depth+ceiling(Depth/2)
  
  y= data[,all.vars(tree$terms)[1]]
  ynew= newdata[,all.vars(tree$terms)[1]]
  Xnew= newdata[,all.vars(tree$terms)[-1]]
  Xnew=as.matrix(Xnew) 
  
  err=rep(0,Depth)
  for (d in 1:Depth) {
    set.seed(221109)
    paramList$MaxDepth=d
    tree <- do.call(ODT, paramList)
    pred <- predict(tree,Xnew)
    
    if(method!="regression"){
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
