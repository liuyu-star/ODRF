#' prune oblique decision tree plot
#' 
#' Draw the error graph of class \code{ODT} at different number of split nodes.
#' 
#' @param ppTree an object of class \code{\link{prune.ODT}}.
#' @param position Position of the curve label.
#' @param main main title
#' @param ... arguments to be passed to methods.
#' 
#' @return Error of test data after each pruning, classification error rate for classification or RPE(MSE/mean((ytest-mean(y))^2)) for regression.
#' 
#' @keywords tree
#' 
#' @seealso \code{\link{ODT}} \code{\link{prune.ODT}}
#' 
#' @examples
#' data(seeds)
#' set.seed(221212)
#' train = sample(1:209,100)
#' train_data = data.frame(seeds[train,])
#' test_data = data.frame(seeds[-train,])
#' 
#' tree = ODT(varieties_of_wheat~.,train_data,type='i-classification')
#' tree = prune(tree,train_data)
#' #oblique decision tree plot (default)
#' plot(tree)
#' #prune oblique decision tree plot
#' class(tree)="prune.ODT"
#' plot(tree)
#' 
#' @aliases plot.prune.ODT
#' @rdname plot.prune.ODT
#' @method plot prune.ODT
#' @export
plot.prune.ODT=function(ppTree,position="topleft",main=paste0("Oblique ",
                        ifelse(ppTree$type=="regression","Regression","Classification")," Tree"),...){
  pruneError=ppTree$pruneError
  
  #par(mfrow = c(1,2))
  #par(plt = c(0.12, 0.88, 0.25, 0.90))
  par(plt = c(0.07, 0.93, 0.1, 0.90))
  #par(adj=0.5)
  minLen=min(6,length(pruneError[,1]))
  x=seq(nrow(pruneError))
  
  plot(x, pruneError[,4],pch = 21, bg = "skyblue", type = "b",lty=1, xlab="The number of split nodes", ylab="Error",main=main,xaxt="n",yaxt="n")#, col = c("black")
  #plot(x, pruneError[,4],pch = 21, bg = "skyblue", type = "p",lty=1, xlab="The number of split nodes", ylab="Error",main=main,xaxt="n",yaxt="n")#, col = c("black")
  axis(1, seq(min(x),max(x),length.out = minLen),round(seq(max(pruneError[,1]),min(pruneError[,1]),length.out = minLen)),cex.lab = 1.5,cex.axis = 1.25)
  axis(2, seq(min(pruneError[,4]),max(pruneError[,4]),length.out = minLen),round(seq(min(pruneError[,4]),max(pruneError[,4]),length.out = minLen),2),cex.lab = 1.5,cex.axis = 1.25)
  abline(h=pruneError[1,4],lwd=1.5,lty=2,col="red")
  
  par(new = T) 
  plot(x, pruneError[,3],pch = 4,type = "p",lty=3,xaxt="n",yaxt="n",ann = F, axes = F)# col = c("red"),
  abline(h=round(seq(1,max(pruneError[,3]),length.out = minLen)),lwd=1.5,lty=2,col="gray")
  #axis(1, seq(1,max(pruneError[,1]),length.out = minLen),round(seq(1,max(pruneError[,1]),length.out = minLen)),cex.lab = 1.5,cex.axis = 1.25)
  axis(4,round(seq(1,max(pruneError[,3]),length.out = minLen)),cex.lab = 1.5,cex.axis = 1.25)
  
  #axis(side = 4)
  mtext("Depth", side = 4, line = 3)
  
  legend(x=position, legend = c("Error", "Depth"),lty = c(1,2),pch = c(21,4), pt.bg = c("skyblue","black") ,col = c("black","black"),bty="n")
  
 
  return(invisible(pruneError))
}
