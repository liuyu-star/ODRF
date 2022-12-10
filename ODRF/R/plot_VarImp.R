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
#' @aliases plot.VarImp
#' @rdname plot.VarImp
#' @method plot VarImp
#' @export
#' 
#' @examples
#' data(iris)
#' Tree.result <- PPTreeclass(Species~., data = iris,"LDA")
#' Tree.result
#' plot(Tree.result,xjust=3)
plot.VarImp <- function(varImp,nvar=min(30, nrow(varImp$varImp)), main=paste0("Oblique ",
                           ifelse(varImp$method=="regression","Regression","Classification")," Forest"), ...) {
  imp=varImp$varImp
  imp=imp[1:nvar,,drop = FALSE]
  ## If there are more than two columns, just use the last two columns.
  op <- par(xaxs = "i")
  dotchart(imp[,2], xlab="Increased error (%)", ylab="",main=main,xaxt="n",
           cex.lab = 1.5,cex.axis = 1.25, bg = "skyblue")
  axis(1, seq(min(imp[,2]),max(imp[,2]),length.out = min(6,nvar)),
       round(seq(min(imp[,2]),max(imp[,2]),length.out = min(6,nvar))*100,ifelse(varImp$method=="regression",0,2)),
       cex.lab = 1.5,cex.axis = 1.25)
  par(op)
  
  return(invisible(imp))
}
