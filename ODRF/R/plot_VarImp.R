#' Variable Importance Plot
#' 
#' Dotchart of variable importance as measured by a Oblique Decision Random Forest.
#' 
#' @param varImp an object of class \code{\link{VarImp}}.
#' @param nvar How many variables to show?
#' @param main plot title.
#' @param ... arguments to be passed to methods.
#' 
#' @return A matrix of importance measure, first column for each predictor variable and second column is Increased error. 
#' classification error rate for classification or RPE(MSE/mean((ytest-mean(y))^2)) for regression.
#' 
#' @seealso \code{\link{ODRF}} \code{\link{VarImp}}
#' 
#' @examples
#' data(breast_cancer)
#' forest = ODRF(diagnosis~.,seeds,type='i-classification')
#' varimp=VarImp(forest,seeds)
#' plot(varimp)
#' 
#' @aliases plot.VarImp
#' @rdname plot.VarImp
#' @method plot VarImp
#' @export
plot.VarImp <- function(varImp,nvar=min(30, nrow(varImp$varImp)), main=paste0("Oblique ",
                           ifelse(varImp$type=="regression","Regression","Classification")," Forest"), ...) {
  imp=varImp$varImp
  imp=imp[1:nvar,,drop = FALSE]
  ## If there are more than two columns, just use the last two columns.
  op <- par(xaxs = "i")
  dotchart(imp[,2], xlab="Increased error (%)", ylab="",main=main,xaxt="n",
           cex.lab = 1.5,cex.axis = 1.25, bg = "skyblue")
  axis(1, seq(min(imp[,2]),max(imp[,2]),length.out = min(6,nvar)),
       round(seq(min(imp[,2]),max(imp[,2]),length.out = min(6,nvar))*100,ifelse(varImp$type=="regression",0,2)),
       cex.lab = 1.5,cex.axis = 1.25)
  par(op)
  
  return(invisible(imp))
}
