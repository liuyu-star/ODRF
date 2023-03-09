#' Variable Importance Plot
#'
#' Dotchart of variable importance as measured by an Oblique Decision Random Forest.
#'
#' @param x An object of class \code{\link{VarImp}}.
#' @param nvar number of variables to show.
#' @param digits Integer indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @param main plot title.
#' @param ... Arguments to be passed to methods.
#'
#' @return The horizontal axis is the increased error of ODRF after replacing the variable, the larger the increased error the more important the variable is.
#'
#' @seealso \code{\link{ODRF}} \code{\link{VarImp}}
#'
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train <- sample(1:569, 200)
#' train_data <- data.frame(breast_cancer[train, -1])
#' forest <- ODRF(train_data[, -1], train_data[, 1], split = "gini",
#'   parallel = FALSE)
#' varimp <- VarImp(forest, train_data[, -1], train_data[, 1])
#' \donttest{
#' plot(varimp)
#' }
#' @keywords forest plot
#' @rdname plot.VarImp
#' @aliases plot.VarImp
#' @method plot VarImp
#' @export
plot.VarImp <- function(x, nvar = min(30,nrow(x$varImp)), digits = NULL, main = NULL, ...) {
  imp <- x$varImp
  imp <- imp[1:nvar, , drop = FALSE]

  if (is.null(main)) {
    main <- paste0("Oblique ", ifelse(x$split == "mse", "Regression", "Classification"), " Forest")
  }

  minErr <- strsplit(as.character(min(imp[, 2])), "")[[1]]
  id <- which(minErr == "e")
  if (x$split != "mse") {
    digits <- 0
  } else if (is.null(digits)) {
    if (length(id) > 0) {
      digits <- sum(as.numeric(paste0(minErr[c(id + 2, length(minErr))])) * c(10, 1))
    } else {
      digits <- which(minErr[-seq(which(minErr == "."))] != 0)[2]
    }
  }

  if(digits==0){
    xlab = paste0("Increased error")
  }else if(digits==2){
    xlab = paste0("Increased error (%)")
  }else{
    xlab = substitute(paste("Increased error (*",10^{-dig},")"),list(dig = digits))
  }

  ## If there are more than two columns, just use the last two columns.
  #op <- par(xaxs = "i") #* 10^digits


  dotchart(sort(imp[, 2]),
    xlab = xlab, ylab = "", main = main, xaxt = "n",
    cex.lab = 1.5, cex.axis = 1.25, bg = "skyblue"
  )
  axis(1, seq(min(imp[, 2]), max(imp[, 2]), length.out = min(6, nvar)),
    round(seq(min(imp[, 2]), max(imp[, 2]), length.out = min(6, nvar)) * 10^digits, 2),
    cex.lab = 1.5, cex.axis = 1.25
  )
  #par(op)

  return(invisible(imp))
}
