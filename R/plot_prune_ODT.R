#' to plot pruned oblique decision tree
#'
#' Plot the error graph of the pruned oblique decision tree at different split nodes.
#'
#' @param x An object of class \code{\link{prune.ODT}}.
#' @param position Position of the curve label, including "topleft" (default), "bottomright", "bottom", "bottomleft", "left", "top", "topright", "right" and "center".
#' @param digits Integer indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @param main main title
#' @param ... Arguments to be passed to methods.
#'
#' @return The leftmost value of the horizontal axis indicates the tree without pruning, while the rightmost value indicates the data without splitting and using the average value as the predicted value.
#'
#' @seealso \code{\link{ODT}} \code{\link{prune.ODT}}
#'
#' @examples
#' data(body_fat)
#' set.seed(221212)
#' train <- sample(1:252, 100)
#' train_data <- data.frame(body_fat[train, ])
#' test_data <- data.frame(body_fat[-train, ])
#'
#' tree <- ODT(Density ~ ., train_data, split = "mse")
#' prune_tree <- prune(tree, test_data[, -1], test_data[, 1])
#' # Plot pruned oblique decision tree structure (default)
#' plot(prune_tree)
#' # Plot the error graph of the pruned oblique decision tree.
#' class(prune_tree) <- "prune.ODT"
#' plot(prune_tree)
#'
#' @keywords tree plot prune
#' @rdname plot.prune.ODT
#' @aliases plot.prune.ODT
#' @method plot prune.ODT
#' @export
plot.prune.ODT <- function(x, position = "topleft", digits = NULL, main = NULL, ...) {
  ppTree <- x
  rm(x)
  pruneError <- ppTree$pruneError

  if (is.null(main)) {
    main <- paste0("Oblique ", ifelse(ppTree$split == "mse", "Regression", "Classification"), " Tree")
  }

  # par(mfrow = c(1,2))
  op <- par(plt = c(0.2, 0.8, 0.2, 0.90))
  # par(plt = c(0.07, 0.93, 0.1, 0.90))
  # par(adj=0.5)
  minLen <- min(6, length(pruneError[, 1]))
  x <- seq_len(nrow(pruneError))
  minErr <- strsplit(as.character(min(pruneError[, 4])), "")[[1]]
  id <- which(minErr == "e")
  if (ppTree$split != "mse") {
    digits <- 0
  } else if (is.null(digits)) {
    if (length(id) > 0) {
      digits <- sum(as.numeric(paste0(minErr[c(id + 2, length(minErr))])) * c(10, 1))
    } else {
      digits <- which(minErr[-seq(which(minErr == "."))] != 0)[2]
    }
  }

  if(digits==0){
    ylab = paste0("Error")
  }else if(digits==2){
    ylab = paste0("Error (%)")
  }else{
    ylab = substitute(paste("Error ("*10^{-dig},")"),list(dig = digits))
  }

  plot(x, pruneError[, 4], pch = 21, bg = "skyblue", type = "b", lty = 1, xlab = "Split node", ylab = ylab, main = main, xaxt = "n", yaxt = "n") # , col = c("black")
  # plot(x, pruneError[,4],pch = 21, bg = "skyblue", type = "p",lty=1, xlab="The number of split nodes", ylab="Error",main=main,xaxt="n",yaxt="n")#, col = c("black")
  axis(1, seq(min(x), max(x), length.out = minLen), round(seq(max(pruneError[, 1]), min(pruneError[, 1]), length.out = minLen)), cex.lab = 1.5, cex.axis = 1.25)
  axis(2, seq(min(pruneError[, 4]), max(pruneError[, 4]), length.out = minLen), round(seq(min(pruneError[, 4]), max(pruneError[, 4]), length.out = minLen) * 10^digits, 2), cex.lab = 1.5, cex.axis = 1.25)
  abline(h = pruneError[1, 4], lwd = 1.5, lty = 2, col = "red")
  on.exit(par(op))

  op <- par(new = T)
  plot(x, pruneError[, 3], pch = 4, type = "p", lty = 3, xaxt = "n", yaxt = "n", ann = F, axes = F) # col = c("red"),
  abline(h = round(seq(1, max(pruneError[, 3]), length.out = minLen)), lwd = 1.5, lty = 2, col = "gray")
  # axis(1, seq(1,max(pruneError[,1]),length.out = minLen),round(seq(1,max(pruneError[,1]),length.out = minLen)),cex.lab = 1.5,cex.axis = 1.25)
  axis(4, round(seq(1, max(pruneError[, 3]), length.out = minLen)), cex.lab = 1.5, cex.axis = 1.25)

  # axis(side = 4)
  mtext("Depth", side = 4, line = 3)

  legend(x = position, legend = c("Error", "Depth"), lty = c(1, 2), pch = c(21, 4), pt.bg = c("skyblue", "black"), col = c("black", "black"), bty = "n")
  on.exit(par(op))

  return(invisible(pruneError))
}
