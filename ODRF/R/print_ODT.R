#' print ODT result
#'
#' Print the oblique decision tree structure.
#' @param x An object of class \code{\link{ODT}}.
#' @param projection Print projection coefficients in each node if TRUE.
#' @param cutvalue Print cutoff values in each node if TRUE.
#' @param verbose Print if TRUE, no output if FALSE.
#' @param ... Arguments to be passed to methods.
#'
#' @references Lee, EK(2017)
#' PPtreeViz: An R Package for Visualizing Projection Pursuit Classification
#' Trees, Journal of Statistical Software <doi:10.18637/jss.v083.i08>
#'
#' @seealso \code{\link{ODT}}
#'
#' @examples
#' data(iris)
#' tree <- ODT(Species ~ ., data = iris)
#' tree
#' print(tree, projection = TRUE, cutvalue = TRUE)
#'
#' @keywords tree
#' @rdname print.ODT
#' @aliases print.ODT
#' @method print ODT
#' @export
print.ODT <- function(x, projection = FALSE, cutvalue = FALSE, verbose = TRUE, ...) {
  ppTree <- x
  rm(x)
  numNode <- length(ppTree$structure$nodeCutValue)
  cutNode <- which(ppTree$structure$nodeCutValue != 0)

  TS <- matrix(0, numNode, 5)
  TS[, 1] <- seq(numNode)
  TS[, 2] <- ppTree[["structure"]][["childNode"]]
  if (ppTree$type != "regression") {
    TS[setdiff(seq(numNode), cutNode), 3] <- max.col(ppTree$structure$nodeNumLabel)[setdiff(seq(numNode), cutNode)]
  } else {
    TS[setdiff(seq(numNode), cutNode), 3] <- round(ppTree$structure$nodeNumLabel[, 1][setdiff(seq(numNode), cutNode)], 3)
  }
  TS[cutNode, 3] <- TS[cutNode, 2] + 1
  TS[cutNode, 4] <- seq_along(cutNode)
  TS[cutNode, 5] <- ppTree[["structure"]][["nodeCutIndex"]][cutNode]
  colnames(TS) <- c("node", "left_node", "right_node/leaf_label", "cut_node", "cut_node_index")
  leaf <- rep(0, numNode)
  leaf[setdiff(seq(numNode), cutNode)] <- seq(numNode - length(cutNode))

  # TS<-ppTree$ppTree.Struct
  # Alpha<-ppTree$projbest.node
  nodeRotaMat <- ppTree[["structure"]][["nodeRotaMat"]]
  Alpha <- matrix(0, length(cutNode), ppTree[["data"]][["p"]])
  for (cn in seq_along(cutNode)) {
    idx <- which(nodeRotaMat[, 2] == cutNode[cn])
    Alpha[cn, nodeRotaMat[idx, 1]] <- nodeRotaMat[idx, 3]
  }

  CutValue <- ppTree$structure$nodeCutValue[cutNode]
  # CutValue<-ppTree$splitCutoff.node
  # gName<-ppTree$Levels
  # gName<-names(table(ppTree$origclass))
  pastemake <- function(k, arg, sep.arg = "") {
    temp <- ""
    for (i in 1:k) {
      temp <- paste(temp, arg, sep = sep.arg)
    }
    return(temp)
  }
  TreePrint <- "1) root"
  i <- 1
  flag.L <- rep(FALSE, nrow(TS))
  keep.track <- 1
  depth.track <- 0
  depth <- 0
  while (sum(flag.L) != nrow(TS)) {
    if (!flag.L[i]) {
      if (TS[i, 2] == 0) {
        flag.L[i] <- TRUE
        n.temp <- length(TreePrint)
        tempp <- strsplit(TreePrint[n.temp], ") ")[[1]]
        temp.L <- paste(tempp[1], ")#", tempp[2], sep = "")
        temp.L <- paste(temp.L, " -> ", "(", "leaf", leaf[i], " = ", ifelse(ppTree$type != "regression", ppTree$Levels[TS[i, 3]], TS[i, 3]), ")", sep = "")
        TreePrint <- TreePrint[-n.temp]
        id.l <- length(keep.track) - 1
        i <- keep.track[id.l]
        depth <- depth - 1
      } else if (!flag.L[TS[i, 2]]) {
        depth <- depth + 1
        emptyspace <- pastemake(depth, "   ")
        temp.L <- paste(emptyspace, "node", TS[i, 2], ")  proj",
          TS[i, 4], "*X < ", round(CutValue[TS[i, 4]], 2),
          sep = ""
        )
        i <- TS[TS[i, 2], 1]
      } else {
        depth <- depth + 1
        emptyspace <- pastemake(depth, "   ")
        temp.L <- paste(emptyspace, "node", TS[i, 3], ")  proj",
          TS[i, 4], "*X >= ", round(CutValue[TS[i, 4]], 2),
          sep = ""
        )
        flag.L[i] <- TRUE
        i <- TS[TS[i, 3], 1]
      }
      keep.track <- c(keep.track, i)
      depth.track <- c(depth.track, depth)
      TreePrint <- c(TreePrint, temp.L)
    } else {
      id.l <- id.l - 1
      i <- keep.track[id.l]
      depth <- depth.track[id.l]
    }
  }
  colnames(Alpha) <- ppTree$data$varName
  rownames(Alpha) <- paste("proj", seq_len(nrow(Alpha)), sep = "")
  # colnames(CutValue)<-paste("Rule",1:ncol(CutValue),sep="")
  names(CutValue) <- paste("CutValue", seq_along(CutValue), sep = "")
  TreePrint.output <-
    paste(
      "=============================================================",
      "\nOblique", ifelse(ppTree$type == "regression", "Regression", "Classification"), "Tree structure",
      "\n=============================================================\n"
    )
  for (i in seq_along(TreePrint)) {
    TreePrint.output <- paste(TreePrint.output, TreePrint[i], sep = "\n")
  }
  TreePrint.output <- paste(TreePrint.output, "\n", sep = "")
  # colnames(Alpha)<-paste(1:ncol(Alpha),":\"",colnames(Alpha),"\"",sep="")
  if (verbose) {
    cat(TreePrint.output)
    if (projection) {
      cat(
        "\nProjection coefficient in each node",
        "\n-------------------------------------------------------------\n"
      )
      print(round(Alpha, 4))
    }
    if (cutvalue) {
      cat(
        "\nCutoff values of each node",
        "\n-------------------------------------------------------------\n"
      )
      print(round(CutValue, 4))
    }
  }
  return(invisible(TreePrint))
}
