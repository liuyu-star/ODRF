#' \code{ODT} as \code{party}
#'
#' To make \code{ODT} object to objects of class \code{party}.
#'
#' @param obj An object of class \code{\link{ODT}}.
#' @param data Training data of class \code{data.frame} is used to convert the object of class \code{ODT},
#' and it must be the training data \code{data} in \code{\link{ODT}}.
#' @param ... Arguments to be passed to methods
#'
#' @return An objects of class \code{party}.
#'
#' @references Lee, EK(2017) PPtreeViz: An R Package for Visualizing Projection Pursuit Classification Trees, Journal of Statistical Software.
#'
#' @seealso \code{\link{ODT}} \code{\link{party}}
#'
#' @examples
#' data(iris)
#' tree <- ODT(Species ~ ., data = iris)
#' tree
#' plot(tree)
#' party.tree <- as.party(tree, data = iris)
#' party.tree
#' plot(party.tree)
#'
#' @keywords tree
#' @import partykit
#' @aliases as.party.ODT
#' @rdname as.party.ODT
#' @method as.party ODT
#' @export
as.party.ODT <- function(obj, data, ...) {
  # if(is.null(data)){
  #  data <- data.frame(y=eval(formula[[2]]),eval(formula[[3]]))
  # }
  ppTree <- obj
  rm(obj)

  if (prod(dim(data)) != ppTree$data$n * (ppTree$data$p+1)) {
    stop("'data' must be the training data 'data' in class ODT.")
  }

  vars <- all.vars(ppTree$terms)
  y <- data[, setdiff(colnames(data), vars[-1])]
  X <- data[, vars[-1]]
  X <- as.matrix(X)

  numNode <- length(ppTree$structure$nodeCutValue)
  cutNode <- which(ppTree$structure$nodeCutValue != 0)

  TS <- matrix(0, numNode, 5)
  TS[, 1] <- seq(numNode)
  TS[, 2] <- ppTree[["structure"]][["childNode"]]
  if (ppTree$split != "mse") {
    TS[setdiff(seq(numNode), cutNode), 3] <- max.col(ppTree$structure$nodeNumLabel)[setdiff(seq(numNode), cutNode)]
  } else {
    TS[setdiff(seq(numNode), cutNode), 3] <- round(ppTree$structure$nodeNumLabel[, 1][setdiff(seq(numNode), cutNode)], 3)
  }
  TS[cutNode, 3] <- TS[cutNode, 2] + 1
  TS[cutNode, 4] <- seq_along(cutNode)
  TS[cutNode, 5] <- ppTree[["structure"]][["nodeCutIndex"]][cutNode]
  colnames(TS) <- c("node", "left_node", "right_node/leaf_label", "cut_node", "cut_node_index")
  TS <- data.frame(TS)
  CutValue <- ppTree$structure$nodeCutValue[cutNode]


  n <- nrow(TS)
  if (n == 1) {
    return(partykit::partynode(as.integer(1)))
  }
  is.leaf <- (TS$cut_node == 0)
  ncomplete <- rep(2, n)
  ncomplete[is.leaf] <- 0
  index <- cumsum(c(1, ncomplete + 1 * (!is.leaf)))
  primary <- numeric(n)
  primary[!is.leaf] <- index[c(!is.leaf, FALSE)]

  nodeRotaMat <- ppTree[["structure"]][["nodeRotaMat"]]
  Alpha <- matrix(0, length(cutNode), ppTree[["data"]][["p"]])
  for (cn in seq_along(cutNode)) {
    idx <- which(nodeRotaMat[, 2] == cutNode[cn])
    Alpha[cn, nodeRotaMat[idx, 1]] <- nodeRotaMat[idx, 3]
  }

  mf <- X %*% t(Alpha)
  rownames(mf) <- rownames(X)
  colnames(mf) <- paste0(paste0("proj", seq_len(ncol(mf)), "X"))
  mf <- data.frame(mf)

  fit <- as.data.frame(matrix(nrow = NROW(mf), ncol = 0))
  # fit <- as.data.frame(matrix(nrow = NROW(mf), ncol = 0))
  # if(ppTree$split!="regression"){
  #  pred=as.numeric(factor(predict_ppCART(ppTree,X),levels = ppTree$Levels))
  # }else{
  #  pred=predict_ppCART(ppTree,X)
  # }
  # fit[["(fitted)"]] <- apply(matrix(pred,ncol=1),1,function(x) which((TS[,3]==x)*is.leaf==1))
  fit[["(fitted)"]] <- predict(ppTree, X, leafnode = TRUE)
  fit[["(response)"]] <- y

  # fitted <- predict_ppCART(tree,X)
  pptree_kids <- function(i) {
    if (is.leaf[i]) {
      return(NULL)
    } else {
      return(c(TS[i, c(3, 2)]))
    }
  }

  pptree_split <- function(j) {
    if (j < 1) {
      return(NULL)
    }
    idj <- as.integer(TS$cut_node[j])
    ret <- partykit::partysplit(
      varid = idj,
      breaks = as.double(CutValue[idj]),
      right = FALSE,
      index = 2L:1L
    )
    ret
  }

  pptree_node <- function(i) {
    if (is.null(pptree_kids(i))) {
      return(partynode(as.integer(i)))
    }

    nd <- partykit::partynode(as.integer(i),
      split = pptree_split(i),
      kids = lapply(pptree_kids(i), pptree_node)
    )

    left <- partykit::nodeids(kids_node(nd)[[1L]], terminal = TRUE)
    right <- partykit::nodeids(kids_node(nd)[[2L]], terminal = TRUE)
    nd$split$prob <- c(0, 0)
    nl <- sum(fit[["(fitted)"]] %in% left)
    nr <- sum(fit[["(fitted)"]] %in% right)
    if (nl > nr) {
      nd$split$prob <- c(1, 0)
    } else {
      nd$split$prob <- c(0, 1)
    }
    nd$split$prob <- as.double(nd$split$prob)
    return(nd)
  }

  node <- pptree_node(1)
  # rval <- pp_party(node = node,
  rval <- partykit::party(
    node = node,
    data = mf, # if (data) mf else mf[0L, ],
    fitted = fit,
    terms = ppTree$terms,
    info = list(method = "ODT")
  )
  class(rval) <- c("constparty", class(rval))
  # class(me) <- append(class(rval),"constparty")
  # rval
  return(rval)
}
