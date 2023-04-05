#' @keywords internal
#' @noRd
#'
predictTree=function(structure, Xnew, split, Levels, ...){
  nodeNumLabel=structure$nodeNumLabel
  if (split %in% c("gini","entropy")) {
    if(all(structure$nodeCutValue == 0)){
      nodeNumLabel=matrix(nodeNumLabel,nrow = 1, ncol = length(Levels))
    }
    nodeLabel <- Levels[max.col(nodeNumLabel)]
    nodeLabel[which(rowSums(nodeNumLabel) == 0)] <- "0"
  } else {
    nodeLabel <- as.character(nodeNumLabel[, 1])
  }

  if (all(structure$nodeCutValue == 0)) {
    pred <- rep(nodeLabel, nrow(Xnew))
    node <- rep(0, nrow(Xnew))
  } else {
    predict_tree <- .Call("_ODRF_predict_ODT",
                  PACKAGE = "ODRF", Xnew, structure$nodeRotaMat,
                  structure$nodeCutValue, structure$childNode, nodeLabel)
    pred=predict_tree$prediction
    node <- as.integer(predict_tree$node)
  }

  if (!split %in% c("gini","entropy")) {
    pred <- as.numeric(pred)
  }

  return(list(prediction=pred,leafnode=node))
}
