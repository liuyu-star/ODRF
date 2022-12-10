#' Print PP.Tree.class result
#' 
#' Print the projection pursuit classification tree result
#' @title Print PP.Tree.class result
#' @param x PPtreeclass object
#' @param coef.print print projection coefficients in each node ifTRUE
#' @param cutoff.print print cutoff values in each node if TRUE
#' @param verbose print if TRUE, no output if FALSE
#' @param ... arguments to be passed to methods
#' @references Lee, EK(2017) 
#' PPtreeViz: An R Package for Visualizing Projection Pursuit Classification 
#' Trees, Journal of Statistical Software <doi:10.18637/jss.v083.i08>
#' 
#' @keywords tree
#' @aliases print.ODT
#' @rdname print.ODT
#' @method print ODT
#' @export
#' 
#' @examples
#' data(iris)
#' Tree.result <- PPTreeclass(Species~.,data = iris,"LDA")
#' Tree.result
#' print(Tree.result,coef.print=TRUE,cutoff.print=TRUE)tree
print.ODT<-function(tree,projection=FALSE,cutvalue=FALSE,verbose=TRUE,...){
  numNode=length(tree$structure$nodeCutValue)
  cutNode=which(tree$structure$nodeCutValue!=0)
  
  TS=matrix(0,numNode,5)
  TS[,1]=seq(numNode)
  TS[,2]=tree[["structure"]][["childNode"]]
  if(tree$method!="regression"){
    TS[setdiff(seq(numNode),cutNode),3]=max.col(tree$structure$nodeNumLabel)[setdiff(seq(numNode),cutNode)]
  }else{
    TS[setdiff(seq(numNode),cutNode),3]=round(tree$structure$nodeNumLabel[,1][setdiff(seq(numNode),cutNode)],3)
  }
  TS[cutNode,3]=TS[cutNode,2]+1
  TS[cutNode,4]=seq(length(cutNode))
  TS[cutNode,5]=tree[["structure"]][["nodeCutIndex"]][cutNode]
  colnames(TS)=c("node","left_node","right_node/leaf_label","cut_node","cut_node_index")
  leaf=rep(0,numNode)
  leaf[setdiff(seq(numNode),cutNode)]=seq(numNode-length(cutNode))
  
  #TS<-tree$Tree.Struct
  #Alpha<-tree$projbest.node
  nodeRotaMat<-tree[["structure"]][["nodeRotaMat"]]
  Alpha=matrix(0,length(cutNode),tree[["data"]][["p"]])
  for (cn in  1:length(cutNode)) {
    idx=which(nodeRotaMat[,2]==cutNode[cn])
    Alpha[cn,nodeRotaMat[idx,1]]=nodeRotaMat[idx,3]
  }
  
  CutValue<-tree$structure$nodeCutValue[cutNode]
  #CutValue<-tree$splitCutoff.node
  #gName<-tree$Levels
  #gName<-names(table(tree$origclass))
  pastemake<-function(k,arg,sep.arg=""){
    temp<-""
    for(i in 1:k)
      temp<-paste(temp,arg,sep=sep.arg)
    return(temp)
  }
  TreePrint<-"1) root"
  i<-1  
  flag.L<-rep(FALSE,nrow(TS))
  keep.track<-1
  depth.track<-0
  depth<-0
  while(sum(flag.L)!=nrow(TS)){
    if(!flag.L[i]){                    
      if(TS[i,2] == 0) {
        flag.L[i]<-TRUE
        n.temp<-length(TreePrint)
        tempp<-strsplit(TreePrint[n.temp],") ")[[1]]
        temp.L<-paste(tempp[1],")#",tempp[2],sep="")
        temp.L<- paste(temp.L," -> ","(","leaf",leaf[i]," = ",ifelse(tree$method!="regression",tree$Levels[TS[i,3]],TS[i,3]),")",sep="")
        TreePrint<-TreePrint[-n.temp]
        id.l<-length(keep.track)-1
        i<-keep.track[id.l]
        depth<-depth -1
      } else if(!flag.L[TS[i,2]]){
        depth<-depth+1
        emptyspace<-pastemake(depth,"   ")
        temp.L<-paste(emptyspace,"node",TS[i,2],")  proj",
                      TS[i,4],"*X < ",round(CutValue[TS[i,4]],2),sep="")
        i<-TS[TS[i,2],1]   
      } else{
        depth<-depth +1
        emptyspace<-pastemake(depth,"   ")          
        temp.L<- paste(emptyspace,"node",TS[i,3],")  proj",
                       TS[i,4],"*X >= ",round(CutValue[TS[i,4]],2),sep="")
        flag.L[i]<-TRUE
        i<-TS[TS[i,3],1]
      } 
      keep.track<-c(keep.track,i)
      depth.track<-c(depth.track,depth)
      TreePrint<-c(TreePrint,temp.L)
    } else{
      id.l<-id.l-1
      i<-keep.track[id.l]
      depth<-depth.track[id.l]
    }
  }
  colnames(Alpha)<-tree$data$varName
  rownames(Alpha)<-paste("proj",1:nrow(Alpha),sep="")  
  #colnames(CutValue)<-paste("Rule",1:ncol(CutValue),sep="")
  names(CutValue)<-paste("CutValue",1:length(CutValue),sep="")
  TreePrint.output<-
    paste("=============================================================",                          
          "\nOblique",ifelse(tree$method=="regression","Regression","Classification"),"Tree result",                           
          "\n=============================================================\n")
  for(i in 1:length(TreePrint))
    TreePrint.output<-paste(TreePrint.output,TreePrint[i],sep="\n")
  TreePrint.output<-paste(TreePrint.output,"\n",sep="")
  colnames(Alpha)<-paste(1:ncol(Alpha),":\"",colnames(Alpha),"\"",sep="")
  if(verbose){
    cat(TreePrint.output)
    if(projection){
      cat("\nProjection Coefficient in each node",
          "\n-------------------------------------------------------------\n")
      print(round(Alpha,4))
    }
    if(cutvalue){
      cat("\nCutoff values of each node",
          "\n-------------------------------------------------------------\n")
      print(round(CutValue,4))
    }    
  }    
  return(invisible(TreePrint)) 
}


   