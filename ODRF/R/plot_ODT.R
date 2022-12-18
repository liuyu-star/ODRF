#' oblique decision tree plot
#' 
#' Draw oblique decision tree with tree structure. It is modified from a function in \code{PPtreeViz} library.
#' 
#' @param ppTree an object of class \code{\link{ODT}}.
#' @param font.size font size of plot
#' @param width.size size of eclipse in each node.
#' @param xadj The size of the left and right movement.
#' @param main main title
#' @param sub sub title
#' @param ... arguments to be passed to methods.
#' @references Lee, EK(2017) PPtreeViz: An R Package for Visualizing Projection Pursuit Classification Trees, Journal of Statistical Software <doi:10.18637/jss.v083.i08>
#' @keywords tree
#' 
#' @seealso \code{\link{ODT}} \code{\link{plot_ODT_depth}}
#' 
#' @examples
#' data(iris)
#' tree <- ODT(Species~., data = iris,type='i-classification')
#' plot(tree)
#' 
#' @aliases plot.ODT
#' @rdname plot.ODT
#' @method plot ODT
#' @export
plot.ODT<-function(ppTree,font.size=17,width.size=1,xadj=0,main=paste0("Oblique ",
                   ifelse(ppTree$type=="regression","Regression","Classification")," Tree"),sub=NULL,...){
  numNode=length(ppTree$structure$nodeCutValue)
  cutNode=which(ppTree$structure$nodeCutValue!=0)
  
  TS=matrix(0,numNode,5)
  TS[,1]=seq(numNode)
  TS[,2]=ppTree[["structure"]][["childNode"]]
  if(ppTree$type!="regression"){
    TS[setdiff(seq(numNode),cutNode),3]=max.col(ppTree$structure$nodeNumLabel)[setdiff(seq(numNode),cutNode)]
  }else{
    TS[setdiff(seq(numNode),cutNode),3]=round(ppTree$structure$nodeNumLabel[,1][setdiff(seq(numNode),cutNode)],3)
  }
  TS[cutNode,3]=TS[cutNode,2]+1
  TS[cutNode,4]=seq(length(cutNode))
  TS[cutNode,5]=ppTree[["structure"]][["nodeCutIndex"]][cutNode]
  colnames(TS)=c("node","left_node","right_node/leaf_label","cut_node","cut_node_index")
  CutValue<-ppTree$structure$nodeCutValue[cutNode]
  
   plotPPtree<-function(ppTree,node.id,xlim,ylim){

      if(TS[node.id,2]==0){
         x<-xlim[1]+0.5  
         y<-ylim[2]-1     
         Final.Node.V<-viewport(x=unit(x,"native"),
                                y=unit(y,"native"),
                                width=unit(1,"native"), 
                                height=unit(1,"native")-unit(2,"lines"),
                                just=c("center","top"))
         pushViewport(Final.Node.V)
         node.terminal.PPtree(ppTree,node.id) 
         upViewport()
         return(NULL)
      }
      nl<-n.final(TS,node.id,"left")
      nr<-n.final(TS,node.id,"right")
      x0<-xlim[1]+nl 
      y0<-max(ylim)-1 
      lf<-ifelse(TS[TS[node.id,2],2]==0,0.5,
                    n.final(TS,TS[node.id,2],"right"))  
      rf<-ifelse(TS[TS[node.id,3],2]==0,0.5,
                    n.final(TS,TS[node.id,3],"left")) 
      x1l<-x0-lf
      x1r<-x0+rf
      y1<-y0-1
      grid.lines(x=unit(c(x0,x1l),"native"),y=unit(c(y0,y1),"native"))
      grid.lines(x=unit(c(x0,x1r),"native"),y=unit(c(y0,y1),"native"))
      node.V<-viewport(x=unit(x0,"native"),
                       y=unit(y0,"native"),
                       width=unit(1,"native"),
                       height=unit(1,"native")-unit(1,"lines"))
      pushViewport(node.V)
      node.inner.PPtree(ppTree,node.id)
      upViewport()
      #ylpos<-y0-0.6
      #yrpos<-y0-0.45
      ylpos<-y0-0.6
      yrpos<-y0-0.4
      xlpos<-x0-(x0-x1l)*0.6 
      xrpos<-x0-(x0-x1r)*0.4 
      #xrpos<-x0-(x0-x1r)*0.6
      LeftEdge.V<-viewport(x=unit(xlpos,"native"),
                           y=unit(ylpos,"native"),
                           width=unit(xlpos-xrpos,"native"),
                           height=unit(1,"lines")*1.2)
      pushViewport(LeftEdge.V)
      edge.lable.PPtree(TS,node.id, left = TRUE)
      upViewport()
      RightEdge.V<-viewport(x=unit(xrpos,"native"),
                            y=unit(yrpos,"native"),
                            width=unit(xlpos-xrpos,"native"),
                            height= unit(1,"lines"))
      pushViewport(RightEdge.V) 
      edge.lable.PPtree(TS,node.id,left=FALSE)
      upViewport()     
      plotPPtree(ppTree,TS[node.id,2],c(xlim[1],x0),c(1,y1+1))
      plotPPtree(ppTree,TS[node.id,3],c(x0,xlim[2]),c(1,y1+1))
   }    

   n.final<-function(TS,node.id,direction){  
      
      n.leaf<-0
      if(direction=="left"){
         keep.id<-TS[node.id,2]
         i<-1
         while(i<=length(keep.id)){
            if(TS[keep.id[i],2]==0){
               n.leaf<-n.leaf+1
               i<-i+1
            } else{  
               keep.id<-c(keep.id,TS[keep.id[i],2:3])
               i<-i+1
            }                               
         }  
      } else if(direction=="right"){
         keep.id<-TS[node.id,3]
         i<-1
         while(i<=length(keep.id)){
            if(TS[keep.id[i],2]==0){
               n.leaf<-n.leaf+1
               i<-i+1
            } else{  
               keep.id<-c(keep.id,TS[keep.id[i],2:3])
               i<-i+1
            }                               
         }  
      }
      return(n.leaf)
   }

   edge.lable.PPtree<-function(TS,node.id,left=TRUE){   
      
      if(left){
         text.t<-paste("< ",round(CutValue[TS[node.id,4]],2),sep="")
         grid.rect(gp=gpar(fill="white",lty=1,col="grey95"), 
                   width=unit(width.size,"strwidth", text.t)*1.2)
         grid.text(text.t,just="center",gp=gpar(fontsize=font.size))
      } else{
         text.t<-paste(">= ",round(CutValue[TS[node.id,4]],2),sep="")
         grid.rect(gp=gpar(fill="white",lty=1,col="grey95"), 
                  width=unit(width.size,"strwidth",text.t)*1.2)
         grid.text(text.t,just="center",gp=gpar(fontsize=font.size))
      } 
   }

   node.inner.PPtree<-function(ppTree,node.id){   
      
      PS<-print(ppTree,verbose=FALSE)
      label1<-rep(NA,length(PS))
      label2<-rep(NA,length(PS))
      ID<-rep(NA,length(PS))
      final.group<-rep(NA,length(PS))
      temp<-strsplit(PS,"->")
      for(i in 1:length(temp)){
         t<-strsplit(temp[[i]][1],")" )
         ID[i]<-as.numeric(substring(t[[1]][1],nchar(t[[1]][1])))
         tt<-strsplit(t[[1]][2]," ")[[1]]
         tt<-tt[tt!="" & tt!="#"]
         label1[i]<-tt[1]
         if(tt[1]!="root")
            label2[i]<-paste(tt[2],tt[3])
         if(length(temp[[i]])==2)
             final.group[i]<-temp[[i]][2]
      }
      label.t<-paste("proj",TS[node.id,4]," * X",sep="")
      Inner.Node.V<-viewport(x=unit(0.5,"npc"),
                             y=unit(0.5,"npc"),
                             width=unit(width.size*1.5,"strwidth",label.t), 
                             height=unit(width.size*2,"lines"))
      pushViewport(Inner.Node.V)
      xell<-c(seq(0,0.2,by=0.01),seq(0.2,0.8,by=0.05),seq(0.8,1,by=0.01))
	    yell<-sqrt(xell*(1-xell))	
      grid.polygon(x=unit(c(xell,rev(xell)),"npc"),
                   y=unit(c(yell,-yell)+0.5,"npc"),
                   gp=gpar(fill="white"))
      grid.text(label.t,y=0.3,gp=gpar(fontsize=font.size),
                just=c("center","bottom"))
      Inner.Node.Id.V<-viewport(x=unit(0.5,"npc"), 
                                y=unit(1,"npc"),
	                              width=max(unit(1,"lines"), 
                                          unit(1.2,"strwidth",
                                               as.character(node.id))),
	                              height=max(unit(1,"lines"), 
                                           unit(1.2,"strheight",
                                                as.character(node.id))),
                                just=c("center","center"),
                                gp=gpar(fontsize=font.size))
      pushViewport(Inner.Node.Id.V)
      grid.rect(gp=gpar(fill="white",lty="solid",fontsize=font.size))
      grid.text(node.id,gp=gpar(fontsize=font.size))
      popViewport()
      upViewport()
   }    

   node.terminal.PPtree<- function(ppTree,node.id){
      
      #gName<-names(table(PPtreeobj$origclass)) print xscale
      gN<-ifelse(ppTree$type!="regression",ppTree$Levels[TS[node.id,3]],TS[node.id,3])
      #gN<-ppTree$Levels[TS[node.id,3]]
      temp<-strsplit(as.character(gN),split="")[[1]]
      gN.width<-length(temp)     
      set.unit<-(sum(tolower(temp)!=temp)*0.65+
                 sum(tolower(temp)==temp)*0.5)/gN.width
      Terminal.Node.V<-viewport(x=unit(0.5,"npc"),   
                                y=unit(0.8,"npc"),   
                                width=unit(1,"lines")*1.5,
                                height=unit(set.unit,"lines")*(gN.width+2),
  	                            just=c("center","top"))
      pushViewport(Terminal.Node.V )	
      grid.rect(gp=gpar(fill="lightgray"))
      grid.text(y=0.05,gN,gp=gpar(fontsize=font.size),rot=90,just="left")
      Terminal.Node.Id.V<-viewport(x=unit(0.5,"npc"), 
                                   y=unit(1,"npc"),
                                   width=max(unit(1,"lines"), 
                                             unit(1.2,"strwidth",
                                                  as.character(node.id))),
                                   height=max(unit(1,"lines"), 
                                              unit(1.2,"strheight",
                                                   as.character(node.id))),
                                   just=c("center","center"),
                                   gp=gpar(fontsize=font.size))
      pushViewport(Terminal.Node.Id.V)
      grid.rect(gp=gpar(fill="lightgray",lty="solid",fontsize=font.size))
      grid.text(node.id,gp=gpar(fontsize=font.size))
      popViewport()
      upViewport()    
   }

   ########################################
   #length(ppTree$Levels)
   nx<-numNode-length(cutNode)
   ny<-max(ppTree$structure$nodeDepth)+1
   #ny<-calc.depth(TS)
   tnex<-1
   node.id<-1
   grid.newpage()
   PPtree.Main.V<-viewport(layout=grid.layout(3,3, 
                                      heights=unit(c(3,1,1),
                                                   c("lines","null","lines")),
    			                            widths=unit(c(1,1,1),
                                                  c("lines","null","lines"))))       
   pushViewport(PPtree.Main.V)
   PPtree.title.V<-viewport(layout.pos.col=2,layout.pos.row=1)
   pushViewport(PPtree.title.V)
   grid.text(y=unit(1,"lines"),
             paste("\n",main,sep=""),
             just="center",gp=gpar(fontsize=font.size))
   upViewport()
   PPtree.Tree.V<-viewport(layout.pos.col=2,layout.pos.row=2,xscale=c(0,nx-xadj),yscale=c(0,ny+1)) 
    			                 #xscale=c(0,nx),yscale=c(0,ny+1))
   pushViewport(PPtree.Tree.V)
   
   
   plotPPtree(ppTree,1,c(0,nx-xadj),ylim=c(1,ny+1)) 
   grid.text(y=unit(1,"lines"),
             sub,
             just="center",gp=gpar(fontsize=font.size*0.7))
}
 