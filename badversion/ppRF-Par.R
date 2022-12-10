
cores=35;Rep=ntrees
cores <-min(c(detectCores(logical = FALSE),cores))
cl <- makeCluster(cores)
library("ppRF1")

chunks <- clusterSplit(cl,seq(Rep))
registerDoParallel(cl, cores=cores)


ico = 1
#N.simu=100
Err = foreach(co=seq(length(chunks)), .combine=list,.multicombine = TRUE,
              .packages = "ppRF1",.noexport ="ppForest") %dopar%
  {

    E=vector("list",length(chunks[[co]]))
    
    parfun=function(ico){
      set.seed(ico)
      
      #ppForestT=vector("list", 6);
      #ppForestT[6]=c(1);
      #names(ppForestT)=c("nodesparseM","nodeCutVar", "nodeCutValue","childnode" ,"nodelabel","oobErr")
      ppForestT=list()
      
      TDindx=sample.int(N-numOOB);
      ppForestT=c(ppForestT,carPPtree(X[TDindx,],y[TDindx],FUN,paramList,catMap, minparent,minleaf,
                                nvartosample,method,weights,Levels));
      
      oobErr=1
      if (numOOB>0){
        NTD = setdiff(1:N,TDindx);
        tree_output = carPPtree_Prediction_Cpp(X[NTD,],ppForestT[1:5][[1]],ppForestT[1:5][[2]],
                                               ppForestT[1:5][[3]],ppForestT[1:5][[4]],ppForestT[1:5][[5]]);
        
        if(tolower(method)%in%c('c','g')){
          oobErr=mean(tree_output!=y[NTD]);
        }else{
          oobErr=mean((tree_output-y[NTD])^2);
        }
      }
      ppForestT=c(ppForestT,oobErr=oobErr)

      return(ppForestT)
    }
    
    for(ico in chunks[[co]]) {
      
      ppForestT = try(parfun(ico), silent=TRUE)
      if ('try-error' %in% class(ppForestT)){
        next
      }
      
      E[[ico-(chunks[[co]][1]-1)]] = ppForestT
    }
    print(ico)
    # return local results
    E
  }

ne=ncol(Err)
err=which(rowSums(Err==1)==ne)
if((length(err)>0)&(length(err)<ne)){
  Err=Err[-err,]
}




