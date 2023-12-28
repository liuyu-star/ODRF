#' @keywords internal
#' @noRd
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet assess.glmnet
#library(glmnet)
linear_split=function(X, y, Xsplit, minleaf=10, lambda=0, numLabels, glmnetParList=list(lambda = NULL)){
  y=c(y)
  n=length(y)
  p=NCOL(X)
  TFclass=ifelse(is.null(glmnetParList$family),FALSE,glmnetParList$family%in%c("binomial","multinomial"))
  bcvar=-1
  bcval=0
  fitL0=fitR0=NULL
  ps=ncol(Xsplit)
  bestval=rep(0,ps)

  glmnetParList$x <- X
  glmnetParList$y <- y
  if(TFclass){
    ty=table(glmnetParList$y)
    if(min(ty)<2|length(ty)!=numLabels){
      return(list(BestCutVar=bcvar, BestCutVal=bcval, BestIndex=bestval, fitL=fitL0, fitR=fitR0))
      #next
      #jdx=which(glmnetParList$y==names(ty)[which.min(ty)])
      #glmnetParList$y=c(glmnetParList$y,glmnetParList$y[jdx])
      #glmnetParList$x=rbind(glmnetParList$x,glmnetParList$x[jdx,])
    }
  }
  if(length(glmnetParList$lambda)==1){
    initfit <- do.call(glmnet, glmnetParList)
    #yhat=predict(initfit, X)
  }else{
    initfit=try(do.call(cv.glmnet, glmnetParList),silent = T)
    if(is(initfit, "try-error")){
      glmnetParList$lambda=0.001
    }else{
      glmnetParList$lambda=initfit$lambda.min
    }
    initfit <- do.call(glmnet, glmnetParList)
    #yhat=predict(initfit, X, s = "lambda.min")
  }
  #glmnetParList$X <- NULL
  #glmnetParList$y <- NULL
  #rss0=sum((y-yhat)^2) #rss0+2*p
  #t=ifelse(lambda == "log",log(n),lambda)
  #if(TFclass)rss0=rss0+2*p else rss0=rss0*(n/(n-t))^2
  rss0 = (1-initfit$dev.ratio)*initfit$nulldev + 2*p
  rss00=rss0

  #n=length(y)
  #p=ncol(X)
  #sps=ceiling(seq(minleaf,n-minleaf,length.out = min(n-2*minleaf+1,100)))#max(ceiling(n/10),100))))
  #ns=length(sps)
  for(sv in seq(ps)){
   #x=Xsplit[,sv]
   Xs=Xsplit[,sv]
   xs=sort(Xs,index.return = TRUE)
   idx=xs$ix
   xs=xs$x
   #ys=y[idx]

   tx=unique(Xs)
   nx=length(tx)
   itx=rep(0,nx)
   for (i in seq(nx)) {
     itx[i]=min(which(Xs==tx[i]))
   }
   #itx=unique(quantile(unique(Xs), (1:100)/100, type=1))
   #if(length(itx)<100)sps=itx else sps=union(sps,itx)
   itx=itx[seq.int(1,nx,length.out = min(nx,50))]
   sps=itx[itx>=minleaf&itx<=n-minleaf]

   #xs0=-Inf
   minrss=Inf
    for(sp in sps){
      #if(xs[sp]==xs0){
       # next
      #}
      #xs0=xs[sp]
      idx0=idx[seq(sp)]
      glmnetParList$x <- X[idx0,]
      glmnetParList$y <- y[idx0]
      if(TFclass){
        ty=table(glmnetParList$y)
        if(min(ty)<2|length(ty)!=numLabels){
          next
          #jdx=which(glmnetParList$y==names(ty)[which.min(ty)])
          #glmnetParList$y=c(glmnetParList$y,glmnetParList$y[jdx])
          #glmnetParList$x=rbind(glmnetParList$x,glmnetParList$x[jdx,])
        }
      }
      fitL=do.call(glmnet, glmnetParList)
      #rssL=c(assess.glmnet(fitL, newx = glmnetParList$x, newy = glmnetParList$y)[[1]])
      #yhat=predict(fitL,glmnetParList$x)
      #rss=sum((glmnetParList$y-yhat)^2)
      rssL = (1-fitL$dev.ratio)*fitL$nulldev + 2*p

      glmnetParList$x <- X[-idx0,]
      glmnetParList$y <- y[-idx0]
      if(TFclass){
        ty=table(glmnetParList$y)
        if(min(ty)<2|length(ty)!=numLabels){
          next
          #jdx=which(glmnetParList$y==names(ty)[which.min(ty)])
          #glmnetParList$y=c(glmnetParList$y,glmnetParList$y[jdx])
          #glmnetParList$x=rbind(glmnetParList$x,glmnetParList$x[jdx,])
        }
      }
      fitR=do.call(glmnet, glmnetParList)
      #rssR=c(assess.glmnet(fitR, newx = glmnetParList$x, newy = glmnetParList$y)[[1]])
      #yhat=predict(fitR,glmnetParList$x)
      #rss=rss+sum((glmnetParList$y-yhat)^2)
      rssR = (1-fitR$dev.ratio)*fitR$nulldev

      #if (lambda=="log"){
      #  tl=log(sp);
      #  tr=log(n-sp);
      #}else{
      #  tl=lambda;
      #  tr=lambda;
      #}
      #if(TFclass){
      #  rss=rssL+2*p+ rssR+2*p
      #}else{
      #  rss=rssL*(sp/(sp-tl))^2+ rssR*((n-sp)/((n-sp)-tr))^2
      #}
      rss=rssL+2*p+ rssR+2*p

      if(rss<rss0){
        rss0=rss;
        bcvar=sv
        bcval=mean(xs[c(0,1)+sp])
        fitL0=fitL
        fitR0=fitR
      }
      if(rss<minrss){
        minrss=rss;
      }
    }

    bestval[sv]=rss00-minrss
  }

  return(list(BestCutVar=bcvar, BestCutVal=bcval, BestIndex=bestval, fitL=fitL0, fitR=fitR0))
}


#Example
#set.seed(10)
#cutpoint=50
#x=matrix(rnorm(100*10),100,10)
#age=sample(seq(20,80),100,replace = TRUE)
#height=sample(seq(50,200),100,replace = TRUE)
#weight=sample(seq(5,150),100,replace = TRUE)
#Xsplit=cbind(age=age,height=height,weight=weight)
#mu=rep(0,100)
#mu[age<=cutpoint]=x[age<=cutpoint,1]+x[age<=cutpoint,2]
#mu[age>cutpoint]=x[age>cutpoint,1]+x[age>cutpoint,3]
#y=mu+rnorm(100)

#treefit=linear_split(X, Xsplit, y, glmnetParList=list(alpha=0.01))
