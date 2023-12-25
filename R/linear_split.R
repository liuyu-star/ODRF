#' @keywords internal
#' @noRd
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet assess.glmnet
#library(glmnet)
linear_split=function(X, y, Xsplit, minleaf=10, lambda=0, numLabels, glmnetParList=list(lambda = NULL)){
  y=c(y)
  n=length(y)
  p=NCOL(X)
  glmnetParList$x <- X
  glmnetParList$y <- y
  if(length(glmnetParList$lambda)==1){
    initfit <- do.call(glmnet, glmnetParList)
    yhat=predict(initfit, X)
  }else{
    initfit <- do.call(cv.glmnet, glmnetParList)
    glmnetParList$lambda=initfit$lambda.min
    yhat=predict(initfit, X, s = "lambda.min")
  }
  #glmnetParList$X <- NULL
  #glmnetParList$y <- NULL
  rss0=sum((y-yhat)^2) #rss0+2*p
  TFclass=ifelse(is.null(glmnetParList$family),FALSE,glmnetParList$family%in%c("binomial","multinomial"))
  t=ifelse(lambda == "log",log(n),lambda)
  if(TFclass)rss0=rss0+2*p else rss0=rss0*(n/(n-t))^2
  rss00=rss0

  bcvar=-1
  bcval=0
  fitL0=fitR0=NULL
  #n=length(y)
  #p=ncol(X)
  ps=ncol(Xsplit)
  bestval=rep(0,ps)
  sps=ceiling(seq(minleaf,n-minleaf,length.out = min(n-2*minleaf+1,max(ceiling(n/10),100))))
  for(sv in seq(ps)){
   #x=Xsplit[,sv]
   xs=sort(Xsplit[,sv],index.return = TRUE)
   idx=xs$ix
   xs=xs$x
   #ys=y[idx]

   xs0=-Inf
   minrss=Inf
    for(sp in sps){
      if(xs[sp]==xs0){
        next
      }
      xs0=xs[sp]
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
      rssL=c(assess.glmnet(fitL, newx = glmnetParList$x, newy = glmnetParList$y)[[1]])
      #yhat=predict(fitL,glmnetParList$x)
      #rss=sum((glmnetParList$y-yhat)^2)

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
      rssR=c(assess.glmnet(fitR, newx = glmnetParList$x, newy = glmnetParList$y)[[1]])
      #yhat=predict(fitR,glmnetParList$x)
      #rss=rss+sum((glmnetParList$y-yhat)^2)

      if (lambda=="log"){
        tl=log(sp);
        tr=log(n-sp);
      }else{
        tl=lambda;
        tr=lambda;
      }
      if(TFclass){
        rss=rssL+2*p+ rssR+2*p
      }else{
        rss=rssL*(sp/(sp-tl))^2+ rssR*((n-sp)/((n-sp)-tr))^2
      }

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
