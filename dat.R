#regression
#https://www.kaggle.com/fedesoriano/body-fat-prediction-dataset
i.data==5
if(i.data==5){
  body_fat = data.frame(read.csv("D:/BaiduNetdiskWorkspace/DataSets/Regression1/bodyfat.csv"))
  X0 = data[,2:15]
  y0 = data[,1]
  p = ncol(X0)
  N = length(y0) #N=252,p=14
  Xy=data.frame(cbind(data[,2:15],data[,1]))
}
save(body_fat,file = "D:/R/Ryuyan/randomForest/Rpackage/ODRF/ODRF/data/body_fat.rda")


#classification
#https://www.kaggle.com/yasserh/breast-cancer-dataset
if(i.data==7){
  breast_cancer=read.csv("D:/BaiduNetdiskWorkspace/DataSets/Classification1/data7.csv")
  y0=as.integer(as.factor(yx[,2]))
  y0[y0==2]=0
  X0=as.matrix(yx[,3:32])
  n=length(y0)
  p=ncol(X0) #n=569,p=30
}
save(breast_cancer,file = "D:/R/Ryuyan/randomForest/Rpackage/ODRF/ODRF/data/breast_cancer.rda")

i.data=6
wd='D:/BaiduNetdiskWorkspace/'
wd0=getwd()
#http://archive.ics.uci.edu/ml/datasets/hill-valley
yx1=read.table(paste0(wd,"/DataSets/Classification1/","data3/Hill_Valley_without_noise_Testing.txt"),header=TRUE,sep=",")
yx2=read.table(paste0(wd,"/DataSets/Classification1/","data3/Hill_Valley_without_noise_Training.txt"),header=TRUE,sep=",")
Hill_Valley=rbind(yx1,yx2)
yx=Hill_Valley
n=dim(yx)[1]
p=dim(yx)[2]-1 #n=1212,p=100
yx[,p+1]=as.integer(as.factor(yx[,p+1]))
X0=as.matrix(yx[,1:p])
y0=rep(0,n)
y0[yx[,p+1]==2]=1
save(Hill_Valley,file = "D:/R/Ryuyan/randomForest/Rpackage/ODRF/ODRF/data/Hill_Valley.rda")

#multiclassification
#https://archive.ics.uci.edu/ml/datasets/seeds
i.data=16
wd='D:/BaiduNetdiskWorkspace/DataSets/ccfs-master/Datasets/classification1'
wd0=getwd()
setwd(wd)

files=dir(pattern=".csv")
files=files[order(sapply(files,function(file){as.numeric(strsplit(file,"-")[[1]][1])}))] 

XY=data.frame(read.csv(files[i.data]))
colnames(XY)=c("area","perimeter","compactness","length_of_kernel","width_of_kernel","asymmetry_coefficient","length_of_kernel_groove","varieties_of_wheat")
seeds=XY

X=XY[,-ncol(XY)]
y=XY[,ncol(XY)]

I = which((rowSums(is.na(X))>0)|is.na(y))#,which(rowSums(X==0)==p))
if(length(I)>0){
  X=X[-I,]
  y=y[-I]
}

#options (warn = 1)
setwd(wd0)


save(seeds,file = "D:/R/Ryuyan/randomForest/Rpackage/ODRF/ODRF/data/seeds.rda")

