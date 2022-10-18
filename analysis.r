library(openxlsx)
library(glmnet)
library(dplyr)
library(pROC)

# load data
rm(list = ls()) 
ori.data = read.xlsx("datafile.xlsx", sheet = "Sheet1")
data=select(ori.data,c("China.MetS","ATP.MetS","IDF.Mets","Age","TC","SBP","DBP","CKDEPI",
                       "BMI","LAP","BRI","CVAI","AVI","BAI","TYG","VAI","Sex"))

# randomize order
set.seed(0)
index = sample(nrow(data))
data=data[index,]

# split by sex
data.male=data[data$Sex==1,]
data.female=data[data$Sex==2,]
data.by.sex=list(data.male,data.female,data)
genders=c('male','female','all')

# generate fixed fold ids
sample.counts=lapply(data.by.sex,nrow)
n.folds=10
foldids=vector(mode='list',length = 3)
for(i in 1:3){
  foldids[[i]]=rep_len(c(1:n.folds),sample.counts[[i]])
}

# assign variables to different groups
target.list=colnames(data)[1:3]
additional.variables=colnames(data)[4:8]
variable.list=colnames(data)[9:16]

# cross validation of base logistic model
n.fold.logistic=function(data,foldid,target,variables){
  val.aucs=c()
  n.fold=max(foldid)
  for(fold in 1:n.fold){
    split.val=(foldid==fold)
    data.train=data[!split.val,]
    data.val=data[split.val,]
    form=as.formula(paste(target,paste(" ~ ",paste(variables, collapse= "+"))))
    f=glm(form,data=data.train,family='binomial')
    pred.f=predict(f,newdata=data.val,type = 'response')
    val.result=as.numeric(pred.f)
    val.aucs[fold]=auc(data.val[,target],val.result)
  }
  list(lambda=NA,mean.val.aucs=mean(val.aucs))
}

# cross validation of elasticnet logistic model with different alpha
func.elasticnet=function(alpha){
  a=alpha
  n.fold.elasticnet=function(data,foldid,target,variables){
    val.aucs=c()
    n.fold=max(foldid)
    f=cv.glmnet(as.matrix(select(data,all_of(variables))),
                data[,target],foldid = foldid,
                family='binomial',alpha = a, type.measure="auc")
    list(lambda=f$lambda.min,mean.val.aucs=f$cvm[f$index[1]])
  }
}

# model building
models.eval=list(n.fold.logistic)
models=c('Base','Ridge (a=0)')
alphas=seq(0,1,0.1)
for(i in 1:length(alphas)){
  models.eval=append(models.eval,func.elasticnet(alphas[i]))
}
for(i in 1:(length(alphas)-2)){
  models=append(models,paste0('Elasticnet a=',alphas[i+1]))
}
models=append(models,'Lasso (a=1)')

# result tables
mean.aucs=data.frame(matrix(NA,nrow=8,ncol=length(models)))
colnames(mean.aucs)=models
rownames(mean.aucs)=variable.list
list1=list(china=mean.aucs,ATP=mean.aucs,IDF=mean.aucs)
list2=list(male=list1,female=list1,all=list1)

# model runing
count=1
for(i in 1:3){
  for(j in 1:8){
    for(sex in 1:3){
      for(model in 1:length(models.eval)){
        print(count)
        count=count+1
        target=target.list[i]
        variables=c(variable.list[j],additional.variables)
        form=as.formula(paste(target,paste(" ~ ",paste(variables, collapse= "+"))))
        #print(form)
        #print(models[model])
        list2[[sex]][[i]][j,model]=models.eval[[model]](data.by.sex[[sex]],foldids[[sex]],target,variables)$mean.val.aucs
      }
    }
  }
}

# result output
for(i in 1:3){
  for(sex in 1:3){
    mean.aucs=list2[[sex]][[i]]
    mean.aucs$best.model=models[apply(mean.aucs[,1:12],1,which.max)]
    mean.aucs$best.result=apply(mean.aucs[,1:12],1,max)
    write.xlsx(mean.aucs,paste0(genders[sex],'_',target.list[i],'.xlsx'),rowNames=T) 
  }
}
                            
