library(plyr)
library(randomForest)
library(e1071)
library(neuralnet)
library(tidyverse)
library(pROC)
library(ape)
library(vegan)
library(scatterplot3d)
library(factoextra)
library(Metrics)
setwd("D:\\microbiomeML")#���ù���·��
#��������
data <- read.csv("S1_table.csv")#OTU��������
newdata<-read.csv("S2_table.csv")#��ˮƽ����

#���ݴ���
newdata1<-t(newdata)
new1<-as_tibble(newdata1)
newdata3<-new1[-1,-1] %>% 
  select(2:37) %>% 
  apply(., 2, as.numeric)
t3<-new1[-1,]$V1
#�������
newdata3.dist<-vegdist(newdata3)
#�����������
data3.pcoa<-pcoa(newdata3.dist,correction = "none")
#��ȡ��������
t2<-data3.pcoa[["vectors"]]
t1<-as.data.frame(cbind(t2,t3))
t1$t3<-as.factor(t1$t3)
eig1<-as.numeric(data3.pcoa$values[,1])
#��ɫ����
mycols<-ifelse(t1$t3=="N","magenta","blue")
mypchs<-ifelse(t1$t3=="N",5,17)
#scatterplot3d(t1$Axis.1,t1$Axis.2,t1$Axis.3,color  = mycols,xlab = "pc1",ylab = "pc2",zlab = "pc3",pch = mypchs,grid = TRUE)
#��ͼ
scatterplot3d(t1$Axis.1,t1$Axis.2,t1$Axis.3,color  = mycols,
              xlab = paste("PCoA1 (",format(100* eig1[1] / sum(eig1), digits=4), "%)", sep=""),
              ylab = paste("PCoA2 (",format(100* eig1[2] / sum(eig1), digits=4), "%)", sep=""),
              zlab = paste("PCoA3 (",format(100* eig1[3] / sum(eig1), digits=4), "%)", sep=""),
              pch = mypchs,grid = TRUE)
#��ȡ����
data.hc<-data %>%
  select(Samples,3:110)%>%
  column_to_rownames(var="Samples")
#�������
demo.dist<-dist(data.hc)
#�����ͼ
mymodel<-hclust(demo.dist)
res<-hcut(dist(data.hc),k=2,hc_func = "hclust",hc_method = "complete",hc_metric = "euclidean")
fviz_dend(res,rect = TRUE,type="rectangle")

#����Ԥ����
data<-data[,-1]
data$Malodour<-factor(ifelse(data$Malodour == "N", 1,2))

#����
CVgroup<-function(k,datasize,seed){
  cvlist<-data.frame()
  set.seed(seed)
  n<-rep(1:k,ceiling(datasize/k))
  temp<-sample(n,datasize)
  x<-1:k
  dataseq<-1:datasize
  cvlist<-lapply(x,function(x) dataseq[temp==x])
  return(cvlist)
}
k<-10
datasize<-90
cvlist<-CVgroup(k=k,datasize=datasize,seed=1234)
#�����洢���ݿ�
svm.p2<-data.frame()
rf.p2<-data.frame()
ne.p2<-data.frame()
#ģ�ͽ���
for (i in 1:10) {
  train<-data[-cvlist[[i]],]
  test<-data[cvlist[[i]],]
  #SVM
  svm.train<-svm(Malodour~.,train,probability=T)
  svm.pred <- predict(svm.train,test,probability=T)
  svm.pred<-attr(svm.pred,"probabilities")[,2]
  svm.p1<-data.frame(cbind(subset(test,select=Malodour),matrix(svm.pred)))
  svm.p2<-rbind(svm.p2,svm.p1)
  #randomforest
  rf.train<-randomForest(Malodour~.,train,ntree=500)
  rf.pred<-predict(rf.train,test,type="prob")[,2]
  rf.p1<-data.frame(cbind(subset(test,select=Malodour),matrix(rf.pred)))
  rf.p2<-rbind(rf.p2,rf.p1)
  #neuralnet
  n <- names(train)
  formula <- as.formula(paste("Malodour ~", paste(n[!n %in% "Malodour"], collapse = " + ")))
  neuralnet.train<-neuralnet(formula, train,hidden =c(500,500),rep=1,err.fct="ce",linear.output = FALSE)
  neuralnet.pred<-predict(neuralnet.train,test)[,2]
  ne.p1<-data.frame(cbind(subset(test,select=Malodour),matrix(neuralnet.pred)))
  ne.p2<-rbind(ne.p2,ne.p1)
}

#׼ȷ��
accuracy(ifelse(svm.p2[,2]<0.5,1,2),svm.p2[,1])
accuracy(ifelse(rf.p2[,2]<0.5,1,2),rf.p2[,1])
accuracy(ifelse(ne.p2[,2]<0.5,1,2),ne.p2[,1])
#AUC
roc1<-roc(svm.p2[,1],svm.p2[,2])
roc1$auc
roc2<-roc(rf.p2[,1],rf.p2[,2])
roc2$auc
roc3<-roc(ne.p2[,1],ne.p2[,2])
roc3$auc
#����ٽ������
coord1=coords(roc1, "best", ret=c("threshold", "specificity", "sensitivity"), transpose = FALSE)
coord2=coords(roc2, "best", ret=c("threshold", "specificity", "sensitivity"), transpose = FALSE)
coord3=coords(roc3, "best", ret=c("threshold", "specificity", "sensitivity"), transpose = FALSE)
#��ͼ
g1<-ggroc(list(svm_model=roc1,randomforest_model=roc2,neuralnet_model=roc3),legacy.axes = TRUE)
g2<-g1+annotate("text",x=.75,y=.45,label=paste("AUC of SVM =", round(roc1$auc,2)))+
  annotate("text",x=.75,y=.35,label=paste("AUC of randomforest =", round(roc2$auc,2)))+
  annotate("text",x=.75,y=.25,label=paste("AUC of neuralnet =", round(roc3$auc,2)))+
  geom_point(aes(x=1-coord1[1,2],y=coord1[1,3]),color="#9F5F9F",size=4,shape=16,alpha= .25)+
  geom_point(aes(x=1-coord2[1,2],y=coord2[1,3]),color="#9F5F9F",size=4,shape=16,alpha= .25)+
  geom_point(aes(x=1-coord3[1,2],y=coord3[1,3]),color="#9F5F9F",size=4,shape=16,alpha= .25)
g2
