### Importing and formatting dataset ###
# This script was retrieved from D'Amen et al. (2017). If you use it, cite:
# D'Amen, Gotelli & Guisan (2017) Disentangling biotic interactions, environmental filters, and dispersal limitation as drivers of species co-occurrence. Ecography, 41 (8): 1233-1244. ##
nperm <- 10000
outpath <- getwd()
coocc.tab.all<-read.table("Cscores01.txt",h=T)
obs.r<-data.frame(coocc.tab.all$obs.C.score) # Observed C-scores
pairs<-coocc.tab.all[,1:2] # Names of significant pairs
CooccProb<-read.table("MatrixPermutations01.txt", h=T)  # Null C-scores
CooccProb.pairs<-cbind(pairs,CooccProb) # Attribute names of pairs to Null values
CooccProb.pairs.mean<-apply(CooccProb.pairs[,3:nperm],1,mean) # Calculate mean of null values
CooccProb.pairs.sd<-apply(CooccProb.pairs[,3:nperm],1,sd) # Calculate standard deviation of null values
CooccProb.stat<-data.frame(pairs,CooccProb.pairs.mean,CooccProb.pairs.sd)  

### Assign each C-score to 22 evenly spaced bins ###
n.bin<-22
colnames(obs.r)<-"Cscore"
obs.r$bin<- cut(obs.r$Cscore, breaks=seq(0,1,by=1/n.bin),include.lowest=TRUE,labels=1:n.bin)  
obs.r<-cbind(obs.r,pairs)

### Calculate how many pairs are in each bin ###
obs.r1<-obs.r[order(obs.r[,2], obs.r[,1]),] 
obs.bin.summary<-as.data.frame(table(obs.r1$bin)) 
obs.bin.meanCS<-aggregate(obs.r1$Cscore, list(obs.r1$bin), mean)
names(obs.bin.summary)<-c("bin","n.pairs")

### Calculate the null expectation of C-scores in each bin ###
matrix.null.bins<-CooccProb 
matrix.null.bin.meanCS<-data.frame(seq(1:n.bin))
names(matrix.null.bin.meanCS)<-"Group.1"

for (i in 1:ncol(CooccProb)){
  null.bin<-cut(CooccProb[,i], breaks=seq(0,1,by=1/n.bin), include.lowest=TRUE, labels=1:n.bin)
  matrix.null.bins[,i]<-null.bin
  null.bin1<-cbind(CooccProb[,i],null.bin)   
  null.bin.meanCS<-aggregate(null.bin1[,1],list(null.bin1[,2]),mean)
  matrix.null.bin.meanCS<-merge(matrix.null.bin.meanCS, null.bin.meanCS, by="Group.1", all.x=T)}

nperm<-10000
tab.summary<-data.frame(matrix(data=NA,byrow=T,nrow=n.bin,ncol=nperm)) #table with the number of pairs in each bin for each null community

for (z in 1:nperm){                                                   
  col.summary<-as.data.frame(table(matrix.null.bins[,z]))
  tab.summary[,z]<-col.summary[,2]
  colnames(tab.summary)[z]<-z} 

matrix.null.bin.meanCS[is.na(matrix.null.bin.meanCS)] <- 0
tab.summary.meanCS<-apply(matrix.null.bin.meanCS[,2:(nperm+1)],1,mean, na.rm=T)
null.meanCS<-cbind(seq(1:n.bin),tab.summary.meanCS)
names(null.meanCS)<-c("bin", "null.meanCS")

### Calculate the mean and 95% confidence interval for each bin ###
myci <- function(t) {                             
  n <- length(t) 
  se <- sd(t)/sqrt(n) 
  m <- mean(t) 
  cv <- qt(0.975,df=n-1) 
  c(m-cv*se,m+cv*se) 
}

mean.pairs.bin<-(apply(tab.summary,1,mean))
CI.pairs.bin<-t(apply(tab.summary,1,myci))
sd.pairs.bin<-apply(tab.summary,1,sd)
stat.pairs.bin<-data.frame(seq(1:n.bin),mean.pairs.bin,CI.pairs.bin)
names(stat.pairs.bin)<-c("bin","MEAN","UP.CI","LOW.CI")
half.bin<-signif(((1/n.bin)/2),2)
bin.mean<-seq(0,1,by=1/n.bin)-half.bin
bin.mean<-bin.mean[2:(n.bin+1)]

# Table with observed number of pairs for each bin #
obs.bin.summary2<-cbind(obs.bin.summary,bin.mean)

# Table with mean and CI null number of pairs for each bin #
stat.pairs.bin2<-cbind(stat.pairs.bin,bin.mean)

### Comparison of observed (black) and null (white) number of pairs in each bin ###
plot(x=stat.pairs.bin2[,5],y=stat.pairs.bin2[,2], xlab="C-score",ylab="Number of pairs")
points(y=obs.bin.summary2[,2],x=obs.bin.summary2[,3], pch=16) 

## Selecting the species pairs of the bins which the observed number of significant pairs
# were more than the mean number of significant pairs
tab.finalBayesM<-data.frame(matrix(ncol=ncol(obs.r1),nrow=0))
names(tab.finalBayesM)<-names(obs.r1)
for(j in 1:n.bin){
  obs.r2<-obs.r1[obs.r1$bin==j,]
  if((obs.bin.summary[j,2]>stat.pairs.bin[j,2])==T){   # if the observed number of species in the bin is higher 
    # than the mean number of species expected by null communities in the same bin
    exp.mean<-as.numeric(null.meanCS[j,2])             # we retain the pairs with a C-score higher than the mean C-score measured from null community in that bin
    obs.r3<-obs.r2[obs.r2$Cscore>exp.mean,]
    tab.finalBayesM<-rbind(tab.finalBayesM, obs.r3) 
  }}

### Retain only pairs statistically significant in an individual test ###
BayesM_merge<-merge(tab.finalBayesM,coocc.tab.all,by.y=c("Sp1","Sp2"))              
sign.BayesM<-BayesM_merge[which(BayesM_merge$pval_less<0.05|BayesM_merge$pval_greater<0.05),]
write.table(sign.BayesM,"Sign.BayesM.txt", sep="\t")