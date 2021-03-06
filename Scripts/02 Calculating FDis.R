# Opening libraries #
library(reshape)
library(funrar)

# Importing and managing the dataset #
data<-read.table("data/trawling data.txt",header=T)
comm<-cast(data[,c(1:3)],trawl~species,value='pres',fun.aggregate=mean)
comm[is.na(comm)]<-0
row.names(comm)<-comm$trawl
comm<-comm[,-1]

eco<-read.table("data/ecomorphological data.txt",header=T,row.names=1)
colnames(comm)<-rownames(eco)
comm<-as.data.frame(comm)
eco[,c(1:6)]<-decostand(eco[,c(1:6)],method="range")
eco$ST<-as.factor(eco$ST)

# COmputing average Functional distinctiveness sensu Grenié et al. (2017) per plot #
divfunc<-funrar(as.matrix(comm),as.matrix(daisy(eco,metric="gower"))
                ,rel_abund=FALSE)
distinct<-data.frame(divfunc$Di)
distinct[is.na(distinct)]<-0
avdist<-as.data.frame(apply(distinct,1,function(x) mean(x)))
colnames(avdist)[1]<-"AvDist"
write.table(avdist,file="./results/avDist.txt")
