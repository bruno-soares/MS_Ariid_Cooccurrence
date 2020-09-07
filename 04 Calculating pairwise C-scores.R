### Pacotes utilizados ###
library(reshape)
library(vegan)
library(ade4)
library(Rmisc)

### Importing dataset and building Site x Sp. matrix ###
data<-read.table("trawling data.txt",header=T)
comm<-cast(data[,c(1:3)],trawl~species,value='pres',fun.aggregate=mean)
comm[is.na(comm)]<-0
row.names(comm)<-comm$trawl
comm<-comm[,-1]

install.packages("ecospat",dependencies=TRUE)
library(ecospat)