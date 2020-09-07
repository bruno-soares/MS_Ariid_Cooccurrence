# This script was adapted from D'Amen et al. (2017). If you use it, cite:
# D'Amen, Gotelli & Guisan (2017) Disentangling biotic interactions, environmental filters, and dispersal limitation as drivers of species co-occurrence. Ecography, 41 (8): 1233-1244. ##

#### Loading and formatting dataset ###
data<-read.table("trawling data.txt",header=T)
presabs<-cast(data[,c(1:3)],trawl~species,value='pres',fun.aggregate=mean)
presabs[is.na(presabs)]<-0
colnames(presabs)[1]<-"plot"

sign.BayesM<-read.table("sign.BayesM.txt",header=T)
pairs<-sign.BayesM[,c(1:2)]
pairs<-cbind(pairs,ifelse((sign.BayesM$obs.C.score-sign.BayesM$exp.C.score)>0,"SEGR","AGGR"))
colnames(pairs)[3]<-"NullModel"

coor<-aggregate(data[,5:6],list(data$trawl),mean)
row.names(coor)<-coor$Group.1
colnames(coor)[1]<-"plot"

Dist<-read.table("avDist.txt",header=T,row.names=1)
Dist<-cbind(Dist,coor$plot)
colnames(Dist)[2]<-"plot"

env<-read.table("env_pcs.txt",header=T)
env<-cbind(env,coor$plot)
colnames(env)[4]<-"plot"

calculate_mode <- function(x) { #function to return the season of each trawl
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
season<-aggregate(data[,10],list(data$trawl),calculate_mode)
row.names(season)<-season$Group.1
colnames(season)[1]<-"plot"


### Preparing preparation of the columns for store the results from the tests  ###
pairs$coor.test<-NA
pairs$r.coor.test<-NA
pairs$env.test<-NA
pairs$r.env.test<-NA
pairs$Dist.test<-NA
pairs$r.Dist.test<-NA
pairs$season.test<-NA
pairs$r.Season.test<-NA
pairs

### Identify the four co-occurrence classes for each pair of species in each site:
# ("co00"= both absent, "co11"=both present, "co01" and "co10"= checkerboard distributions)
# and perform the importance of spatial configuration, environment, functional distinctiveness and seasonality

for (i in 1:nrow(pairs)) {
  sp1<-as.character(pairs[i,1])
  sp2<-as.character(pairs[i,2])
  
  tab<-presabs[,c(sp1,sp2)]                
  x<-pairs[pairs$Sp1==sp1&pairs$Sp2==sp2,] 
  
  tab$id<-apply(tab,1,sum)     # identification of the four classes
  tab$id[tab$id == 2] <- "co11"
  tab$id[tab$id == 0] <- "co00"
  tab$id[tab[[1]]==1&tab[[2]]==0] <- "co10"
  tab$id[tab[[1]]==0&tab[[2]]==1] <- "co01"
  
  tab<-cbind(presabs[1],tab)  # association of the site id to the smaller table  
  
  if (x$NullModel=="SEGR") {                             # identification of the pattern of segregation or
    tab1<-tab[tab[[4]]=="co10"|tab[[4]]=="co01",]          # aggregation for the considered species pairs
  }else if (x$NullModel=="AGGR") {                       # to select the sites to test in the following lines
    tab1<-tab[tab[[4]]=="co11"|tab[[4]]=="co00",]
  }
  
  #### Test of spatial arrangement: MANOVA on the coordinates
  coor.tab1<-merge(coor, tab1, by.y="plot",all = FALSE)
  fac1=factor(tab1$id)
  xx=as.matrix(coor.tab1[,2:3])
  m.coor=manova(xx~fac1)
  sum.coor=summary(m.coor)
  pairs$coor.test[i]<-sum.coor$stats[1,6] 
  
  #### Test of functional arrangement: ANOVA on the distinctiveness
  Dist.tab1<-merge(Dist, tab1, by.y="plot",all = FALSE)
  fac1=factor(Dist.tab1$id)
  xx=as.matrix(Dist.tab1[,2])
  anova.Dist<-aov(xx~fac1)
  sum.Dist<-summary(anova.Dist)
  boxplot(xx~fac1)
  pairs$Dist.test[i]<-sum.Dist[[1]][[1,"Pr(>F)"]]
  pairs$r.Dist.test[i]<-summary(lm(xx~fac1))$r.squared
  
  #### Test of environmental filter: MANOVA on the PCA axes
  env.tab1<-merge(comp.env,tab1, by.y="plot",all = FALSE) 
  fac1=factor(env.tab1$id)
  yy<-as.matrix(env.tab1[,2])
  m.env<-aov(yy~fac1)
  sum.env=summary(m.env)
  boxplot(yy~fac1)
  pairs$env.test[i]<-sum.env[[1]][[1,"Pr(>F)"]]
  pairs$r.env.test[i]<-summary(lm(yy~fac1))$r.squared
} 

pairs
write.table(pairs, "ResultBlois.txt", sep = "\t")