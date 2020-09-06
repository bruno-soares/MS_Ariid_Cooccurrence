### Pacotes utilizados ###
library(reshape)
library(vegan)
library(ade4)
library(Rmisc)

### Dados brutos e construção das planilhas ###
data<-read.table("database.txt",header=T)
comm<-cast(data[,c(1:3)],trawl~species,value='pres',fun.aggregate=mean)
comm[is.na(comm)]<-0
row.names(comm)<-comm$trawl
comm<-comm[,-1]

### Reading function: Pairwise C-score ###
## Format required: a plots (rows) x species (columns) matrix of presences/absences
## The function c.score calculates the C-score matrix to detect species association, for the whole community and for species pairs
## Randomization: column sum is fixed
## It returns the C-score index for the observed community (ObsCscoreTot), p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the working directory where the same 
## metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version)
## NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002 - Ecology)
# ecospat.Cscore01(data.in, 100 ,outpath)

ecospat.Cscore01 <- function(data.in,npermut,outpath)
  
{
  
  # C-coef Observed matrix
  cat("Computing observed co-occurence matrix", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)  
  
  
  spec.occ <- data.matrix(data.in)
  
  ####  C score #############	
  coocc<-t(spec.occ)%*%spec.occ  # nb of checkboard units
  n.spec=dim(coocc)[1]			 
  mat1<-array(apply(spec.occ,MAR=2,sum),dim=c(n.spec,n.spec)) 
  mat2<-t(array(apply(spec.occ,MAR=2,sum),dim=c(n.spec,n.spec)))
  mat.obs.c.coef <- ((mat1 - coocc)*(mat2 - coocc))/(mat1*mat2) # observed c score 
  df.obs.c.coef <- data.frame(Col = rep(1:ncol(mat.obs.c.coef),each=ncol(mat.obs.c.coef)),Row = rep(1:nrow(mat.obs.c.coef),
                                                                                                    nrow(mat.obs.c.coef)),Sp1 = rep(colnames(mat.obs.c.coef),each=ncol(mat.obs.c.coef)),
                              Sp2 = rep(rownames(mat.obs.c.coef),nrow(mat.obs.c.coef)),Co.Occ = c(mat.obs.c.coef)) # dataframe with cscore for each species pair
  v.diago.inf <- c(rownames(df.obs.c.coef)[df.obs.c.coef[,1]>df.obs.c.coef[,2]],rownames(df.obs.c.coef)[df.obs.c.coef[,1]==df.obs.c.coef[,2]])# Remove identical combinations of species 
  df.obs.c.coef <- df.obs.c.coef[-as.numeric(v.diago.inf),]
  CscoreTot<-mean(df.obs.c.coef$Co.Occ)
  
  # Matrix to store the permuations
  mat.perm <- matrix(0,nrow(df.obs.c.coef),npermut, dimnames = list(c(paste(df.obs.c.coef[,3],df.obs.c.coef[,4])),c(1:npermut)))
  
  
  # Permutations C-coef 
  
  cat("Computing permutations", "\n",append = F)
  cat(".............", "\n",append = F)
  
  for (i in 1:npermut)
  {
    if (i == 1)
    {
      cat(npermut ," permutations to go", "\n",append = F)
      cat(".............", "\n",append = F)			
    }	
    if (i == npermut / 2)
    {
      cat(npermut / 2," permutations to go", "\n",append = F)
      cat(".............", "\n",append = F)
    }	
    
    
    spec.occ.perm1<-data.matrix(data.in)
    spec.occ.perm1 <- permatswap(spec.occ.perm1,fixedmar="both",mtype="prab",time=1) # row/column sums are preserved
    # time=1 : separate swapping sequence that always begins with the original matrix
    spec.occ.perm <- as.matrix(spec.occ.perm1[[3]][[1]] )
    
    coocc.perm <- t(spec.occ.perm)%*%spec.occ.perm 
    mat1.perm <- array(apply(spec.occ.perm,MAR=2,sum),dim=c(n.spec,n.spec))
    mat2.perm <- t(array(apply(spec.occ.perm,MAR=2,sum),dim=c(n.spec,n.spec)))
    
    mat.obs.c.coef.perm <- ((mat1.perm - coocc.perm)*(mat2.perm - coocc.perm))/(mat1.perm*mat2.perm)
    
    df.obs.c.coef.perm <- data.frame(Col = rep(1:ncol(mat.obs.c.coef.perm),each=ncol(mat.obs.c.coef.perm)),Row = rep(1:nrow(mat.obs.c.coef.perm),
                                                                                                                     nrow(mat.obs.c.coef.perm)),Sp1 = rep(colnames(mat.obs.c.coef),each=ncol(mat.obs.c.coef.perm)), 
                                     Sp2 = rep(rownames(mat.obs.c.coef),nrow(mat.obs.c.coef.perm)),Co.Occ = c(mat.obs.c.coef.perm))
    
    # Remove identical combinations of species (same Co-occ coef) and the diagonal (Co-occ coeff = 0)	
    df.obs.c.coef.perm <- df.obs.c.coef.perm[-as.numeric(v.diago.inf),]
    
    # Store result of permuation
    mat.perm[,i] <- df.obs.c.coef.perm[,5]
    
  }
  
  
  ## for the whole community
  vec.CScore.tot<-as.vector(apply(mat.perm,MAR=2,mean)) # C-score for all null communities (mean on the columns)
  SimulatedCscore<-mean(vec.CScore.tot) # mean of Simulation C-score: Simulated C-score
  sd.SimulatedCscore<-sd(vec.CScore.tot) # standard deviation of null communities
  Zscore<-(CscoreTot-SimulatedCscore)/sd.SimulatedCscore # standardized effect size
  
  randtest.less<-as.randtest(vec.CScore.tot, CscoreTot, alter="less")
  pval.less<-randtest.less$pvalue
  randtest.greater<-as.randtest(vec.CScore.tot, CscoreTot, alter="greater")
  pval.greater<-randtest.greater$pvalue
  plot(randtest.greater, xlab= "Simulated C-scores",main=paste("", sep=""))
  
  
  # Calculate P-values based on random distribution
  mat.pval <- matrix(0,nrow(mat.perm),4,dimnames = list(rownames(mat.perm),c("Obs.Co.Occ","Zscore","pval_less","pval_greater")))
  mat.pval[,1] <- df.obs.c.coef[,5]
  
  cat("Computing P-values", "\n",append = F)
  cat(".............", "\n",append = F)
  
  for (k in 1:nrow(mat.perm))
  {
    
    mat.pval[k,2]<-	(df.obs.c.coef[k,5]-mean(mat.perm[k,]))/sd(mat.perm[k,])
    
    randtest<-as.randtest(sim=mat.perm[k,], obs=df.obs.c.coef[k,5], alter="less")
    mat.pval[k,3]<-randtest$pvalue
    randtest<-as.randtest(sim=mat.perm[k,], obs=df.obs.c.coef[k,5], alter="greater")
    mat.pval[k,4]<-randtest$pvalue        
  }
  
  # Exporting Co-occ matrix
  cat("Exporting dataset", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)  
  
  hist(as.vector(mat.pval[,2]), xlab="Zscore", main = paste(""))
  abline(v=c(2,-2),col = "red")
  
  mat.pval.names<-data.frame(df.obs.c.coef[,3:4],mat.pval,df.obs.c.coef.perm[,5])
  mat.pval.names2<-data.frame(mat.pval.names[,1:3],mat.pval.names[,7],mat.pval.names[,4:6])
  names(mat.pval.names2)[3]<-"obs.C-score"
  names(mat.pval.names2)[4]<-"exp.C-score"
  write.table(mat.pval.names2,file=paste(outpath,"\\Cscores01.txt", sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  
  tab<-mat.pval.names2
  v<-c(0)
  for (i in 1:nrow(tab)){
    if (tab[i,6]<=0.05||tab[i,7]<=0.05){
      v<-c(v,i)
    }
  }
  m<-data.frame()
  for(j in 1:length(v)){
    m<-rbind(m,tab[v[j],])
  }
  
  m1<-na.omit(m)
  
  write.table(m1,file=paste(outpath,"\\Sign_Cscores01.txt",sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  write.table(mat.perm, file=paste(outpath,"\\MatrixPermutations01.txt",sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  
  l<-list(ObsCscoreTot=CscoreTot, SimCscoreTot=SimulatedCscore, PVal.less=pval.less,PVal.greater=pval.greater,Z.score=Zscore)
  
  return(l)
  
  cat("Computations finished!", "\n",append = F)
  
  # END FUNCTION
}


### C-score para os nossos dados ###
data.in <- comm
nperm <- 10000
outpath <- getwd()
res <- ecospat.Cscore01(comm,nperm,outpath)
print(res)



## Bayesian identification of significant species pairs ##
coocc.tab.all<-read.table("Cscores01.txt",h=T)
obs.r<-data.frame(coocc.tab.all$obs.C.score)                          # Values of observed C-scores
pairs<-coocc.tab.all[,1:2]                                            # species pairs names
CooccProb<-read.table("MatrixPermutations01.txt", h=T)     		# matrix with C-scores values for null communities
CooccProb.pairs<-cbind(pairs,CooccProb)                               # attribute species pairs name to these values
CooccProb.pairs.mean<-apply(CooccProb.pairs[,3:nperm],1,mean)         # calculate for each pair the mean and sd C-score across null pairs
CooccProb.pairs.sd<-apply(CooccProb.pairs[,3:nperm],1,sd)
CooccProb.stat<-data.frame(pairs,CooccProb.pairs.mean,CooccProb.pairs.sd)  
## 2 Assign each pairwise C-score to one of 22 evenly spaced bins spanning the interval from 0 to 1.
n.bin<-22
colnames(obs.r)<-"Cscore"
obs.r$bin<- cut(obs.r$Cscore, breaks=seq(0,1,by=1/n.bin),include.lowest=TRUE,labels=1:n.bin)  
obs.r<-cbind(obs.r,pairs)
## calculates how many pairs are present in each bin for oberved communities
obs.r1<-obs.r[order(obs.r[,2], obs.r[,1]),] 
obs.bin.summary<-as.data.frame(table(obs.r1$bin)) 
obs.bin.meanCS<-aggregate(obs.r1$Cscore, list(obs.r1$bin), mean)
names(obs.bin.summary)<-c("bin","n.pairs")

## 3.## Calculation of the average number of species pairs from the null communities with different scores in each bin (table null.meanCS)
## This average represents the null expectation of the C-score for species pairs in each bin. 
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

## 4.## Within each bin, calculation of the mean and the 95% confidence limit 
myci <- function(t) {                             # Funtion to calculate the confidence intervals (CI)
  n <- length(t) # n is the sample size
  se <- sd(t)/sqrt(n) # Find the standard error of the sample
  m <- mean(t) # Find the sample mean
  cv <- qt(0.975,df=n-1) # cv is a critical value for the t distribution. P( t > cv ) = 0.025 = P( t < -cv )
  c(m-cv*se,m+cv*se) # Return the 95% confidence interval
}

mean.pairs.bin<-(apply(tab.summary,1,mean))
CI.pairs.bin<-t(apply(tab.summary,1,myci))
sd.pairs.bin<-apply(tab.summary,1,sd)

stat.pairs.bin<-data.frame(seq(1:n.bin),mean.pairs.bin,CI.pairs.bin)
names(stat.pairs.bin)<-c("bin","MEAN","UP.CI","LOW.CI")

half.bin<-signif(((1/n.bin)/2),2)
bin.mean<-seq(0,1,by=1/n.bin)-half.bin
bin.mean<-bin.mean[2:(n.bin+1)]

#table with OBSERVED n pairs for each bin
obs.bin.summary2<-cbind(obs.bin.summary,bin.mean)
#table with MEAN and CI n pairs for each bin
stat.pairs.bin2<-cbind(stat.pairs.bin,bin.mean)

## Graph comparing how many pairs are present in each bin in average across the null communities (white cirles)
## and the number of observed pairs in the same bin (black circle)
plot(x=stat.pairs.bin2[,5],y=stat.pairs.bin2[,2], xlab="C-score",ylab="Number of pairs")
points(y=obs.bin.summary2[,2],x=obs.bin.summary2[,3], pch=16) 

## 5. ## Within each bin, order the species pairs by their scores and retain the species pairs with the largest scores that
## place them above the mean for the number of species pairs expected from the simulated distribution (Bayes M criterion).
## This step can be modified to apply the Bayes CL criterion using the values calculated above for the 95% confidence limit
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


## 6. ## Further reduce the set of significant pairs by retaining only those that 
## are statistically significant in an individual test (simple CL criterion). 

BayesM_merge<-merge(tab.finalBayesM,coocc.tab.all,by.y=c("Sp1","Sp2"))              
sign.BayesM<-BayesM_merge[which(BayesM_merge$pval_less<0.05|BayesM_merge$pval_greater<0.05),]

write.table(sign.BayesM,"Sign.BayesM.txt", sep="\t")



## Blois framework ##
###	Em todas as planilhas, a coluna com os sítios deve se chamar plot
#### LOAD NECESSARY DATA
presabs<-cast(data[,c(1:3)],trawl~species,value='pres',fun.aggregate=mean)
presabs[is.na(presabs)]<-0
colnames(presabs)[1]<-"plot"
sign.BayesM$obs.C.score-sign.BayesM$exp.C.score
ifelse((sign.BayesM$obs.C.score-sign.BayesM$exp.C.score)>0,"SEGR","AGGR")

pairs<-sign.BayesM[,c(1:2)]
pairs<-cbind(pairs,ifelse((sign.BayesM$obs.C.score-sign.BayesM$exp.C.score)>0,"SEGR","AGGR"))
colnames(pairs)[3]<-"NullModel"
# column 1 ="Sp1", column 2 ="Sp2", 
# column 3 = "NullModel" containing null model result, either "SEGR" or "AGGR" for each pair
#SEGR = obs > exp; AGGR = obs < exp

env<-aggregate(data[,7:9],list(data$trawl),mean)
row.names(env)<-env$Group.1
colnames(env)[1]<-"plot"

coor<-aggregate(data[,5:6],list(data$trawl),mean)
row.names(coor)<-coor$Group.1
colnames(coor)[1]<-"plot"

Dist<-read.table("avDist.txt",header=T,row.names=1)
Dist<-cbind(Dist,coor$plot)
colnames(Dist)[2]<-"plot"

#### PCA for the environmental variables
pc.env<-princomp(~prof+sal+temp,cor=TRUE,data=env)
bstick(pc.env)
eigenvals(pc.env) #just first axis
summary(pc.env)
comp.env<-data.frame(env[,1],pc.env$scores[,1]) #Modifiquei para data.frame ao inv?s de cbind
colnames(comp.env)[1]<-"plot"

#### preparation of the columns for store the results from the tests
pairs$coor.test<-NA
pairs$env.test<-NA
pairs$r.env.test<-NA
pairs$Dist.test<-NA
pairs$r.Dist.test<-NA
pairs

#### For each species pairs identification of the four co-occurrence classes in each site ("co00"= both absent, "co11"=both present
#### "co01" and "co10"= checkerboard distributions) and spatial configuration and environmental tests

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
  
  
  #### Test of functional arrangement: ANOVA on the distinctiveness
  Dist.tab1<-merge(Dist, tab1, by.y="plot",all = FALSE)
  fac1=factor(Dist.tab1$id)
  xx=as.matrix(Dist.tab1[,2])
  anova.Dist<-aov(xx~fac1)
  sum.Dist<-summary(anova.Dist)
  boxplot(xx~fac1)
  pairs$Dist.test[i]<-sum.Dist[[1]][[1,"Pr(>F)"]]
  pairs$r.Dist.test[i]<-summary(lm(xx~fac1))$r.squared
  
  #### Test of spatial arrangement: MANOVA on the coordinates
  coor.tab1<-merge(coor, tab1, by.y="plot",all = FALSE)
  fac1=factor(tab1$id)
  xx=as.matrix(coor.tab1[,2:3])
  m.coor=manova(xx~fac1)
  sum.coor=summary(m.coor)
  pairs$coor.test[i]<-sum.coor$stats[1,6] 
  
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