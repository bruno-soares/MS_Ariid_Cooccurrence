### Pacotes utilizados ###
library(reshape)
library(vegan)
library(ade4)
library(Rmisc)

### Importing dataset and building Site x Sp. matrix ###
data<-read.table("data/trawling data.txt",header=T)
comm<-cast(data[,c(1:3)],trawl~species,value='pres',fun.aggregate=mean)
comm[is.na(comm)]<-0
row.names(comm)<-comm$trawl
comm<-comm[,-1]

### Reading C-score index function ###
# This function was retrieved from D'Amen et al. (2017). If you use this function, cite:
# D'Amen, Gotelli & Guisan (2017) Disentangling biotic interactions, environmental filters, and dispersal limitation as drivers of species co-occurrence. Ecography, 41 (8): 1233-1244. ##
ecospat.Cscore01 <- function(data.in,npermut,outpath){
  
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


### Calculating C-Score ###
data.in <- comm
nperm <- 10000
outpath <- "C:/Users/soare/OneDrive/Atual/Artigos em andamento/ariidae/MS_Ariid_Cooccurrence/results"
res <- ecospat.Cscore01(comm,nperm,outpath)
print(res)
