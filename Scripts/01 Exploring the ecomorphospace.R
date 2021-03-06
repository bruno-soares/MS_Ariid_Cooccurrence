# Opening libraries #
library(cluster)
library(vegan)
library(ggplot2)
library(gridExtra)

# Importing the dataset #
eco<-read.table("data/ecomorphological data.txt",header=T,row.names=1)
eco[,c(1:6)]<-decostand(eco[,c(1:6)],method="range")
eco$ST<-as.factor(eco$ST)

# Generating the PCoA #
pcoa<-cmdscale(daisy(eco,metric="gower"),k=ncol(eco)-1,eig=TRUE,add=TRUE)
scores<-as.data.frame(pcoa$points)

# Computing % of explanation of each axis #
(eigenvals(pcoa)/sum(eigenvals(pcoa)))*100

# Fitting ecomorphological indices to ordination #
efit<-envfit(pcoa,eco, permutations = 999, na.rm = TRUE,choices=c(1:3))
efit2<-as.data.frame(efit$vectors$arrows*sqrt(efit$vectors$r))
efit2<-rbind(efit2,efit$factors$centroids*sqrt(efit$factors$r))
efit

# Plotting Supplementary Figure 1 #
fig1<-ggplot()+
  geom_text(mapping=aes(x=scores[,1],y=scores[,2],label=rownames(scores)),size=5)+
  xlab("PC1 (33.21%)")+xlim(c(min(scores[,1])-0.08,max(scores[,1])+0.08))+
  ylab("PC2 (26.87%)")+ylim(c(min(scores[,2])-0.08,max(scores[,2])+0.08))+
  geom_segment(data=efit2,aes(x=0,xend=efit2[,1]/2.5,y=0,yend=efit2[,2]/2.5),
               arrow = arrow(length = unit(0.5, "cm")),colour="blue",
               size=0.8,alpha=0.2,inherit.aes=FALSE)+
  geom_text(data=efit2,aes(x=(efit2[,1]/2.5),y=(efit2[,2]/2.5+0.015),label=row.names(efit2)),
            size=4,color="blue",fontface="bold",alpha=0.5)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(colour = "black", size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 12))

fig2<-ggplot()+
  geom_text(mapping=aes(x=scores[,1],y=scores[,3],label=rownames(scores)),size=5)+
  xlab("PC1 (33.21%)")+xlim(c(min(scores[,1])-0.08,max(scores[,1])+0.08))+
  ylab("PC3 (19.95%)")+ylim(c(min(scores[,3])-0.08,max(scores[,3])+0.08))+
  geom_segment(data=efit2,aes(x=0,xend=efit2[,1]/2.5,y=0,yend=efit2[,3]/2.5),
               arrow = arrow(length = unit(0.5, "cm")),colour="blue",
               size=0.8,alpha=0.2,inherit.aes=FALSE)+
  geom_text(data=efit2,aes(x=efit2[,1]/2.5,y=efit2[,3]/2.5+0.02,label=row.names(efit2)),
            size=4,color="blue",fontface="bold",alpha=0.5)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(colour = "black", size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 12))

fig3<-ggplot()+
  geom_text(mapping=aes(x=scores[,2],y=scores[,3],label=rownames(scores)),size=5)+
  xlab("PC2 (26.87%)")+xlim(c(min(scores[,2])-0.08,max(scores[,2])+0.08))+
  ylab("PC3 (19.95%)")+ylim(c(min(scores[,3])-0.08,max(scores[,3])+0.08))+
  geom_segment(data=efit2,aes(x=0,xend=efit2[,2]/2.5,y=0,yend=efit2[,3]/2.5),
               arrow = arrow(length = unit(0.5, "cm")),colour="blue",
               size=0.8,alpha=0.2,inherit.aes=FALSE)+
  geom_text(data=efit2,aes(x=efit2[,2]/2.5+0.015,y=efit2[,3]/2.5+0.015,label=row.names(efit2)),
            size=4,color="blue",fontface="bold",alpha=0.5)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(colour = "black", size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 12))

sfig1<-grid.arrange(fig1,fig2,fig3,nrow=2)
sfig1
ggsave("Figures/Supplfig1.png",sfig1,dpi=600,height=18,width=20,units=c("cm"))
dev.off()
