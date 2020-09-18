# Loading packages #
library(ggplot2)
library(gridExtra)
library(reshape)
library(vegan)

# Importing the dataset #
trawlings<-read.table("data/trawling data.txt",header=T)

# Renaming species for further plotting #
levels(trawlings$species)[1]<-"Aph"
levels(trawlings$species)[2]<-"Aqu"
levels(trawlings$species)[3]<-"Aru"
levels(trawlings$species)[4]<-"Bba"
levels(trawlings$species)[5]<-"Car"
levels(trawlings$species)[6]<-"Ngr"
levels(trawlings$species)[7]<-"Sco"
levels(trawlings$species)[8]<-"Spa"
levels(trawlings$species)[9]<-"Spr"

# Average differences in sites that each species occurred #
anova(lm(log(trawlings$prof)~trawlings$species)) # significant p-value
anova(lm(log(trawlings$sal)~trawlings$species)) # significant p-value
anova(lm(log(trawlings$temp)~trawlings$species)) # non-significant p-value

# Plotting Figure 3 #
Depth<-ggplot()+
  geom_boxplot(aes(x=trawlings$species,y=trawlings$prof),
               notch=TRUE,fill="grey")+
  xlab("Species")+
  ylab("Depth (m)")+
  scale_y_log10(limits = c(5,60))+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 16))

Salinity<-ggplot()+
  geom_boxplot(aes(x=trawlings$species,y=trawlings$sal),
               notch=TRUE,fill="grey")+
  ylab("Salinity")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 16))

Temperature<-ggplot()+
  geom_boxplot(aes(x=trawlings$species,y=trawlings$temp),
               notch=TRUE,fill="grey")+
  xlab("Species")+
  ylab("Temperature (ºC)")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 16))

figure3<-grid.arrange(Depth,Salinity,Temperature,nrow=3)
ggsave("Figures/Figure 3.png",figure3,dpi=600,height=24,width=16,units=c("cm"))

# PCA for environmental variables #
env<-aggregate(trawlings[,7:9],list(trawlings$trawl),mean)
row.names(env)<-env$Group.1
env<-env[,-1]
pc.env<-princomp(decostand(env,method="standardize"),cor=TRUE)
eigenvals(pc.env) #Selecting first axis by Kaiser-Guttman criterion
summary(pc.env) #First axis explains 53.4% of total variation
pc.env$loadings #Depth and salinity are positively related to the first PC
scores<-pc.env$scores
loadings<-as.data.frame(pc.env$loadings[c(1:3),])
write.table(scores,"results/env_pcs.txt")

# Plotting PCA #
Suppl.Fig2<-ggplot()+
  geom_point(mapping=aes(x=scores[,1],y=scores[,2]),size=2,alpha=0.15)+
  xlab("PC1 (53.41%)")+  ylab("PC2 (32.93%)")+
  geom_segment(aes(x=0,xend=loadings[,1]*2,y=0,yend=loadings[,2]*2),
               arrow = arrow(length = unit(0.5, "cm")),colour="blue",
               size=0.8,alpha=0.5,inherit.aes=FALSE)+
  geom_text(aes(x=(loadings[,1]*2+0.3),y=(loadings[,2]*2+0.16),label=c("Depth","Salinity","Temperature")),
            size=4,color="blue",fontface="bold",alpha=0.5)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 14, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 10))

ggsave("results/Suppl. Fig. 2.png",Suppl.Fig2,dpi=600,height=10,width=10,units=c("cm"))
