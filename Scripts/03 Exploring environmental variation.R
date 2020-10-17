# Loading packages #
library(ggplot2)
library(gridExtra)
library(reshape)
library(vegan)

# Importing the dataset #
trawlings<-read.table("data/trawling data.txt",header=T)
avg_salinity<-read.table("data/avg. salinity.txt",header=T)
trawlings<-merge(trawlings,avg_salinity)

# Renaming species for further plotting #
trawlings$species<-as.factor(trawlings$species)
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
anova(lm(log(trawlings$sal_avg)~trawlings$species)) # significant p-value
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
  geom_boxplot(aes(x=trawlings$species,y=trawlings$sal_avg),
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
env<-aggregate(trawlings[,c(7,9,11)],list(trawlings$trawl),mean)
row.names(env)<-env$Group.1
env<-env[,-1]
pc.env<-princomp(decostand(env,method="standardize"),cor=TRUE)
eigenvals(pc.env) #Selecting first axis by Kaiser-Guttman criterion
summary(pc.env) #First axis explains 53.4% of total variation
pc.env$loadings #Depth and salinity are positively related to the first PC
scores<-pc.env$scores
loadings<-as.data.frame(pc.env$loadings[c(1:3),])
write.table(scores,"results/env_pcs.txt")


# Is environment spatially constrained? #
# PCA for environmental variables #
env<-aggregate(trawlings[,c(6,5,7:9)],list(trawlings$trawl),mean)
row.names(env)<-env$Group.1
env<-env[,-1]

library(geosphere)
distlin<-distm(env[,1:2])/1000
library(vegan)
mantel(distlin,vegdist(decostand(env,method="standardize"),method="euclidean"))
