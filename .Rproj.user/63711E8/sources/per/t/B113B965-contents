# Loading packages #
library(ggplot2)
library(gridExtra)
library(reshape)
library(vegan)

# Importing the dataset #
trawlings<-read.table("trawling data.txt",header=T)

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
  geom_boxplot(aes(x=trawlings$species,y=log(trawlings$prof)),
               notch=TRUE,fill="grey")+
  xlab("Species")+
  ylab("Log of Depth (m)")+
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
  ylab("Salinity (?)")+
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
ggsave("Figure 3.png",figure3,dpi=600,height=24,width=16,units=c("cm"))

# PCA for environmental variables #
env<-aggregate(trawlings[,7:9],list(trawlings$trawl),mean)
row.names(env)<-env$Group.1
env<-env[,-1]
pc.env<-princomp(decostand(env,method="standardize"),cor=TRUE)
eigenvals(pc.env) #Selecting first axis by Kaiser-Guttman criterion
summary(pc.env) #First axis explains 53.4% of total variation
pc.env$loadings #Depth and salinity are positively related to the first PC

# Plotting PCA #