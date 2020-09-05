library(reshape)
library(ggplot2)
library(funrar)
library(cluster)
library(vegan)
library(ape)
library(gridExtra)
abiot<-read.table("abiot.txt",header=T)
abiot

levels(abiot$NAME_OF_SPECIE)
levels(abiot$NAME_OF_SPECIE)[3]<-"Aru"
levels(abiot$NAME_OF_SPECIE)[6]<-"Ngr"
levels(abiot$NAME_OF_SPECIE)[1]<-"Aph"
levels(abiot$NAME_OF_SPECIE)[5]<-"Car"
levels(abiot$NAME_OF_SPECIE)[8]<-"Spa"
levels(abiot$NAME_OF_SPECIE)[9]<-"Spr"
levels(abiot$NAME_OF_SPECIE)[2]<-"Aqu"
levels(abiot$NAME_OF_SPECIE)[4]<-"Bba"
levels(abiot$NAME_OF_SPECIE)[7]<-"Sco"
plot(abiot$Prof~abiot$NAME_OF_SPECIE)
plot(abiot$Sal.Fund~abiot$NAME_OF_SPECIE)
plot(abiot$Temp.Fund~abiot$NAME_OF_SPECIE)

Depth<-ggplot()+
  geom_boxplot(aes(x=abiot$NAME_OF_SPECIE,y=log(abiot$Prof)),
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
  geom_boxplot(aes(x=abiot$NAME_OF_SPECIE,y=abiot$Sal.Fund),
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
  geom_boxplot(aes(x=abiot$NAME_OF_SPECIE,y=abiot$Temp.Fund),
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

library(gridExtra)
figure3<-grid.arrange(Depth,Salinity,Temperature,nrow=3)
ggsave("Figure 3.png",figure2,dpi=600,height=24,width=16,units=c("cm"))

anova(lm(log(abiot$Prof)~abiot$NAME_OF_SPECIE))
anova(lm(log(abiot$Sal.Fund)~abiot$NAME_OF_SPECIE))
anova(lm(log(abiot$Temp.Fund)~abiot$NAME_OF_SPECIE))