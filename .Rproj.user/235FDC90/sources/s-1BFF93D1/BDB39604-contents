### Pacotes utilizados ###
library(ggplot2)
library(gridExtra)

### Figure 4 ###
tab0<-presabs[,c("Aquadriscutis","Sparkeri")]
x<-pairs[pairs$Sp1=="Aquadriscutis"&pairs$Sp2=="Sparkeri",]
tab0$id<-apply(tab0,1,sum)     # identification of the four classes
tab0$id[tab0$id == 2] <- "Double presence"
tab0$id[tab0$id == 0] <- "Double absence"
tab0$id[tab0[[1]]==1&tab0[[2]]==0] <- "co10"
tab0$id[tab0[[1]]==0&tab0[[2]]==1] <- "co01"
tab0<-cbind(presabs[1],tab0)  # association of the site id to the smaller table  
Dist.tab0<-merge(Dist,tab0, by.y="plot",all = FALSE)
Dist.tab0<-droplevels(Dist.tab0[!Dist.tab0$id=='co01',])
Dist.tab0<-droplevels(Dist.tab0[!Dist.tab0$id=='co10',])
boxplot(Dist.tab0$AvDist~Dist.tab0$id)

AquSpa<-ggplot()+
  geom_boxplot(aes(x=Dist.tab0$id,y=Dist.tab0$AvDist),fill="grey")+
  ylab("FDis")+
  xlab("")+
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


tab1<-presabs[,c("Arugispinis","Cspixii")]
x<-pairs[pairs$Sp1=="Arugispinis"&pairs$Sp2=="Cspixii",]
tab1$id<-apply(tab1,1,sum)     # identification of the four classes
tab1$id[tab1$id == 2] <- "Double presence"
tab1$id[tab1$id == 0] <- "Double absence"
tab1$id[tab1[[1]]==1&tab1[[2]]==0] <- "co10"
tab1$id[tab1[[1]]==0&tab1[[2]]==1] <- "co01"
tab1<-cbind(presabs[1],tab1)  # association of the site id to the smaller table  
Dist.tab1<-merge(Dist,tab1, by.y="plot",all = FALSE)
Dist.tab1<-droplevels(Dist.tab1[!Dist.tab1$id=='co01',])
Dist.tab1<-droplevels(Dist.tab1[!Dist.tab1$id=='co10',])
boxplot(Dist.tab1$AvDist~Dist.tab1$id)

AruCar<-ggplot()+
  geom_boxplot(aes(x=Dist.tab1$id,y=Dist.tab1$AvDist),fill="grey")+
  ylab("")+
  xlab("Co-occurrence patterns")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 16))


figure4<-grid.arrange(AquSpa,AruCar,nrow=2)
figure4
ggsave("Figure 4.png",figure4,dpi=600,height=9,width=8,units=c("cm"))
