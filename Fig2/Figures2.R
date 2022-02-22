#libraries
library(ggplot2)
require(scales)
library(grid)
library(data.table)
#import data
#setwd("~/Documents/Deletion_paper/msa_0106/Figure_Proteins")
#############################
#all calculations
#############################
#####################################################
#input pre-calculated in Fig1
#####################################################
df0<-fread("~/Documents/Deletion_paper/msa_0106/Fig1/count_per_protein_pos.txt", sep="\t")
allp0<-meta%>% mutate(AA.Substitutions= strsplit(as.character(meta$AA.Substitutions), ",")) %>% 
  unnest(AA.Substitutions)

colnames(df0)<-c("V1","POS.AA","AA","Freq","Protein")
df0$Protein<-tolower(df0$Protein)
df<-df0[df0$V1 %like% "del",]
ins<-df0[df0$V1 %like% "ins",]
colnames(ins)<-c("V1","POS.AA","AA","Freq_ins","Protein")
subs<-df0[!(df0$V1 %like% "ins" | df0$V1 %like% "del"),]
colnames(subs)<-c("V1","POS.AA","AA","Freq_subs","Protein")


#df<- read.delim("Mutations_Protein.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
df1<- read.csv("Structure_stat.csv", sep = ",", header=TRUE) #COPY DATA DA proteins both up and down
df2<- read.csv("Structure_Models.csv", sep = ",", header=TRUE) #COPY DATA DA proteins both up and down
df12<-merge(df1,df2,by.x="X.pdbid",
            by.y="representative.PDB.chain", all.x=T )

anno<- read.delim("Annotations.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down

#spike Fig
len<-1273
annoS<-anno[anno$Protein=="S",]
S0<-df[df$Protein=="s",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
ymax<-max(S[!(is.na(S))])*100
ES<- read.delim("SE.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))

#add structure stat
SS<-df12[df12$Protein == "S",]
SS$Ex<-as.numeric((ifelse(SS$X..exposed.in.an.isolated.chain >= 30,ymax-(ymax/3),"NA")))
SS$Bu<-as.numeric((ifelse(SS$X..exposed.in.an.isolated.chain <= 20,ymax-(ymax/3),"NA")))
SS$Residue..protein<-as.numeric(SS$Residue..protein)
SSf<-merge(Sl, SS, by.x = "POS", by.y = "Residue..protein", all.x = T )

#bind
SF<-cbind(S,ES,SSf[,-1])
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
#ins subs merged
insS<-ins[ins$Protein=="s",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)

subsS<-subs[subs$Protein=="s",]
subsS<-aggregate(Freq_subs~POS.AA,subsS,sum)
SF<-merge(SF,subsS, by.x="POS", by.y="POS.AA", all.x=T)

#plot
PSpike<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=0.7, col="red")+
  geom_point(aes(y = Freq_ins), size=0.8, col="blue")+
  geom_point(aes(y = Freq_subs), size=0.6, col="darkgreen")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  #geom_line(aes(y=Ex), size=2, col="darkred")+
  #geom_line(aes(y=Bu), size=2, col="#0070C0")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(100, len, by = 100),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                     labels = c(1, 10, 100, 1000, 10000, 100000,1000000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
                     )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                             size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                   size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  annotate("segment", 
           x = 1, xend = len, y = ymax-(ymax*0.49), yend = ymax-(ymax*0.49),
           size = 2, col="black"
           )+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = c(rep(ymax-(ymax*0.49), 4),rep(ymax-(ymax*0.94), 2)), 
           yend = c(rep(ymax-(ymax*0.49), 4),rep(ymax-(ymax*0.94), 2)),
           colour = c("darkred","lightblue3","red",
                      "darkblue","darkorange2","purple"), 
           size =c(2,rep(5,5))
           )+
  annotate("segment", 
           x = 1, xend = 12, 
           y = ymax-(ymax*0.49), yend = ymax-(ymax*0.49),
           size = 2.5,
           col="black"
           #arrow=arrow(length = unit(3,'mm')
           )+
  annotate("text", 
         x = (annoS$xs+((annoS$xe-annoS$xs)/2)),
         hjust = 0.5,
         y = c(rep(ymax-(ymax*0.81), 4),rep(ymax-(ymax*0.98), 2)),
         colour = c("white","lightblue4","red",
                    "darkblue","darkorange2","purple") ,
         size = 7, label=annoS$Annotation,
         fontface = c(rep("bold",4),rep("bold.italic",2)),
         )+
  annotate("segment", 
           x = c(69,144,210,253),
           xend = c(69,144,210,253), 
           y = 1, yend = ymax,
           colour = "pink", 
           size=c(5,11,6,13),
           alpha=0.3)+
  annotate("text", 
           x = c(69,144,206,258),
           y= rep(c(ymax-(ymax*0.96),ymax-(ymax*0.98)),2),
           label = c("S-HVR1","S-HVR2","S-HVR3","S-HVR4"), 
           colour = "black", 
           fontface="bold",
           size=4
           )

PSpike
pdf("Spike.pdf", width = 16, height = 5)
PSpike
dev.off()

#NSP1
len<-180
annoS<-anno[anno$Protein=="nsp1",]
S0<-df[df$Protein=="nsp1",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
ymax<-max(S[!(is.na(S))])*100
ES<- read.delim("NSP1E.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
#add structure stat
SS<-df12[df12$Protein == "nsp1",]
SS$Ex<-as.numeric((ifelse(SS$X..exposed.in.an.isolated.chain >= 50,ymax-(ymax/3),"NA")))
SS$Bu<-as.numeric((ifelse(SS$X..exposed.in.an.isolated.chain <= 20,ymax-(ymax/3),"NA")))
SS$Residue..protein<-as.numeric(SS$Residue..protein)
SSf<-merge(Sl, SS, by.x = "POS", by.y = "Residue..protein", all.x = T )

#bind
SF<-cbind(S,ES,SSf[,-1])
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]

#ins merged
insS<-ins[ins$Protein=="nsp1",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)

#plot
PNSP1<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=0.9, col="red")+
  geom_point(aes(y = Freq_ins), size=0.9, col="blue")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(10, len, by = 10),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000),
                     labels = c(1, 10, 100, 1000, 10000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+
  annotate("segment", 
           x = 1, xend = len, 
           y = ymax-(ymax*0.49), yend = ymax-(ymax*0.49),
           size = 2,
           col="white"
           #arrow=arrow(length = unit(3,'mm'))
  )+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = ymax-(ymax*0.49),
           yend = ymax-(ymax*0.49),
           lty=c(1,5,1),
           colour = c("lightblue3","royalblue1",
                     "palevioletred3"), 
           size = c(5,2,5)
  )+
  annotate("text", 
           x = (annoS$xs+((annoS$xe-annoS$xs)/2)),
           hjust = 0.5,
           y =ymax-(ymax*0.82),
           colour = c("lightblue4","royalblue3","palevioletred4") ,
           size = 7,
           label=c("N-terminal head domain","Linker","C-terminal plug domain"),
           fontface = c('bold',"bold","bold"),
  )+
  annotate("segment", 
           x = c(84,142),
           xend = c(84,142), 
           y = 1, yend = ymax,
           colour = "pink", 
           size=c(13,8),
           alpha=0.3)+
  annotate("text", 
           x = c(84,142),
           y= ymax-(ymax*0.97),
           label = c("NSP1-HVR1","NSP1-HVR2"), 
           colour = "black", 
           fontface="bold",
           size=5
  )

PNSP1
pdf("NSP1.pdf", width = 16, height = 5)
PNSP1
dev.off()

#NSP6
len<-290

annoS<-anno[anno$Protein=="nsp6",]
S0<-df[df$Protein=="nsp6",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )

ymax<-20000000
ES<- read.delim("NSP6E.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
#merge
SF<-cbind(S,ES)
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
#ins merged
insS<-ins[ins$Protein=="nsp6",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)


PNSP6<- ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=0.9, col="red")+
  geom_point(aes(y = Freq_ins), size=0.9, col="blue")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(30, len, by = 30),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                     labels = c(1, 10, 100, 1000, 10000, 100000,1000000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
annotate("segment", 
           x = 1, xend = len, y = 9000000, yend = 9000000,
           colour ="black", size =1)+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = c(rep(9000000, 17)), 
           yend = c(rep(9000000, 17)),
           colour = c(rep("deeppink",5),
                      rep("blue",4),
                      rep("darkorange2",8)), 
           size =c(rep(6,5),
                   rep(6,4),
                   rep(8,8))
  )+
  annotate("segment", 
           x = 106,
           xend = 106, 
           y = 1, yend = ymax,
           colour = "pink", 
           size=22,
           alpha=0.3)+
  annotate("text", 
           x = 106,
           y= (max(top$Freq)+(max(top$Freq)/1.1)),
           label = c("NSP6-HVR"), 
           colour = "black", 
           fontface="bold",
           size=6)

PNSP6
pdf("NSP6.pdf", width = 16, height = 5)
PNSP6
dev.off()

##Fig2
library("gridExtra")
pdf("Fig2.pdf", width = 16, height = 15)
grid.arrange(PNSP1,PNSP6,PSpike ,                       # First row with one plot spaning over 2 columns
             ncol = 1, # Second row with 2 plots in 2 different columns
             nrow = 3)  
dev.off()


#ORF3a
len<-275
annoS<-anno[anno$Protein=="3a",]
S0<-df[df$Protein=="3a",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
ymax<-max(S[!(is.na(S))])*100
ES<- read.delim("ORF3aE.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
#bind
SF<-cbind(S,ES)
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
#ins merged
insS<-ins[ins$Protein=="3a",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)

PORF3a<- ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.5, col="red")+
  geom_point(aes(y = Freq_ins), size=1.5, col="blue")+
  
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(20, len, by = 20),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000,10000),
                     labels = c(1, 10, 100, 1000,10000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.1,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  annotate("segment", 
           x = 1, xend = len, y = 90000, yend = 90000,
           colour ="black", size =1)+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = c(rep(90000, 7)), 
           yend = c(rep(90000, 7)),
           colour = c(rep("deeppink",2),
                      rep("blue",2),
                      rep("darkorange2",3)), 
           size =c(rep(6,2),
                   rep(6,2),
                   rep(8,3))
  )+
  annotate("segment", 
           x = c(104,258),
           xend = c(103,258), 
           y = 1, yend = ymax,
           colour = "pink", 
           size=c(20,13),
           alpha=0.3)+
  annotate("text", 
           x = c(104,258),
           y= (max(top$Freq)+(max(top$Freq)/1.1)),
           label = c("ORF3a-HVR1","ORF3a-HVR2"), 
           colour = "black", 
           fontface="bold",
           size=5.5)

PORF3a
pdf("ORF3a.pdf", width = 16, height = 5)
PORF3a
dev.off()

#ORF8
len<-121
annoS<-anno[anno$Protein=="8b",]
S0<-df[df$Protein=="8b",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
ymax<-max(S[!(is.na(S))])*100

ES<- read.delim("ORF8E.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
SF<-cbind(S,ES)

mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
#ins merged
insS<-ins[ins$Protein=="8b",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)

PORF8<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.1, col="red")+
  geom_point(aes(y = Freq_ins), size=1.1, col="blue")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(10, len, by = 10)))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000,100000,1000000,10000000),
                     labels = c(1, 10, 100, 1000, 10000,100000,1000000,10000000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+
  annotate("segment", 
           x = 1, xend = len, 
           y = ymax-(ymax*0.49), yend = ymax-(ymax*0.49),
           size = 1,
           col="black"
           #arrow=arrow(length = unit(3,'mm'))
  )+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = ymax-(ymax*0.49),
           yend = ymax-(ymax*0.49),
           colour ="darkred", 
           size=3,
           arrow=arrow(length = unit(4,'mm'))
  )+
  
  annotate("text", 
           x = (annoS$xs+((annoS$xe-annoS$xs)/2)),
           hjust = 0.5,
           y =ymax-(ymax*0.9),
           colour = "darkred" ,
           size = 7,
           label="SP",
           fontface = "bold",
  )+
  annotate("segment", 
           x = c(66,120),
           xend = c(66,120), 
           y = 1, yend = ymax,
           colour = "pink", 
           size=c(12,12),
           alpha=0.3)+
  annotate("text", 
           x = c(66,120),
           y= ymax-(ymax*0.88),
           label = c("ORF8-HVR1","ORF8-HVR2"), 
           colour = "black", 
           fontface="bold",
           size=5.5)

PORF8
pdf("ORF8.pdf", width = 16, height = 5)
PORF8
dev.off()


#ORF7a
len<-121
annoS<-anno[anno$Protein=="7a",]
S0<-df[df$Protein=="7a",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
ymax<-max(S[!(is.na(S))])*100
ES<- read.delim("ORF7aE.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
SF<-cbind(S,ES)

mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]

#ins merged
insS<-ins[ins$Protein=="7a",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)


PORF7a<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.5, col="red")+
  geom_point(aes(y = Freq_ins), size=1.5, col="blue")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(10, len, by = 10)))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c(1, 10, 100, 1000, 10000, 100000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+
  annotate("segment", 
           x = 1, xend = len, 
           y = ymax-(ymax*0.49), yend = ymax-(ymax*0.49),
           size = 1,
           col="black"
           #arrow=arrow(length = unit(3,'mm'))
  )+
  annotate("segment", 
           x = annoS$xs, xend = annoS$xe, 
           y = c(ymax-(ymax*0.49),ymax-(ymax*0.49),ymax-(ymax*0.97),ymax-(ymax*0.49)), 
           yend = c(ymax-(ymax*0.49),ymax-(ymax*0.49),ymax-(ymax*0.97),ymax-(ymax*0.49)),
           colour =c("darkred","mediumpurple1",
                     "darkorange2","maroon"),
           size=c(1,2,2,2)
  )+
    annotate("segment", 
             x = 1, xend = 15, y = ymax-(ymax*0.49), yend = ymax-(ymax*0.49),
             colour ="darkred",
             arrow=arrow(length = unit(4,'mm'))
             )+
  annotate("text", 
           x = (annoS$xs+((annoS$xe-annoS$xs)/2)),
           hjust = 0.5,
           y =c(ymax-(ymax*0.8),ymax-(ymax*0.8),ymax-(ymax*0.99),ymax-(ymax*0.8)),
           colour = c("darkred","mediumpurple2",
                      "darkorange2","maroon") ,
           size = c(6,6,6,0),
           label=annoS$Annotation,
           fontface = c("bold","bold","bold.italic","bold"))+
  annotate("text", 
           x = 112,
           hjust = 0.5,y =ymax-(ymax*0.8),
           colour = "maroon" ,
           size = 5, label="ER retention motif",
           fontface = "bold")+
  annotate("segment", 
           x = 65,
           xend = 65, 
           y = 1, yend = ymax,
           colour = "pink", 
           size=85,
           alpha=0.3)+
  annotate("text", 
           x = c(14,60),
           y= ymax-(ymax*0.98),
           label = c("","ORF7a-HVR"), 
           colour = "black", 
           fontface="bold",
           size=5.5)

PORF7a
pdf("ORF7a.pdf", width = 16, height = 5)
PORF7a
dev.off()

#NSP3 Fig
len<-1945
annoS<-anno[anno$Protein=="nsp3",]
S0<-df[df$Protein=="nsp3",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
ymax<-max(S[!(is.na(S))])*100

ES<- read.delim("NSP3E.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
SF<-cbind(S,ES)

mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]

#ins merged
insS<-ins[ins$Protein=="nsp3",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)

PNSP3<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.5, col="red")+
  geom_point(aes(y = Freq_ins), size=1.5, col="blue")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(200, len, by = 200),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000,10000,100000),
                     labels = c(1, 10, 100, 1000,10000,100000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  annotate("segment", 
           x = 1, xend = len, y = 3000000, yend = 3000000,
           size = 2, col="black"
  )+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = c(3000000, 400000,rep(3000000, 4),400000,
                 3000000, 3000000,3000000,100000,100000), 
           yend = c(3000000, 400000,rep(3000000, 4),400000,
                    3000000, 3000000,3000000,100000,100000),
           colour =c("maroon","darkolivegreen3","firebrick1","purple","mediumpurple1",
                     "darkslategray2","royalblue3","goldenrod1","magenta","darkseagreen3",
                     "darkorange2","darkorange2"), 
           size =3)+
  annotate("text", 
           x = (annoS$xs+((annoS$xe-annoS$xs)/2)),
           hjust = c(rep(0.5,11),0),
           y = c(1200000, 200000,rep(1200000, 4),200000,
                 1200000, 1200000,500000,55000,55000
                 ),
           colour = c("maroon","darkolivegreen4","firebrick1","purple","mediumpurple1",
                      "darkslategray3","royalblue3","goldenrod1","magenta","darkseagreen3",
                      "darkorange2","darkorange2"),
           size = 6,
           label=annoS$Annotation,
           fontface = c(rep("bold",10),rep("bold.italic",2))
  )+
  annotate("segment", 
           x = c(200,400,1263),
           xend = c(200,400,1263), 
           y = 1, yend = ymax,
           colour = c("white", "white","pink"), 
           size=15,
           alpha=c(0,0,0.3))+
  annotate("text", 
           x = c(200,400,1263),
           y= 10000,
           label = c("","","NSP3-HVR"), 
           colour = "black", 
           fontface="bold",
           size=5.5)

PNSP3
pdf("NSP3.pdf", width = 16, height = 5)
PNSP3
dev.off()

#N
len<-419
annoS<-anno[anno$Protein=="N",]
S0<-df[df$Protein=="n",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
ymax<-max(S[!(is.na(S))])*100

ES<- read.delim("NE.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
SF<-cbind(S,ES)

#ins merged
insS<-ins[ins$Protein=="n",]
insS<-aggregate(Freq_ins~POS.AA,insS,sum)
SF<-merge(SF,insS, by.x="POS", by.y="POS.AA", all.x=T)

PN<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.5, col="red")+
  geom_point(aes(y = Freq_ins), size=1.5, col="blue")+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(50, len, by = 50),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000),
                     labels = c(1, 10, 100, 1000, 10000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+
  annotate("segment", 
           x = 1, xend = len, 
           y = 90000, yend = 90000,
           size = 2,
           col="white"
           #arrow=arrow(length = unit(3,'mm'))
  )+
  annotate("segment", 
           x = annoS$xs, 
           xend = annoS$xe, 
           y = c(900000, 900000, 200000,200000),
           yend = c(900000, 900000, 200000,200000),
           colour = c("lightblue3","palevioletred3", "red", "orange"), 
           size = 5
  )+
  annotate("text", 
           x = (annoS$xs+((annoS$xe-annoS$xs)/2)),
           hjust = 0.5,
           y =c(450000, 450000, 100000, 100000),
           colour = c("lightblue4","palevioletred4", "red", "orange") ,
           size = 7,
           label=annoS$Annotation,
           fontface = "bold",
  )+
  annotate("segment", 
           x = c(31,208),
           xend = c(31,208), 
           y = 1, yend = ymax,
           colour = "pink", 
           size=13,
           alpha=0.3)+
  annotate("text", 
           x = c(31,208),
           y= (max(top$Freq)+(max(top$Freq)/1.1)),
           label = c("N-HVR1","N-HVR2"), 
           colour = "black", 
           fontface="bold",
           size=5.5
  )

PN
pdf("N.pdf", width = 16, height = 5)
PN
dev.off()



##FigS2 cont
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
pdf("FigS2.c.pdf", width = 16, height = 25)
grid.arrange(
             PNSP3,
             PN,
             PORF3a,
             PORF7a,
             PORF8,
             ncol = 1, # Second row with 2 plots in 2 different columns
             nrow = 5)  
dev.off()

#NSP2
ymax<-200000
len<-638
ES<- read.delim("NSP2E.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
annoS<-anno[anno$Protein=="nsp2",]
S0<-df[df$Protein=="nsp2",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
SF<-cbind(S,ES)
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
PNSP2<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.5)+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(100, len, by = 100),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000),
                     labels = c(1, 10, 100, 1000, 10000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+
  annotate("segment", 
           x = 1, xend = len, 
           y = 90000, yend = 90000,
           size = 1,
           col="black"
  )+
  annotate("segment", 
           x = 268,
           xend =268, 
           y = 1, yend = ymax,
           colour = "pink", 
           size=10,
           alpha=0.3)+
  annotate("text", 
           x = 268,
           y= (max(top$Freq)+(max(top$Freq)/1.1)),
           label = "NSP2-HVR", 
           colour = "red", 
           fontface="bold",
           size=5.5)

PNSP2
pdf("NSP2.pdf", width = 16, height = 5)
PNSP2
dev.off()

#ORF6
ymax<-200000
len<-61
ES<- read.delim("ORF6E.txt", sep = "\t", header=TRUE) #COPY DATA DA proteins both up and down
ES$E<-as.numeric((ifelse(ES$Assignment=="E",ymax,"NA")))
annoS<-anno[anno$Protein=="6",]
S0<-df[df$Protein=="6",]
S1<-aggregate(Freq~POS.AA, S0, sum)
S2<-S1[S1$Freq >1,]
Sl<-data.frame(POS=c(1:len))
S<-merge(Sl, S2, by.x = "POS", by.y = "POS.AA", all.x = T )
SF<-cbind(S,ES)
mt <- SF[!(is.na(SF$Freq)), ]
top <-mt[rev(order(mt$Freq)),]
PORF6<-ggplot(SF, aes(x = POS)) +
  geom_point(aes(y = Freq), size=1.5)+
  geom_line(aes(y=E), size=2, col="seagreen3")+
  theme_classic()+
  scale_x_continuous(breaks = c(1,seq(10, len, by = 10),len))+
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1,ymax),
                     breaks = c(1, 10, 100, 1000, 10000),
                     labels = c(1, 10, 100, 1000, 10000)
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     #labels = trans_format("log10", math_format(10^.x))
  )+
  theme(axis.line.y = element_blank(), axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Number of genomes\n")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 22, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 22, face = "bold"))+ 
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 22, face = "plain"))+
  annotate("segment", 
           x = 1, xend = len, 
           y = 90000, yend = 90000,
           size = 1,
           col="black"
           #arrow=arrow(length = unit(3,'mm'))
  )+
  annotate("segment", 
           x = c(2,26),
           xend = c(2,26), 
           y = 1, yend = ymax,
           colour = "pink", 
           size=c(10,70),
           alpha=0.3)+
  annotate("text", 
           x = c(2,26),
           y= (max(top$Freq)+(max(top$Freq)/1.1)),
           label = c("ORF6-HVR1","ORF6-HS"), 
           colour = "red", 
           fontface="bold",
           size=5.5)

PORF6
pdf("ORF6.pdf", width = 16, height = 5)
PORF6
dev.off()
