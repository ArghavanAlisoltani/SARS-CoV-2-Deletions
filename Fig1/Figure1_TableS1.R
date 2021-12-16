# load pacakges
library(ggtree)
library(cowplot)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)
library("readxl")
library(viridis)
library("ggsci")
library(data.tree)
library(ggstance)
library(data.table)
library(ggtreeExtra)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(dplyr)

#####################################################
#code chunk1: count indels
#####################################################
setwd("~/Documents/Coronavirus/deletion_paper/Figure1_increased")
meta<-data.frame(fread("sars2_track_ns_ins.tsv", sep= "\t",header = T))
metanoN0<-meta[meta$number_of_Ns<300,]#300=29903*0.01%
rm(meta)
allp<-metanoN0%>% mutate(AA.Substitutions= strsplit(as.character(protein_variants.mutations.), ",")) %>% 
  unnest(AA.Substitutions)
allp$AA.Substitutions<-gsub("[ ]","",allp$AA.Substitutions)
allp<-allp[!(allp$AA.Substitutions %like% "fs"),]
agg<-aggregate(AA.Substitutions~gisaid_epi_isl,allp,toString, collapse = ",")
metanoN<-merge(metanoN0, agg, by="gisaid_epi_isl")
rm(agg,metanoN0)
#count indels tree
library(stringr)
metanoN$countdel<-str_count(metanoN$AA.Substitutions,coll("del)"))
table(metanoN$countdel)
metanoN$countdel<- ifelse(metanoN$countdel >0&metanoN$countdel <6, metanoN$countdel,0)
table(metanoN$countdel)
meta$countins<-str_count(meta$AA.substitutions..GISAID.,coll("ins"))
table(meta$countins)
my<-meta[meta$AA.substitutions..GISAID. %like% "Spike_ins",]
#meta$countins<-str_count(meta$protein.variants.mutations.,coll("_ins"))
fwrite(metanoN, "Out_meta_tree_IFDs.txt", sep="\t")
###############################################
#no_singletone del sustitution counts
##############################################
del<-allp[allp$AA.Substitutions %like% "del",]
sub<-allp[!(allp$AA.Substitutions %like% "del"),]
delc<-data.frame(table(del$AA.Substitutions))
subc<-data.frame(table(sub$AA.Substitutions))
delc$protein<-gsub("[_].*", "",delc$Var1)
subNS<-subc[subc$Freq >1,]
delNS<-delc[delc$Freq >1,]
proteincount<-data.frame(table(delNS$protein))
fwrite(proteincount,"proteincount_nosingleton.txt",sep= "\t")
fwrite(delNS,"deletionscount_nosingleton.txt",sep= "\t")
fwrite(subNS,"subscount_nosingleton.txt",sep= "\t")
rm(list = ls())
###############################################
#binomial proteins
##############################################  
count<-data.frame(fread("proteincount_nosingleton.txt", sep= "\t",header = T))
length<-data.frame(fread("Protein_length.txt", sep= "\t",header = T))
df<-merge(length,count, by.x = "Protein",by.y="Var1")
df$lengthall<-rep(sum(df$Length.of.Protein),25)
df$mutall<-rep(sum(df$Freq),25)

df$Odds.ratio_No_Singleton<-mapply (function (C, L, TC ,LP) {
  Oddsratio <- ((C/L)/(TC/LP))
  return (Oddsratio)
}, df$Freq, df$Length.of.Protein, df$mutall, df$lengthall)

df$pvalue_No_Singleton<- mapply (
  function (C, L, TC, LP) {
    binomR <- binom.test (C, TC, L/LP,
                          alternative = "two.sided")
    return (binomR$p.value)}
  ,df$Freq, df$Length.of.Protein, df$mutall, df$lengthall)

df$qvalue_No_Singleton<-p.adjust(df$pvalue_No_Singleton, method = "fdr")

#ORF1ab as background
nsps<-df[df$Protein %like% "NSP",]
df$nspsL<-rep(sum(nsps$Length.of.Protein),25)
df$nspsFreq<-rep(sum(nsps$Freq),25)

df$Odds.ratio_NSPs<-mapply (function (C, L, TC ,LP) {
  Oddsratio <- ((C/L)/(TC/LP))
  return (Oddsratio)
}, df$Freq, df$Length.of.Protein, df$nspsFreq, df$nspsL)

df$pvalue_NSPs<- mapply (
  function (C, L, TC, LP) {
    binomR <- binom.test (C, TC, L/LP,
                          alternative = "two.sided")
    return (binomR$p.value)}
  ,df$Freq, df$Length.of.Protein, df$nspsFreq, df$nspsL)

df$qvalue_NSPs<-p.adjust(df$pvalue_NSPs, method = "fdr")

fwrite(df[,c(2,3,4,7,9,12,14)], "binome_results.txt",sep= "\t")
rm(list = ls())
#######################################################################################
#sequence of numbers (repeat events for each genomic position) + for homoplasy
#######################################################################################
##substitutions
df<-data.frame(fread("subscount_nosingleton.txt", 
                     sep= "\t",header = T))
df$sub<-gsub(".*[(]","",df$Var1)
df$start<-as.numeric(gsub("[:].*","",df$sub))

fwrite(df, "subtitutions_pos.txt",row.names = T, 
       sep = "\t")
#dels
df<-data.frame(fread("deletionscount_nosingleton.txt", 
                     sep= "\t",header = T))
df$sub<-gsub(".*[:]","",df$Var1)
df$start<-as.numeric(gsub("[-].*","",df$sub))
df$end<-gsub(".*[-]","",df$sub)
df$end<-as.numeric(gsub("[del)]","",df$end))
df$dellen<-1+(df$end-df$start)
df$dellenP<-df$dellen/3
df$sub2<-gsub("[:].*","",df$Var1)
df$pstart<-as.numeric(gsub(".*[_]","",df$sub2))
df$pend<-df$pstart+df$dellenP
df$name<-paste(df$Var1,"|",df$Freq)
row.names(df)<-df$name
##genomic
dft<-t(df[,c(5,6)])
library(tidyr)
library(dplyr)
library(reshape2)
m<-seq.int(200,300)
class(m)
list1<-list()
for (i in colnames(dft)){
  list1[[i]]<-seq.int(dft[1,i],dft[2,i])
}
bind<-do.call(rbind, lapply(list1, data.frame, stringsAsFactors=FALSE))
bind$name<-gsub("[|].*","",row.names(bind))
bind$freq<-gsub(".*[|]","",row.names(bind))
bind$freq<-gsub("[.].*","",bind$freq)
colnames(bind)[1]<-"POS"
##protein
dft1<-t(df[,c(10,11)])
list1<-list()
for (i in colnames(dft1)){
  list1[[i]]<-seq.int(dft1[1,i],dft1[2,i])
}
bind1<-do.call(rbind, lapply(list1, data.frame, stringsAsFactors=FALSE))
bind1$name<-gsub("","",row.names(bind1))
colnames(bind1)[1]<-"POS"
bind1$name<-gsub("[|].*","",row.names(bind1))
bind1$Freq<-gsub(".*[|]","",row.names(bind1))
bind1$Freq<-gsub("[.].*","",bind1$Freq)
bind1$protein<-gsub("[_].*","",bind1$name)
fwrite(bind1, "countdel_per_protein_pos.txt",row.names = T, sep="\t")

#############################
#all plots and Fig1
#############################
#####################################################
#code chunk1:  panelA
#####################################################
# meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
# tree<-read.newick('nextstrain.nwk')
# metadata<-fread("hcov_global.tsv", sep= "\t",header = T)
# merged<-merge(metadata, meta,by.x = "gisaid_epi_isl", by.y = "gisaid_epi_isl",
#               all.x=T )# meta is obtained based on Lukasz's data
# tips<-data.frame(tree[["tip.label"]])
# merged2<-merge(tips, merged,  by.x ="tree...tip.label...", by.y="strain",
#              all.x=T )
# fwrite(merged2, "merged2.txt",sep= "\t")
merged2<-fread("merged2.txt", sep= "\t",header = T)
tree<-read.newick('nextstrain.nwk')

table(merged2$countdel)
#tree IFDs
p<-ggtree(tree,
          size=0.1, 
          mrsd="2021-11-17",
          #layout="circular",
          as.Date=T,
          #xlim = c(0,1.5),
          #aes(color=as.character(newpang)))
          aes(color="grey")) %<+% merged2[,] +
  geom_tippoint(aes(color=as.character(countdel),
                    size=as.character(countdel),
                    alpha=as.character(countdel),
                    shape=as.character(countdel)),
                #alpha=1,
                #size=0.01,
                show.legend = F
  ) +
  theme_tree2()+
  scale_x_date(date_labels = "%b-%y",date_breaks = "2 month")+
  
  scale_colour_manual(na.translate = T, 
                      na.value = "grey",
                      name=" ",
                      values=c("#B3D5E9" ,"#FFD1CB","#FEA298","#FE7464","#F81B02","#BA1402","grey"),
                      labels = c('No IFD','1 IFD', '2 IFD', '3 IFD','4 IFD','5 IFD','NA'))+  
  scale_shape_manual(na.value = 0,
                     values = c("0"=19,
                                "1"=19,
                                "2"=19,
                                "3"=19,
                                "4"=19,
                                "5"=19
                     ))+
  scale_size_manual(na.value = 0,
                    values = c(
                      "0"=0.3,
                      "1"=0.3,
                      "2"=0.4,
                      "3"=0.45,
                      "4"=0.5,
                      "5"=0.55
                    ))+
  scale_alpha_manual(na.value = 0,
                     values = c("0"=1,
                                "1"=1,
                                "2"=1,
                                "3"=1,
                                "4"=1,
                                "5"=1
                     ))+
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black", 
                                   size = 10,
                                   angle = 90))+
  
  theme(legend.position = c(0.15,0.75))+
  theme(legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9, angle = 0),
        #legend.key.height = unit(5, "lines")
        legend.key.size =unit(0.3, "cm"))+
  guides(fill = guide_legend(ncol = 1))+
  guides(fill = guide_legend(override.aes = list(shape = 22)),
         color = guide_legend(override.aes = list(size = 8)))+
  theme( plot.margin = unit( c(0.1,0.1,0.1,0.1) , units = "lines" ) )

p
######################################
#stackedbars indels
######################################
library(ggplot2)
library(tidyverse)
library(lubridate)
library(data.table)
meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
meta$Date1<-as.IDate(meta$date_collected)
meta$month<- round_date(meta$Date1, "2 months")
head(meta$month)
tbl<-data.frame(table(meta$month,meta$countdel))
tbl$time<-gsub("(.*)[-].*","\\1",tbl$Var1)
tbl$Var2<-ifelse(tbl$Var2 =="6","5",tbl$Var2)

p1<-ggplot(tbl, aes(fill=as.character(Var2), y=Freq, x=time)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity")+#"fill" gives percent and "stack" gives freq
  ylab("Genomes (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_fill_manual(values = c("#B3D5E9" ,"#FFD1CB","#FEA298","#FE7464","#F81B02","#BA1402","darkred"),
    labels = c('No IFD','1 IFD', '2 IFD', '3 IFD','4 IFD','5 IFD','6 IFD'))+  
  theme(legend.title = element_blank(), legend.text=element_text(size=8), 
        legend.key.size = unit(0.3,"cm"))+
  theme(legend.position = "top", legend.margin=margin(0,0,0,0))+
  theme(legend.direction="horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme( plot.margin = unit( c(0.1,0.1,1.1,0.1) , units = "lines" ))+
  theme(axis.title.x=element_blank())
  
p1
tbl1<-data.frame(table(meta$month,meta$countdel,meta$pangolin_lineage))
tbl1<-aggregate(countdel~pangolin_lineage+month,meta, median)
tbl1$time<-gsub("(.*)[-].*","\\1",tbl1$month)
tbl1$Freq<-rep(1, nrow(tbl1))
tbl1$countdel1<- ifelse(tbl1$countdel <1 , "0 IFD",
                           ifelse(tbl1$countdel >=1 &tbl1$countdel <2,"1 IFD",
                                  ifelse(tbl1$countdel >=2 &tbl1$countdel <3,"2 IFD",
                                         ifelse(tbl1$countdel >=3 &tbl1$countdel <4,"3 IFD",
                                                ifelse(tbl1$countdel >=4 &tbl1$countdel <5,"4 IFD",
                                                       ifelse(tbl1$countdel >=5& tbl1$countdel <10,"5 IFD", "0 IFD"
                                                       ))))))
                           
table(tbl1$countdel1)
p2<-ggplot(tbl1, aes(fill=countdel1, y=Freq, x=time)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity")+#"fill" gives percent and "stack" gives freq
  ylab("PANGO Lineages (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_fill_manual(values = 
                      c("#B3D5E9" ,"#FFD1CB","#FEA298","#FE7464","#F81B02","#BA1402","darkred"),
                    labels = c('No IFD','1 IFD', '2 IFD', '3 IFD','4 IFD','5 IFD','6 IFD'))+  
  theme(legend.title = element_blank(), legend.text=element_text(size=8), 
        legend.key.size = unit(0.3,"cm"))+
  theme(legend.position = "top", legend.margin=margin(0,0,0,0))+
  theme(legend.direction="horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme( plot.margin = unit( c(0.1,0.1,1.1,0.1) , units = "lines" ))+
  theme(axis.title.x=element_blank())
p2

###############################################################
#plot genome
###############################################################
#Plots for windows
bind<-fread("bind.txt",sep="\t")
df<-fread("subtitutions_pos.txt",sep="\t")
lab0<-29500
mylabels=(c(1:(lab0/100))*100)
my<-aggregate(as.numeric(freq)~POS,bind,sum)
my<-my[my$`as.numeric(freq)`>2,]
my2<-aggregate(as.numeric(Freq)~start,df,sum)
my2<-my2[my2$`as.numeric(Freq)`>2,]

a0<-data.frame(table(cut(as.numeric(my$POS),breaks=(lab0/100),
                           ordered = TRUE,labels = mylabels)))

a0$Var1<-(as.numeric(as.character(a0$Var1)))+200
a0$man<-ifelse(as.numeric(a0$Var1)%%3==0,a0$Var1," ")
a00<-c(1,0,1)
a01<-c(100,0," ")
aff<-c(29800,0," ")
aff1<-c(29900,0,29900)
a0<-rbind(a00,a01, a0,aff,aff1)
a0$Var1<-as.numeric(a0$Var1)
a1<-data.frame(table(cut(as.numeric(my2$start),breaks=(lab0/100),
                         ordered = TRUE,labels = mylabels)))
a1$Var1<-(as.numeric(as.character(a1$Var1)))+200
a00<-c(1,0)
a01<-c(100,0)
aff<-c(29800,0)
aff1<-c(29900,0)
a1<-rbind(a00,a01, a1,aff,aff1)
colnames(a0)<-c("Var1","dels","man")
colnames(a1)<-c("Var1","substit")

merged<-merge(a0,a1,by.x = "Var1",by.y="Var1")
merged$dels<-as.numeric(merged$dels)
#merged$man<-as.numeric(merged$man)

anno<-fread("Annotations_proteins.txt",sep="\t")
y=c(rep(c(140,110,125),9),125,140)
y1=c(rep(c(140,110,125),9),125,140)+5

##plot
p3<-ggplot(merged, aes(x = Var1)) +
  geom_line(col="red", aes(y = dels), group=1)+
  geom_line(col="blue", aes(y = substit), group=1)+
  #geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 10, face = "plain", angle = 0))+
  scale_x_continuous(breaks = c(1,5000,10000,15000,20000,25000,29000),
                   labels = c(1,5000,10000,15000,20000,25000,29000))+
  scale_y_continuous(labels = c("1","50","100"," "))+
  theme(axis.title.x=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank())+
  xlab(" ")+
  ylab("Distinct mutations (per 100 nt)")+
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 10, face = "plain"))+ 
  theme(axis.title.y = element_text(color = "grey20",hjust=0.3,
                                    size = 10, face = "bold"))+
  annotate("segment", 
           x = anno$S,
           xend = anno$E, 
           y = y,
           yend = y,
           colour = "black", 
           size=1)+
  annotate("text", 
           x = anno$S+((anno$E-anno$S)/2),
           y = y1, 
           label = anno$Label,
           colour = "black", 
           fontface="bold",
           size=2.2
  )

p3
######################################
##panels Fig1
######################################
lay <- rbind(c(1,2,NA),
             c(1,3,NA),
             c(4,4,NA))
fig1<-grid.arrange(p, p1, p2, p3,
             #grobs = ,
             layout_matrix = lay,
             widths = c(1.5, 2)
)
fig1
pdf("Fig1.pdf", width = 9, height = 9)
print(fig1)
dev.off()



