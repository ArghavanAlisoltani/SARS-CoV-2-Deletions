# load pacakges
library(ggtree)
library(cowplot)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)
library("readxl")
library(viridis)
library(ggplot2)
library("ggsci")
library(data.tree)
library(ggstance)
library(data.table)
library(ggtreeExtra)
library(gridExtra)
library(grid)
library(lattice)
library(dplyr)

#####################################################
#code chunk1: count indels
#####################################################
setwd("~/Downloads/msa_0106/Fig1")
meta<-data.frame(fread("sars2_track_ns_ins.tsv", sep= "\t",header = T))
metanoN0<-meta[meta$number_of_Ns<300,]#300=29903*0.01%
#metaGI<-data.frame(fread("metadata.tsv", sep= "\t",header = T))
#metanoN01<-merge(metaGI[,1:13], metanoN0, by.x="Accession.ID",
#                by.y="gisaid_epi_isl", all.y=T)
#fwrite(metanoN0,"metanoN0.tsv",sep= "\t")
rm(meta)
allp<-metanoN0%>% mutate(AA.Substitutions= strsplit(as.character(protein_variants.mutations.), ",")) %>% 
  unnest(AA.Substitutions)
rm(metanoN0)
allp$AA.Substitutions<-gsub("[ ]","", allp$AA.Substitutions)
allmut<-data.frame(table(allp$AA.Substitutions))
fwrite(allmut,"allmut.txt",sep= "\t")
rm(allmut)
#allp$AA.Substitutions1<-gsub("fs\\(",":_frmshft("
#                             ,allp$AA.Substitutions)
indels<-allp[allp$AA.Substitutions %like% "del"|
               allp$AA.Substitutions %like% "ins",]
fwrite(indels,"indels.tsv",sep= "\t")

rm(allp)
indels<-data.frame(fread("indels.tsv", sep= "\t",header = T))

agg<-aggregate(AA.Substitutions~gisaid_epi_isl,indels,toString, collapse = ",")

metaG<-data.frame(fread("metanoN0.tsv", sep= "\t",header = T))
metanoN<-merge(metaG, agg, by.x="Accession.ID",by.y="gisaid_epi_isl", all.x=T)
rm(agg,indels,metaG)
#count indels tree
library(stringr)
metanoN$countdel<-str_count(metanoN$AA.Substitutions,coll("del)"))
table(metanoN$countdel)
metanoN$countdel<- ifelse(metanoN$countdel >6,6,metanoN$countdel)
table(metanoN$countdel)
metanoN$countins<-str_count(metanoN$AA.Substitutions,coll("ins"))
table(metanoN$countins)
metanoN$countins<- ifelse(metanoN$countins >2,2,metanoN$countins)
table(metanoN$countins)
metanoN[is.na(metanoN)] <- 0

fwrite(metanoN, "Out_meta_tree_IFDs.txt", sep="\t")

#######################################################################################
#sequence of numbers (repeat events for each genomic position) 
#######################################################################################
rm(list = ls())
df<-data.frame(fread("allmut.txt", sep= "\t",header = T))
df<-df[df$Freq>1,]
df<-df[!(df$Var1 %like% "fs\\("),]
df$mut<-gsub(".*[(]","",df$Var1)
df$start<-as.numeric(gsub("[:].*","",df$mut))
# end of all mutations
df$delend<-gsub(".*[-]","",df$mut)
df$delend<-as.numeric(gsub("[del)]","",df$delend))

df$ins1<-gsub("[>].*","",df$mut)
df$ins1<-gsub(".*[:]","",df$ins1)
df$ins1<-nchar(df$ins1)
df$ins2<-gsub(".*[>]","",df$mut)
df$ins2<-gsub("[)]","",df$ins2)
df$ins2<-nchar(df$ins2)
df$end<-df$start+(df$ins2-df$ins1)
df$end<-ifelse(df$Var1 %like% "del", df$delend,df$end)

#protein
df$sub2<-gsub("[:].*","",df$Var1)
df$pstart<-as.numeric(gsub(".*[_]","",df$sub2))
df$len<-ifelse(df$Var1 %like% "ins",df$end-df$start ,(1+(df$end-df$start)))
df$lenP<-ifelse(df$Var1 %like% "ins"| df$Var1 %like% "del", df$len/3,df$len)
df$pend<-ifelse(df$Var1 %like% "ins",(df$pstart+df$lenP),
                ifelse(df$Var1 %like% "del", ((df$pstart+df$lenP)-1),
                       df$pstart)
)
df$start<-ifelse(df$Var1 %like% "ins",df$start+1,df$start)
df$pstart<-ifelse(df$Var1 %like% "ins",df$pstart+1,df$pstart)

df$name<-paste(df$Var1,"|",df$Freq)
row.names(df)<-df$name
df1<-df[,-c(3,5:7,9,11,12)]
fwrite(df1,"mutations_coordinates.tsv",sep="\t")
##genomic
dft<-t(df[,c("start","end")])
library(tidyr)
library(dplyr)
library(reshape2)
#m<-seq.int(200,300)#test
list1<-list()
for (i in colnames(dft)){
  list1[[i]]<-seq.int(dft[1,i],dft[2,i])
}
bind<-do.call(rbind, lapply(list1, data.frame, stringsAsFactors=FALSE))
bind$name<-gsub("[|].*","",row.names(bind))
bind$freq<-gsub(".*[|]","",row.names(bind))
bind$freq<-gsub("[.].*","",bind$freq)
colnames(bind)[1]<-"POS"
fwrite(bind,"muts_genomic.tsv")
##protein
dft1<-t(df[,c("pstart","pend")])
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
fwrite(bind1, "count_per_protein_pos.txt", row.names = T, sep="\t")

###############################################
#binomial proteins
##############################################  
rm(list = ls())
count<-data.frame(fread("allmut.txt", sep= "\t",header = T))
indels<-count[count$Var1 %like% "del"|
                count$Var1 %like% "ins",]
indels<-indels[!(indels$Var1 %like% "fs\\("),]
indels<-indels[indels$Freq>2,]#4976200*0.00001
indels<-indels[!duplicated(indels$Var1),]
indels$protein<-gsub("[_].*","",indels$Var1,)
proteins<-data.frame(table(indels$protein))
length<-data.frame(fread("Protein_length.txt", sep= "\t",header = T))
df<-merge(length,proteins, by.x = "Protein",by.y="Var1")
df$lengthall<-rep(sum(df$Length),nrow(df))
df$mutall<-rep(sum(df$Freq),nrow(df))

length<-data.frame(fread("Protein_length.txt", sep= "\t",header = T))
df<-merge(length,proteins, by.x = "Protein",by.y="Var1")
df$lengthall<-rep(sum(df$Length),25)
df$mutall<-rep(sum(df$Freq),25)

df$Odds.ratio_No_Singleton<-mapply (function (C, L, TC ,LP) {
  Oddsratio <- ((C/L)/(TC/LP))
  return (Oddsratio)
}, df$Freq, df$Length, df$mutall, df$lengthall)

df$pvalue_No_Singleton<- mapply (
  function (C, L, TC, LP) {
    binomR <- binom.test (C, TC, L/LP,
                          alternative = "two.sided")
    return (binomR$p.value)}
  ,df$Freq, df$Length, df$mutall, df$lengthall)

df$qvalue_No_Singleton<-p.adjust(df$pvalue_No_Singleton, method = "fdr")

#ORF1ab as background
nsps<-df[df$Protein %like% "NSP",]
df$nspsL<-rep(sum(nsps$Length),25)
df$nspsFreq<-rep(sum(nsps$Freq),25)

df$Odds.ratio_NSPs<-mapply (function (C, L, TC ,LP) {
  Oddsratio <- ((C/L)/(TC/LP))
  return (Oddsratio)
}, df$Freq, df$Length, df$nspsFreq, df$nspsL)

df$pvalue_NSPs<- mapply (
  function (C, L, TC, LP) {
    binomR <- binom.test (C, TC, L/LP,
                          alternative = "two.sided")
    return (binomR$p.value)}
  ,df$Freq, df$Length, df$nspsFreq, df$nspsL)

df$qvalue_NSPs<-p.adjust(df$pvalue_NSPs, method = "fdr")
df1<-df[,c(2,3,4,7,9,12,14)]
fwrite(df1, "binome_results.txt",sep= "\t")

#######################################################################################
#for homoplasy combinations 
#######################################################################################
rm(list = ls())
df0<-data.frame(fread("sars2_track_ns_ins.tsv", sep= "\t",header = T))
df<-data.frame(fread("muts_genomic.tsv", sep= ",",header = T))
dftop<-df[df$freq>4976200*0.0001,]
dfdel<-dftop[dftop$name %like% "del",]
dfins<-dftop[dftop$name %like% "ins",]
dfindels<-rbind(dfdel,dfins)
dfindels<-dfindels[!duplicated(dfindels$name),]
tree<-read.newick('global.tree')
tips<-data.frame(tree[["tip.label"]])
rm(tree)
merged2<-merge(tips, df0, by.x ="tree...tip.label...", 
               by.y= "gisaid_epi_isl", all.x=T)
indels<-merged2[merged2$protein_variants.mutations. %like% "del"|
                  merged2$protein_variants.mutations.  %like% "ins",]

rm(merged2, df0, tips)
allp0<-indels%>% 
  mutate(AA.Substitutions= 
           strsplit(as.character(indels$protein_variants.mutations.), ",")) %>% 
  unnest(AA.Substitutions)

rm(indels)
allp0$AA.Substitutions1<-gsub("[ ]","",allp0$AA.Substitutions)
allp<-allp0[allp0$AA.Substitutions1 %in% dfindels$name,]
head(allp[1:3,1:3])
my<-data.frame(table(allp$AA.Substitutions1))
fwrite(allp,"homoplasy_indels.tsv",sep="\t")
rm(list=ls())

allp<-fread("homoplasy_indels.tsv",sep="\t")
ff<-data.frame(table(allp$tree...tip.label...,allp$AA.Substitutions1))
head(ff[1:3,1:3])
rm(allp)
sprd<-spread(ff, key = "Var2", value = "Freq")
tree<-read.newick('global.tree')
tips<-data.frame(tree[["tip.label"]])
merged2<-merge(tips, sprd, by.x ="tree...tip.label...", 
               by.y= "Var1", all.x=T)
m<-merged2[,-1]
row.names(m)<-merged2$tree...tip.label...
rm(tree, tips, ff)
head(m[1:3,1:3])
m[is.na(m)] <- 0
head(m[1:3,1:3])
tm<-t(m)
tm<-data.frame(tm)
tm$snp<-row.names(tm)
head(tm)
fwrite(tm, "tm_homoplasy.tsv", sep = "\t")
rm(list = ls())
#merge coordinates
tm<-data.frame(fread("tm_homoplasy.tsv", sep= "\t",header = T))
df<-data.frame(fread("mutations_coordinates.tsv", sep= "\t"))
head(tm[1:3,1:3])
tbl<-data.frame(table(tm$snp))
df<-df[df$Var1 %in% tm$snp,]
my0<-merge(df,tm,by.x="Var1", by.y= "snp")
rm(tm)
#my<-my0[my0$Freq>10000,]
my<-my0[,8:ncol(my0)]
my[my > 1] <- 1
head(my[1:3,1:8])
my1<-cbind(my0[,3:4],my)
my2<-my0[,1:4]
head(my1[1:3,470121:470122])
head(my1[1:3,1:4])
m2<-data.frame(my1[,c(1,2,2683923)])
m1<-data.frame(my0[,c(3,4,2683928)])

fwrite(my1, "Indels_new.csv", sep = ",")
fwrite(my2, "Anno_homoplasy_new.csv", sep = ",")

#############################################
#homoplasy co-occ
############################################
hm<-fread("Indels_new.csv", sep=",")
Anno<-fread("Anno_homoplasy_new.csv", sep=",")
head(hm[1:3,1:4])
#sprdpango<-read.delim(pipe('pbpaste'),header=F,)
sprdpango<-fread("co-occ_combinations.txt", sep="\t", header = F)
Anno$new<-ifelse(Anno$Var1 %in% sprdpango$V2, 1, 0)# selected/manupulated previous one in excel
hm1<-hm[c(37,50,70,77,79,84,89,95,97,98,101,104),]
Anno1<-Anno[c(37,50,70,77,79,84,89,95,97,98,101,104),]
head(hm1[1:3,1:4])
rm(hm)
thm1<-t(hm1)
thm1<-data.frame(thm1)
head(thm1[1:3,])
thm1$X1_X2<-thm1$X1+thm1$X2
thm1$X1_X3<-thm1$X1+thm1$X3
thm1$X1_X4<-thm1$X1+thm1$X4
thm1$X1_X5<-thm1$X1+thm1$X5
thm1$X1_X6<-thm1$X1+thm1$X6
thm1$X1_X7<-thm1$X1+thm1$X7
thm1$X1_X8<-thm1$X1+thm1$X8
thm1$X1_X9<-thm1$X1+thm1$X9
thm1$X1_X10<-thm1$X1+thm1$X10
thm1$X1_X11<-thm1$X1+thm1$X11
thm1$X1_X12<-thm1$X1+thm1$X12

thm1$X2_X3<-thm1$X2+thm1$X3
thm1$X2_X4<-thm1$X2+thm1$X4
thm1$X2_X5<-thm1$X2+thm1$X5
thm1$X2_X6<-thm1$X2+thm1$X6
thm1$X2_X7<-thm1$X2+thm1$X7
thm1$X2_X8<-thm1$X2+thm1$X8
thm1$X2_X9<-thm1$X2+thm1$X9
thm1$X2_X10<-thm1$X2+thm1$X10
thm1$X2_X11<-thm1$X2+thm1$X11
thm1$X2_X12<-thm1$X2+thm1$X12

thm1$X3_X4<-thm1$X3+thm1$X4
thm1$X3_X5<-thm1$X3+thm1$X5
thm1$X3_X6<-thm1$X3+thm1$X6
thm1$X3_X7<-thm1$X3+thm1$X7
thm1$X3_X8<-thm1$X3+thm1$X8
thm1$X3_X9<-thm1$X3+thm1$X9
thm1$X3_X10<-thm1$X3+thm1$X10
thm1$X3_X11<-thm1$X3+thm1$X11
thm1$X3_X12<-thm1$X3+thm1$X12

thm1$X4_X5<-thm1$X4+thm1$X5
thm1$X4_X6<-thm1$X4+thm1$X6
thm1$X4_X7<-thm1$X4+thm1$X7
thm1$X4_X8<-thm1$X4+thm1$X8
thm1$X4_X9<-thm1$X4+thm1$X9
thm1$X4_X10<-thm1$X4+thm1$X10
thm1$X4_X11<-thm1$X4+thm1$X11
thm1$X4_X12<-thm1$X4+thm1$X12

thm1$X5_X6<-thm1$X5+thm1$X6
thm1$X5_X7<-thm1$X5+thm1$X7
thm1$X5_X8<-thm1$X5+thm1$X8
thm1$X5_X9<-thm1$X5+thm1$X9
thm1$X5_X10<-thm1$X5+thm1$X10
thm1$X5_X11<-thm1$X5+thm1$X11
thm1$X5_X12<-thm1$X5+thm1$X12

thm1$X6_X7<-thm1$X6+thm1$X7
thm1$X6_X8<-thm1$X6+thm1$X8
thm1$X6_X9<-thm1$X6+thm1$X9
thm1$X6_X10<-thm1$X6+thm1$X10
thm1$X6_X11<-thm1$X6+thm1$X11
thm1$X6_X12<-thm1$X6+thm1$X12

thm1$X7_X8<-thm1$X7+thm1$X8
thm1$X7_X9<-thm1$X7+thm1$X9
thm1$X7_X10<-thm1$X7+thm1$X10
thm1$X7_X11<-thm1$X7+thm1$X11
thm1$X7_X12<-thm1$X7+thm1$X12

thm1$X8_X9<-thm1$X8+thm1$X9
thm1$X8_X10<-thm1$X8+thm1$X10
thm1$X8_X11<-thm1$X8+thm1$X11
thm1$X8_X12<-thm1$X8+thm1$X12

thm1$X9_X10<-thm1$X9+thm1$X10
thm1$X9_X11<-thm1$X9+thm1$X11
thm1$X9_X12<-thm1$X9+thm1$X12

thm1$X10_X11<-thm1$X10+thm1$X11
thm1$X10_X12<-thm1$X10+thm1$X12

thm1$X11_X12<-thm1$X11+thm1$X12

head(thm1[1:3,])
thm2<-thm1
head(thm2)
thm2[thm2==1]<-0
thm2[thm2==2]<-1
head(thm2)
hmcomb<-t(thm2)
head(hmcomb[1:3,1:3])
hmcomb<-data.frame(hmcomb)
Anno_combi_homop<-hmcomb[,c(1,2)]#manupulated in excel
hmcomb<-data.frame(hmcomb)
hmcomb1<-hmcomb[13:78,-c(1,2)]
head(hmcomb1[1:3,1:3])
my0<-fread("Anno_combination_homoplasy.txt")
my<-cbind(my0[,10:11],hmcomb1)
head(my[1:3,1:3])
max(hmcomb1)
sum(thm2[,57])
fwrite(my, "Indels_combination.csv", sep = ",")

#############################
#all plots and Fig1
#############################
#####################################################
#code chunk1:  indels on tree
#####################################################
meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
table(meta$countdel)
tree<-read.newick('nextstrain.nwk')
metadata<-fread("hcov_global.tsv", sep= "\t",header = T)
merged<-merge(metadata, meta,by.x = "gisaid_epi_isl",
              by.y = "Accession.ID",all.x=T )# meta is obtained based on Lukasz's data
tips<-data.frame(tree[["tip.label"]])
merged2<-merge(tips, merged, by.x ="tree...tip.label...", by.y="strain",
               all.x=T )

table(merged2$countdel)

#tree IFDs
p0<-ggtree(tree,
           size=0.1, 
           mrsd="2021-12-27",
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
                      values=c("#B3D5E9" ,"#FFD1CB","#FEA298","#FE7464","#F81B02","#BA1402","darkred","grey"),
                      labels = c('No D','1 D', '2 D', '3 D','4 D','5 D',"6 D or more",'NA'))+  
  scale_shape_manual(na.value = 0,
                     values = c("0"=19,
                                "1"=19,
                                "2"=19,
                                "3"=19,
                                "4"=19,
                                "5"=19,
                                "6"=19
                     ))+
  scale_size_manual(na.value = 0,
                    values = c(
                      "0"=0.3,
                      "1"=0.3,
                      "2"=0.4,
                      "3"=0.45,
                      "4"=0.5,
                      "5"=0.55,
                      "6"=0.6
                    ))+
  scale_alpha_manual(na.value = 0,
                     values = c("0"=1,
                                "1"=1,
                                "2"=1,
                                "3"=1,
                                "4"=1,
                                "5"=1,
                                "6"=1
                     ))+
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black", 
                                   size = 10,
                                   angle = 90))+
  
  theme(legend.position = c(0.3,0.75))+
  theme(legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9, angle = 0),
        #legend.key.height = unit(5, "lines")
        legend.key.size =unit(0.1, "cm"))+
  guides(fill = guide_legend(ncol = 1))+
  guides(fill = guide_legend(override.aes = list(shape = 22)),
         color = guide_legend(override.aes = list(size = 6)))+
  theme( plot.margin = unit( c(0.1,0.1,0.1,0.1) , units = "lines" ) )
# geom_taxalink("SouthAfrica/NHLS-UCT-GS-A894/2021",
#               "hCoV-19/SouthAfrica/NHLS-UCT-GS-B895/2021", color="blue3", curvature=-.9) 
p0

#tree IFIs
table(merged2$countins)
pang<-data.frame(table(merged2$Nextstrain_clade))
pang<-c("21J (Delta)","21I (Delta)","20I (Alpha, V1)",	
        "21K (Omicron)","20J (Gamma, V3)", "21A (Delta)",
        "20H (Beta, V2)","21H (Mu)")
merged2$clade<-ifelse(merged2$Nextstrain_clade %in% pang,merged2$Nextstrain_clade,"Others" )
merged2$clade<-ifelse(merged2$clade %like% "Delta","Delta", merged2$clade)

pang1<-data.frame(table(merged2$clade))

p1<-ggtree(tree,
           size=0.1, 
           mrsd="2021-12-27",
           #layout="circular",
           as.Date=T,
           #show.legend = F,
           #xlim = c(0,1.5),
           aes(color=as.character(clade))) %<+% merged2[,] +
  #aes(color="grey")) %<+% merged2[,] +
  geom_tippoint(aes(color=as.character(countins),
                    size=as.character(countins),
                    alpha=as.character(countins),
                    shape=as.character(countins)),
                #alpha=1,
                #size=0.01,
                show.legend = F
  ) +
  theme_tree2()+
  scale_x_date(date_labels = "%b-%y",date_breaks = "2 month")+
  scale_colour_manual(na.translate = T, 
                      na.value = "grey",
                      name=" ",
                      values=c(
                        "#B3D5E9","#FEA298","red","orange","plum","blue",
                        "yellow","skyblue","lightgreen","grey"),
                      labels = c('No I','1 I','2 I or more',"Beta",
                                 "Alpha","Gamma", "Mu","Omicron","Delta","Others"))+  
  scale_shape_manual(na.value = 0,
                     values = c("0"=19,
                                "1"=19,
                                "2"=19,
                                "3"=19,
                                "4"=19,
                                "5"=19,
                                "6"=19
                     ))+
  scale_size_manual(na.value = 0,
                    values = c(
                      "0"=0.3,
                      "1"=0.6,
                      "2"=0.7,
                      "3"=0.6,
                      "4"=0.65,
                      "5"=0.55,
                      "6"=0.6
                    ))+
  scale_alpha_manual(na.value = 0,
                     values = c("0"=1,
                                "1"=1,
                                "2"=1,
                                "3"=1,
                                "4"=1,
                                "5"=1,
                                "6"=1
                     ))+
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black", 
                                   size = 10,
                                   angle = 90))+
  
  theme(legend.position = c(0.2,0.75))+
  theme(legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9, angle = 0),
        #legend.key.height = unit(2, "lines"),
        legend.key.size =unit(0.1, "cm"))+
  guides(fill = guide_legend(ncol = 2))+
  guides(colour = guide_legend(ncol = 2))+
  guides(fill = guide_legend(override.aes = list(shape = 22)),
         color = guide_legend(override.aes = list(size = 6)))+
  theme( plot.margin = unit( c(0.1,0.1,0.1,0.1) , units = "lines" ) )
# geom_taxalink("SouthAfrica/NHLS-UCT-GS-A894/2021",
#               "hCoV-19/SouthAfrica/NHLS-UCT-GS-B895/2021", color="blue3", curvature=-.9) 
p1
######################################
#stackedbars indels
######################################
library(ggplot2)
library(tidyverse)
library(lubridate)
library(data.table)
#meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
#meta<-meta[meta$Host=="Human",]
meta$Date1<-as.IDate(meta$Collection.date)
meta$month<- round_date(meta$Date1, "2 months")
head(meta$month)
tbl<-data.frame(table(meta$month,meta$countdel))
tbl$time<-gsub("(.*)[-].*","\\1",tbl$Var1)

p2<-ggplot(tbl, aes(fill=as.character(Var2), y=Freq, x=time)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity")+#"fill" gives percent and "stack" gives freq
  ylab("Genomes (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_fill_manual(values = c("#B3D5E9" ,"#FFD1CB","#FEA298","#FE7464","#F81B02","#BA1402","darkred"),
                    labels = c('No D','1 D', '2 D', '3 D','4 D','5 D','6 D or more'))+  
  theme(legend.title = element_blank(), legend.text=element_text(size=8), 
        legend.key.size = unit(0.3,"cm"))+
  theme(legend.position = "top", legend.margin=margin(0,0,0,0))+
  theme(legend.direction="horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme( plot.margin = unit( c(0.1,0.1,1.1,0.1) , units = "lines" ))+
  theme(axis.title.x=element_blank())

p2
tbl1<-data.frame(table(meta$month,meta$countdel,meta$Pango.lineage))
tbl1<-aggregate(countdel~Pango.lineage+month,meta, median)
tbl1$time<-gsub("(.*)[-].*","\\1",tbl1$month)
tbl1$Freq<-rep(1, nrow(tbl1))
tbl1$countdel1<- ifelse(tbl1$countdel <1 , "0 D",
                        ifelse(tbl1$countdel >=1 &tbl1$countdel <2,"1 D",
                               ifelse(tbl1$countdel >=2 &tbl1$countdel <3,"2 D",
                                      ifelse(tbl1$countdel >=3 &tbl1$countdel <4,"3 D",
                                             ifelse(tbl1$countdel >=4 &tbl1$countdel <5,"4 D",
                                                    ifelse(tbl1$countdel >=5 &tbl1$countdel <6,"5 D",
                                                           ifelse(tbl1$countdel >=6& tbl1$countdel <10,"6 D ", "0 D"
                                                           )))))))

table(tbl1$countdel1)
p3<-ggplot(tbl1, aes(fill=countdel1, y=Freq, x=time)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity")+#"fill" gives percent and "stack" gives freq
  ylab("PANGO Lineages (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_fill_manual(values = 
                      c("#B3D5E9" ,"#FFD1CB","#FEA298","#FE7464","#F81B02","#BA1402","darkred"),
                    labels = c('No D','1 D', '2 D', '3 D','4 D','5 D','6 D or more'))+  
  theme(legend.title = element_blank(), legend.text=element_text(size=8), 
        legend.key.size = unit(0.3,"cm"))+
  theme(legend.position = "top",
        legend.margin=margin(0,0,0,0))+
  theme(legend.direction="horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme( plot.margin = unit( c(0.1,0.1,1.1,0.1) , units = "lines" ))+
  theme(axis.title.x=element_blank())
p3

############################################################
#Insersions
######################################
tbl2<-data.frame(table(meta$month,meta$countins))
tbl2$time<-gsub("(.*)[-].*","\\1",tbl2$Var1)

p4<-ggplot(tbl2, aes(fill=as.character(Var2), y=Freq, x=time)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity")+#"fill" gives percent and "stack" gives freq
  ylab("Genomes (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_fill_manual(values = c("#B3D5E9","#FEA298","red"),
                    labels = c('No I','1 I', '2 I or more', '3 I','3 I','5 D','6 D'))+  
  theme(legend.title = element_blank(), legend.text=element_text(size=8), 
        legend.key.size = unit(0.3,"cm"))+
  theme(legend.position = "top", legend.margin=margin(0,0,0,0))+
  theme(legend.direction="horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme( plot.margin = unit( c(0.1,0.1,1.1,0.1) , units = "lines" ))+
  theme(axis.title.x=element_blank())

p4
tbl3<-data.frame(table(meta$month,meta$countins,meta$Pango.lineage))
tbl3<-aggregate(countins~Pango.lineage+month,meta, median)
tbl3$time<-gsub("(.*)[-].*","\\1",tbl3$month)
tbl3$Freq<-rep(1, nrow(tbl3))
tbl3$countins<- ifelse(tbl3$countins <1, "0 D",
                       ifelse(tbl3$countins >=1 &tbl3$countins <2,"1 I",
                              ifelse(tbl3$countins >=2 & tbl3$countins <10,"2 I", "0 D"
                              )))

table(tbl3$countins)
p5<-ggplot(tbl3, aes(fill=countins, y=Freq, x=time)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity")+#"fill" gives percent and "stack" gives freq
  ylab("PANGO Lineages (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_fill_manual(values = 
                      c("#B3D5E9","#FEA298","red"),
                    labels = c('No I','1 I', '2 I or more', '3 D','4 D','5 D','6 D'))+  
  theme(legend.title = element_blank(), legend.text=element_text(size=8), 
        legend.key.size = unit(0.3,"cm"))+
  theme(legend.position = "top",
        legend.margin=margin(0,0,0,0))+
  theme(legend.direction="horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme( plot.margin = unit( c(0.1,0.1,1.1,0.1) , units = "lines" ))+
  theme(axis.title.x=element_blank())
p5


###############################################################
#plot genome
###############################################################
mut0<-fread("muts_genomic.tsv",sep = ",", header = T)
mut0<-mut0[mut0$freq >4976200*0.00007,]
POS<-c(1:29900)
window<-rep(c(1:299),each=100)*100
final<-data.frame(cbind(POS,window))
mut<-merge(final,mut0,by="POS",all.x = T)
mut$del<-ifelse(mut$name %like% "del"& mut$freq>400,1,0)
mut$del<-ifelse(mut$name %like% "7A"& mut$freq<1400,0,mut$del)

mut$ins<-ifelse(mut$name %like% "ins",1,0)
mut$sub<-ifelse(mut$name %like% "ins"|
                  mut$name %like% "del",0,1)
mut$sub<-as.numeric(ifelse(is.na(mut$name),"0",mut$sub))
mut1<-mut[,c(2,5,6,7)]
merged<-aggregate(.~window,mut1,sum)

merged$man<-ifelse(as.numeric(merged$window)%%3==0,
                   merged$window," ")
a00<-c(1,0,0,0,1)
#aff1<-c(30000,0,0,0,30000)
merged1<-rbind(a00,merged)

#merging
anno<-fread("Annotations_proteins.txt",sep="\t")
y=c(rep(c(130,100,115),9),115,130)
y1=c(rep(c(130,100,115),9),115,130)+5

##plot
p6<-ggplot(merged1, aes(x = window)) +
  geom_line(col="red" ,aes(y = del), group=1)+
  geom_line(col="darkgreen", aes(y = sub), group=1)+
  geom_line(col="blue", aes(y = ins), group=1)+
  #geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 10, face = "plain", angle = 0))+
  scale_x_continuous(breaks = c(1,5000,10000,15000,20000,25000,29000),
                     labels = c(1,5000,10000,15000,20000,25000,29000))+
  scale_y_continuous(breaks = c(1,20,40,60,80,90,100),
                     labels = c("1","20","40","60","80","",""))+
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

p6
######################################
##panels tree
######################################
lay <- rbind(c(1,2,3,3),
             c(1,2,4,4),
             # c(1,2,4,4),
             #c(1,2,4,4),
             c(5,5,5,5)
)
fig1<-grid.arrange(p0, p1,
                   p3, p5,
                   p6,
                   #grobs = ,
                   layout_matrix = lay
                   #widths = c(1, 1)
)
fig1

lay2 <- rbind(c(1,1),
              c(2,2))
figS1<-grid.arrange(
  p2, p4,
  #grobs = ,
  layout_matrix = lay2
  #widths = c(1, 1)
)
figS1

#################################
#each indels on tree
#################################
#####################################################
#code chunk1:  indels on tree
#####################################################
meta<-fread("sars2_track_ns_ins.tsv", sep= "\t",header = T)
tree<-read.newick('nextstrain.nwk')
metadata<-fread("hcov_global.tsv", sep= "\t",header = T)
merged<-merge(metadata, meta,by.x = "gisaid_epi_isl",
              by.y = "gisaid_epi_isl",all.x=T )# meta is obtained based on Lukasz's data
tips<-data.frame(tree[["tip.label"]])
merged2<-merge(tips, merged, by.x ="tree...tip.label...", by.y="strain",
               all.x=T )

table(merged2$countdel)

#tree IFIs
# table(merged2$countins)
# pang<-data.frame(table(merged2$Nextstrain_clade))
# pang<-c("21J (Delta)","21I (Delta)","20I (Alpha, V1)",
#         "21K (Omicron)","20J (Gamma, V3)", "21A (Delta)",
#         "20H (Beta, V2)","21H (Mu)")
# merged2$clade<-ifelse(merged2$Nextstrain_clade %in% pang,merged2$Nextstrain_clade,"Others" )
# merged2$clade<-ifelse(merged2$clade %like% "Delta","Delta", merged2$clade)
# 
# pang1<-data.frame(table(merged2$clade))
# fwrite(merged2, "ForEach_tree_data.tsv",sep="\t" )

merged2<-data.frame(fread("ForEach_tree_data.tsv", sep="\t"))
allp<-merged2%>% mutate(AA.Substitutions= strsplit(as.character(protein_variants.mutations.), ",")) %>% 
  unnest(AA.Substitutions)
allp$AA.Substitutions<-gsub("[ ]","", allp$AA.Substitutions)

df<-data.frame(fread("muts_genomic.tsv", sep= ",",header = T))
dftop<-df[df$freq>4976200*0.0001,]
dfdel<-dftop[dftop$name %like% "del",]
dfins<-dftop[dftop$name %like% "ins",]
dfindels<-rbind(dfdel,dfins)
dfindels<-dfindels[!duplicated(dfindels$name),]
allp<-allp[allp$AA.Substitutions %in% dfindels$name,]
m<-data.frame(table(allp$tree...tip.label...,allp$AA.Substitutions))
sprd<-spread(m, key = "Var2", value = "Freq")
#names(sprd)<-gsub("[()].*","",names(sprd))
#sprd<-sprd[col]
merged3<-merge(merged2, sprd, by.x ="tree...tip.label...", by.y="Var1",
               all.x=T )
#merged3[is.na(merged3)] <- 0
merged3$nsp1_8<-ifelse(merged3$protein_variants.mutations. %like% "NSP1_8",1,0)
merged3$co_occ<-merged3$`S_69:69-70del(21765:21765-21770del)`+
  merged3$`NSP6_106:106-108del(11288:11288-11296del)`
merged3$co_occ<-ifelse(merged3$co_occ==2,1,0)

table(merged3$nsp1_8)
table(merged3$co_occ)
tree<-read.newick('nextstrain.nwk')
library("foreach")
library("doParallel")
registerDoParallel(cores=4)
foreach(i=names(merged3)[102]) %dopar%{
  my_col <- sym(i)
p<-ggtree(tree,
           size=0.1, 
           mrsd="2021-12-27",
           #layout="circular",
           as.Date=T,
           #show.legend = F,
           #xlim = c(0,1.5),
           aes(color=as.character(clade))) %<+% merged3[,] +
  #aes(color="grey")) %<+% merged2[,] +
  geom_tippoint(aes(color=as.character(!!my_col),
                    size=as.character(!!my_col),
                    alpha=as.character(!!my_col),
                    shape=as.character(!!my_col)),
                #alpha=1,
                #size=0.01,
                show.legend = F
  ) +
  theme_tree2()+
  scale_x_date(date_labels = "%b-%y",date_breaks = "2 month")+
  scale_colour_manual(na.translate = T, 
                      na.value = "grey",
                      name=" ",
                      values=c(
                        "white","red","orange","plum","blue",
                        "yellow","skyblue","lightgreen","grey"),
                      labels = c(' ',gsub("[()].*","",paste(i)),"Beta",
                                 "Alpha","Gamma", "Mu","Omicron","Delta","Others"))+  
  scale_shape_manual(na.value = 0,
                     values = c("0"=19,
                                "1"=19
                     ))+
  scale_size_manual(na.value = 0,
                    values = c(
                      "0"=0,
                      "1"=0.5
                    ))+
  scale_alpha_manual(na.value = 0,
                     values = c("0"=0,
                                "1"=1
                     ))+
  theme(axis.text.x = element_text(face = "bold",
                                   color = "black", 
                                   size = 6,
                                   angle = 90))+
  #theme(legend.position = c(0.4,0.75))+
  theme(legend.position = "none")+
  theme(legend.title = element_text(size = 11), 
        legend.text = element_text(size = 9, angle = 0),
        #legend.key.height = unit(2, "lines"),
        legend.key.size =unit(0.1, "cm"))+
  guides(fill = guide_legend(ncol = 2))+
  guides(colour = guide_legend(ncol = 2))+
  guides(fill = guide_legend(override.aes = list(shape = 22)),
         color = guide_legend(override.aes = list(size = 6)))+
  theme( plot.margin = unit( c(0.1,0.1,0.1,0.1) , units = "lines" ) )+
  #ggtitle(gsub("[()].*","",paste(i)))+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.7))
ggsave(paste0("tree_",i,".pdf"),p,
       width = 4, height = 10,
       units = "cm",limitsize = FALSE)
}

#####################################################
#code chunk:  Heatmaps and combinations
#####################################################
library(ggplot2)
library(tidyverse)
library(lubridate)
library(data.table)
meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
meta<-meta[meta$Host=="Human",]
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.1.529"|meta$Pango.lineage %like% "BA.", "Omicron","Others")
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.1.7"|meta$Pango.lineage %like% "Q.", "Alpha",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.617.2"|meta$Pango.lineage %like% "AY.", "Delta",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.351", "Beta",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.6211", "Mu",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "P.1", "Gamma",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "C.37", "Lambda",meta$newpango)

table(meta$newpango)

#####################################################
#code chunk:  combinations
#####################################################
tblcomb<-data.frame(table(meta$AA.Substitutions))
tblcomb<-tblcomb[tblcomb$Freq>5000,]
tblcomb<-tblcomb[-1,]
meta2<-meta[meta$AA.Substitutions %in% tblcomb$Var1,]
tblcomb2<-data.frame(table(meta2$AA.Substitutions, meta2$newpango))
tblcomb2<-tblcomb2[tblcomb2$Freq>0,]
agg<-aggregate(Var2~Var1,tblcomb2,toString, sep=",")
merged<-merge(tblcomb,agg,by="Var1")
#meta$Date1<-as.IDate(meta$Collection.date)
#meta$month<- round_date(meta$Date1, "1 months")
meta$region<-sub("[ /].*","", meta$Location)
meta$timelab<-sub("-[^-]+$","", meta$Collection.date)

meta<-meta%>%
  filter(!meta$timelab %in% c("2020", "2021","2021-1"), )
tbl<-data.frame(table(meta$timelab))

allp<-meta%>% 
  mutate(AA.Substitutions= strsplit(as.character(AA.Substitutions), ",")) %>% 
  unnest(AA.Substitutions)

allp$AA.Substitutions<-gsub("[ ]","",
                            allp$AA.Substitutions)
nsp1<-allp[allp$AA.Substitutions %like% "NSP6_",]
tbl<-data.frame(table(nsp1$newpango,nsp1$AA.Substitutions))
################################################
#recurrence
###############################################
##recurrence time
time<-data.frame(table(allp$AA.Substitutions, allp$timelab))
sprdtime<-spread(time, key = Var2, value=Freq)
#sprdtime<-sprdtime[,-1]
timecount<-data.frame(sprdtime [,1],rowSums(sprdtime [,-1]>5))

##recurrence region
region<-data.frame(table(allp$AA.Substitutions,allp$region))
sprdregion<-spread(region, key = Var2, value=Freq)
#sprdregion<-sprdregion[-1,]
regioncount<-data.frame(sprdregion [,1],rowSums(sprdregion [,-1]>5))

##recurrence PANGO
pango<-data.frame(table(allp$AA.Substitutions,allp$Pango.lineage))
pangocl<-data.frame(table(allp$Pango.lineage))

sprdpango<-spread(pango, key = Var2, value=Freq)
#sprdpango<-sprdpango[-1,]
pangocount<-data.frame(sprdpango [,1],rowSums(sprdpango [,-1]>5))

##recurrence GISAID clades
gisad<-data.frame(table(allp$AA.Substitutions,allp$Clade))
gisadcl<-data.frame(table(allp$Clade))
sprdgisad<-spread(gisad, key = Var2, value=Freq)
#sprdgisad<-sprdgisad[-1,]
gisadcount<-data.frame(sprdgisad [,1],rowSums(sprdgisad [,-1]>5))

Freqall<-data.frame(table(allp$AA.Substitutions))
reccurrent<-cbind(Freqall, pangocount, gisadcount, regioncount, timecount)
fwrite(reccurrent, "reccurrent.tsv",sep="\t")
#####merge homoplasy and recurrent results
df1<-fread("Anno_homoplasy_new.csv", sep= ",")
df2<-fread("consistencyIndexReport_11-01-22.txt", sep= "\t")
df3<-fread("reccurrent.tsv", sep= "\t")
homo<-cbind(df1,df2)
merged<-merge(homo,df3, by="Var1")
fwrite(merged, "Table_homoplasy.tsv",sep="\t")

##########################################
#co-occure
##########################################
#tbl$time<-gsub("(.*)[-].*","\\1",tbl$Var1)
counts<-data.frame(table(allp$AA.Substitutions))
top<-counts[counts$Freq>1000,]
top<-top[-1,]
allptop<-allp[allp$AA.Substitutions %in% top$Var1,]
fwrite(allptop, "for_heatmap.tsv", sep="\t")
#re-import
allptop<-fread("for_heatmap.tsv",sep="\t")
#allptop<-allptop[allptop$newpang %like% "Omicron",]
pango<-data.frame(table(allptop$AA.Substitutions, allptop$Accession.ID))
sprdpango<-spread(pango, key = Var2, value=Freq)
head(sprdpango[1:3,1:3])
row.names(sprdpango)<-sprdpango$Var1
sprdpango<-sprdpango[,-1]
fwrite(sprdpango,"for_co-occure.tsv",sep="\t", row.names = T)
rm(sprdpango,allptop,pango)
sprdpango<-data.frame(fread("for_co-occure.tsv",sep="\t"),row.names=1)
#sprdpango<-read.delim(pipe('pbpaste'),header=TRUE,row.names = 1)
#head(sprdpango[1:3,1:3])
#ttest<-data.frame(t(sprdpango))
#cor2<-cor(ttest, method="pearson")
#cortest<-cor.test()

mx<-data.matrix(sprdpango)
head(mx[1:3,1:3])
mxt<-t(mx)
library(Hmisc)
pval<-rcorr(mxt,type="spearman")
p<-round(pval[["P"]],4)
cor<-round(pval[["r"]],4)

#library(ltm)
#corb<-biserial.cor(ttest$mu1,ttest$mu2)
#cor<-cor(mxt)
#cortest<-cor.test(mxt)

co_oc<-crossprod(mxt)
write.table(co_oc, "co_occure.tsv", sep="\t")
write.table(cor, "cor.tsv", sep="\t")
write.table(p, "pval_cor.tsv", sep="\t")

2016################################
#ggcooccur heatmap
################################
library(ggplot2)
library(reshape2)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
met<-  read.delim("co_occure_final_input_ggheatmap.txt",
                  sep = "\t", header=TRUE,
                  check.names = F, row.names = 1) #COPY DATA
cor<-  read.delim("cor_final_for_heat.txt",
                  sep = "\t", header=TRUE,
                  check.names = F, row.names = 1) #COPY DATA
V<-data.matrix(met)
V<-round(data.matrix(cor), 2)

cormat <- reorder_cormat(V)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#melted_cormat[,order(colnames(melted_cormat))]

#elted_cormat <- melt(df1, na.rm = TRUE)

# Create a ggheatmap
library(RColorBrewer)
cols<-rev(brewer.pal(11,"Spectral"))
#cols<-rev(heat.colors(16))
#cols <- cm.colors(20)
#cols<- rev(topo.colors(20))
max(cor)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  # scale_fill_gradient2(low = "lavender", high = "red", mid = "orange", 
  #                     midpoint = 1, limit = c(0,6000),
  #                     breaks=seq(0,1,100,1000,5000,200000), 
  #                     space = "Lab", 
  #                      name="Co-occurrence") +
  scale_fill_gradientn(colors = cols[-c(1,11)],
                       #colors = c("lavender","pink","orange"),
                       limits=c(-1,1),
                       na.value = "indianred1",
                       #breaks=seq(0,200000, length.out=100)
  ) +
  #scale_y_discrete(expand=c(0,500))+
  #scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"),
  #   na.value = "red")+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 40, vjust = 1, face = "bold",
                                   size = 8, hjust = 1
  ))+
  theme(axis.text.y = element_text(angle = 0,hjust = 1,  face = "bold", 
                                   size = 8))


p1<-ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", 
            size = 2.5, angle = 0,) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    legend.position = "right",
    legend.direction = "vertical"
    #legend.title =element_text("Genomes") 
  )+
  scale_y_discrete(position = "right")+
  labs(fill = "Correlation")+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.9), "cm"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 6,
                               title.position = "top", title.hjust = 0.5))

# Print the heatmap
print(p1)

###############################
#heatmap fig4
###############################
allptop<-fread("for_heatmap.tsv",sep="\t")

tbl0<-data.frame(table(allptop$newpango,allptop$AA.Substitutions))
pango<-allptop[!duplicated(allptop$Accession.ID),]
tbl1<-data.frame(table(pango$newpango))
merged<-merge(tbl0,tbl1,by="Var1",all.x = T)
merged$average<-(merged$Freq.x/merged$Freq.y)*100
merged<-merged[,-c(3,4)]
sprd<-spread(merged, key = Var1, value=average)
sprd<-sprd[,-8]
sprd$sum<-rowSums(sprd[,-1])
sprd<-sprd[!(sprd$Var2 %like% "fs"),]
sprd<-sprd[order(sprd$sum),]
sprd<-sprd[-(1:45),-8]

mat<- data.matrix(sprd[,2:ncol(sprd)]) #labeling numerical data as matrix
rnames <- sprd[,1] #where the name for each row is
rnames<- gsub("[(].*","",rnames)
#rnames<- sub("[_].*.[:]","_",rnames)
rnames<- sub("[:].*.[-]","-",rnames)
rownames(mat) <-rnames
mat<-t(mat)
#complex heatmaps
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)
library(data.table)
Heatmap(mat[,1:ncol(mat)],
        name = "IFD (%Genomes)",
        row_km =1, 
        row_names_gp = grid::gpar(fontsize = 11),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_side = "right",
        width = unit(6, "cm"),
        col =colorRamp2(c(0,1,2,4,8,16,32,64,80),
                        c("#FCFDBFFF",
                          #"#FA815FFF",
                          "#F4685CFF",
                          "#E85362FF",
                          "#D6456CFF",
                          "#C03A76FF",
                          "#AB337CFF",
                          "#952C80FF",
                          "#400F73FF",
                          "#000004FF"
                        )),
        show_row_names = T,
        show_column_names = T, 
        cluster_rows = T,
        cluster_columns = T,
        #row_dend_gp = 0.5,
        )
theme( plot.margin = unit( c(0.1,0.1,0.1,0.1) , units = "lines" ) )
#margin(t, r, l, b)

 ###pango percentage mutations
meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.1.529"|meta$Pango.lineage %like% "BA.", "Omicron","Others")
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.1.7"|meta$Pango.lineage %like% "Q.", "Alpha",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.617.2"|meta$Pango.lineage %like% "AY.", "Delta",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.351", "Beta",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.6211", "Mu",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "P.1", "Gamma",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "C.37", "Lambda",meta$newpango)

table(meta$newpango)

allp<-meta%>% 
  mutate(AA.Substitutions= strsplit(as.character(AA.Substitutions), ",")) %>% 
  unnest(AA.Substitutions)

allp$AA.Substitutions<-gsub("[ ]","",
                            allp$AA.Substitutions)
counts<-data.frame(table(allp$AA.Substitutions))
top<-counts[counts$Freq>300,]
top<-top[-1,]
top<-top[!(top$Var1 %like% "fs"),]
rm(meta, counts)
allp<-allp[allp$AA.Substitutions %in% top$Var1,]

tbl0<-data.frame(table(allp$Pango.lineage,allp$AA.Substitutions))
pango<-allp[!duplicated(allp$Accession.ID),]
tbl1<-data.frame(table(pango$Pango.lineage))

merged<-merge(tbl0,tbl1,by="Var1",all.x = T)
merged$average<-(merged$Freq.x/merged$Freq.y)*100
merged<-merged[,-c(3,4)]
sprd<-spread(merged, key = Var2, value=average)
fwrite(sprd, "for_supplementary_percent_PANGO.tsv", sep="\t")
##co-occure each lineage

fwrite(allp, "allptop300.tsv",sep="\t")
###########co-occure each VOC
allp <-fread("allptop300.tsv",sep="\t")
allp <-fread("allptop.tsv",sep="\t")
VOC<-allp[allp$newpango=="Delta",]
rm(allp,VOC)
pango<-data.frame(table(VOC$AA.Substitutions, VOC$Accession.ID))
sprdpango<-spread(pango, key = Var2, value=Freq)
head(sprdpango[1:3,1:3])
row.names(sprdpango)<-sprdpango$Var1
sprdpango<-sprdpango[,-1]
mx<-data.matrix(sprdpango)
rm(sprdpango,pango)
mxt<-t(mx)
co_oc<-data.frame(crossprod(mxt), check.names = F)
m<-round(max(co_oc)/1900)
co_oc<-co_oc[rowSums(co_oc)>m,]
co_oc<-co_oc[,colSums(co_oc)>m]

library(ggplot2)
library(reshape2)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
V<-data.matrix(co_oc)
cormat <- reorder_cormat(V)
upper_tri <- get_upper_tri(cormat)
#

write.table(upper_tri, "co_occure_fortableS_Delta.tsv", sep="\t")

#Spike combinations
allptop<-fread("for_heatmap.tsv",sep="\t")
spike<-allptop[allptop$AA.Substitutions %like% "S_",]
spike<-spike[!(spike$AA.Substitutions %like% "fs"),]
spike$order<-sub("[:].*","",spike$AA.Substitutions)
spike$order<-as.numeric(sub(".*[_]","",spike$order))
head(spike$order)
spike<-spike[order(spike$order),]
aggs<-aggregate(AA.Substitutions~Accession.ID,spike,toString,collapse = ",")
merged<-merge(spike[,c(1:17,28:30)], aggs, by="Accession.ID", all.y=T)
tbl<-data.frame(table(merged$AA.Substitutions))
tbl<-tbl[tbl$Freq>2000,]
tbl1<-data.frame(table(merged$AA.Substitutions, merged$newpang))
tbl1<-tbl1[tbl1$Freq>10,]
tbl1<-aggregate(Var2~Var1,tbl1,toString,collapse = ",")
merged2<-merge(tbl, tbl1, by="Var1", all.x=T)

#genome length
#metalength<-meta[meta$`Sequence length`>29400&meta$`Sequence length` < 30100,]
meta<-fread("Out_meta_tree_IFDs.txt", sep= "\t",header = T)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.1.529"|meta$Pango.lineage %like% "BA.", "7-Omicron","0-Others")
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.1.7"|meta$Pango.lineage %like% "Q.", "1-Alpha",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.617.2"|meta$Pango.lineage %like% "AY.", "4-Delta",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.351", "2-Beta",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "B.1.6211", "6-Mu",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "P.1", "3-Gamma",meta$newpango)
meta$newpango<-ifelse(meta$Pango.lineage %like% "C.37", "5-Lambda",meta$newpango)

table(meta$newpango)

pdf('boxplot_length_lineage.pdf', width = 8, height = 6)
boxplot(meta$Sequence.length~meta$newpango,
        las=2, xlab = "",
        ylab = "",
        ylim=c(29500,30000),
        outline=F)
dev.off()
boxplot(meta$`Sequence length`~meta$Clade, 
        las=2, xlab = "",
        ylab = "",
        ylim=c(29500,30000),
        outline=F)


#anova length
res.aov2 <- aov(Sequence.length ~ newpango, data = meta)
summary(res.aov2)
#wilcox.test omicron
meta$new2 <- ifelse(meta$newpango == "Omicron",
                    meta$newpango,"Others")
data.frame(table(meta$new2))

wilcox.test(meta$Sequence.length~meta$new2)

#time split
meta$timelab<-sub("-[^-]+$","", meta$Collection.date)
meta<-meta%>%
  filter(!meta$timelab %in% c("2020", "2021","2021-1"), )
tbl<-data.frame(table(meta$timelab))

meta$recent <- ifelse(meta$timelab %like% "2019-"|
                        meta$timelab %like%"2020-",
                      "2019/2020","2021")
data.frame(table(meta$recent))
wilcox.test(meta$Sequence.length~ meta$recent)

boxplot(meta$Sequence.length~meta$timelab, 
        las=2, xlab = "",
        ylab = "",
        ylim=c(29600,30000),
        outline=F)

ggplot(meta, aes(x=newpango, y=Sequence.length, fill=newpango)) + 
  geom_violin(
    scale = "area"
  )+
  ylim(29600,30000) +
  ylab("")+
  xlab("")+
  theme_classic()+
  scale_fill_manual(name=" ",
                    values=rainbow(7),
                     labels = c('Others','Alpha', 'Beta', 'Gamma','Delta','Lambda',"Omicron"))+  
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1,
                                   size=20))+
  theme(axis.text.y =element_text(size=18))+
  scale_x_discrete( labels = c('Others','Alpha', 'Beta', 'Gamma','Delta','Lambda',"Omicron"))

#quantity indels
df<-data.frame(fread("mutations_coordinates.tsv", sep= "\t"))
df$protein<-gsub("[_].*","",df$name)
df<-df[!(df$Var1 %like% "fs\\("),]
df<-df[df$Var1 %like% "del\\(" | df$Var1 %like% "ins",]
df$length<-(df$end-df$start)+1
df$plength<-df$length/3
m<-data.frame(table(df$plength))
m<-m[m$Freq>20,]
boxplot(m$Freq~m$Var1, outline=F, las=2)
agg<-aggregate(Freq~plength,df,sum)
agg<-agg[agg$Freq>200,]
boxplot(agg$Freq~agg$plength, outline=F)






