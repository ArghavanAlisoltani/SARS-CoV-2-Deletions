# neural network
install.packages("neuralnet")
library(neuralnet)
met<-read.delim(file="2.A.proteinFamily.txt", sep="\t")
met2<-read.delim(file="2.B.proteinFamily.txt", sep="\t")
mytable<-data.frame(table(met$X..PRTases_typeI))

mytable2<-data.frame(table(met2$))
write.csv(mytable2, "virgo_anno_Interpro_table.csv")
