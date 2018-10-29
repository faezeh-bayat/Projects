
setwd("/Users/faezehbayat/Documents/new experiments/gene expression analysis")
rm(list=ls())
library(ggplot2)
data <- read.csv("total results.csv")
ggplot(data,aes(x=Dataset,y=Pearson.Correlation,fill=Method))+geom_bar(position = "dodge2")
x=data$Dataset
y=data$Pearson.Correlation
grp=data$Method


ggplot(data, aes(x, y, group = grp)) + geom_col(aes(fill = grp), position = "dodge2") +theme_bw()
ggplot(data, aes(x, y, group = grp)) + geom_col(aes(fill = grp), position = "dodge2") +theme_bw()+
scale_fill_manual("legend", values = c("Variance stabilization" = "red", "Fold enrichment" = "orange", "log(fold enrichment)" = "blue", "log(p-value)"="pink","P-value"="green"))
ggplot(data, aes(x, y, group = grp)) + geom_col(aes(fill = grp), position = "dodge2") +theme_bw()+
scale_fill_brewer(palette="Dark2")+labs(x = "\n\nChIP-seq datasets",y="\nPearson correlation")

