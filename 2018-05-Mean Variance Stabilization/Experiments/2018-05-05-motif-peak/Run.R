setwd("H:/Experiments/2018-05-05-motif-peak/CTCF")
rm(list=ls())
library(pracma)
source("/Users/faezehbayat/Documents/Experiments/2018-05-22-MAnorm evaluation/MAnorm/Source.R")

rep1<- read.table("rep1.txt") #replicate1 fold enrichment data
rep2<- read.table("rep2.txt")
motifs <- read.table("motifs.txt")


meanvar_scores <- Mean_Variance(rep1,rep2)
rep1_scores <- rep1_meanvar_score(meanvar_scores,rep1)
rep2_scores <- rep2_meanvar_score(meanvar_scores,rep2)
rep1_occurance <- rep1_evaluation(motifs,rep1_scores)
rep2_occurance <- rep2_evaluation(motifs,rep2_scores)





