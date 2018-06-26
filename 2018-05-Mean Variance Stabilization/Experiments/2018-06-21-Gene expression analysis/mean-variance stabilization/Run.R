setwd("/Users/faezehbayat/Documents/Experiments/2018-06-21-Gene expression analysis/mean-variance stabilization")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(affy)
library(R.basic)
source("/Users/faezehbayat/Documents/Experiments/2018-05-22-MAnorm evaluation/Source.R")




#sample1_rep1_signals <- read.table("sample1_rep1.bed")
#sample2_rep1_signals <- read.table("sample2_rep1.bed")

load("sample1_rep1.Rdata")
rep1_scores[,4] <- NULL
load("sample1_rep2.Rdata")
rep2_scores[,4] <- NULL
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_expression <- as.matrix(gene_expression)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)

#Results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
Sample1_rep2_results <- gene_expression_analysis(rep2_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
plot(Sample1_rep2_results[,11],Sample1_rep2_results[,12])
