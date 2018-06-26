setwd("/Users/faezehbayat/Documents/Experiments/2018-06-21-Gene expression analysis/fold enrichment")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(affy)
library(R.basic)
source("/Users/faezehbayat/Documents/Experiments/2018-05-22-MAnorm evaluation/Source.R")




sample1_rep1_signals <- read.table("sample1_rep1.bed")
sample2_rep1_signals <- read.table("sample2_rep1.bed")

gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_expression <- as.matrix(gene_expression)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)


Results <- gene_expression_analysis(sample1_rep1_signals,gene_expression,gene_coordinates) 
save(Results,file="Results.Rdata")
plot(Results[,11],Results[,12])
