setwd("/Users/faezehbayat/Documents/Experiments/2018-05-22-MAnorm evaluation/fold enrichment")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(affy)
library(R.basic)
source("Projects/2018-05-Mean Variance Stabilization/Source/Source.R")


#call Run.sh for getting the genomic position intersections
sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
table_MA <-read.table("MAnorm.bed",header=FALSE)
MAnorm_result <- MAnorm(sample1_rep1_signals,sample2_rep1_signals,table_MA)#Calling MAnorm function for calculating M and A values for each peak

MAnorm_result <- as.data.frame(fread("MAnorm_result.xls"))
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_expression <- as.matrix(gene_expression)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)

Results <- classify(MAnorm_result,gene_expression,gene_coordinates) #Calling classify function to classify peaks based on their M-value and finds their nearby peaks



