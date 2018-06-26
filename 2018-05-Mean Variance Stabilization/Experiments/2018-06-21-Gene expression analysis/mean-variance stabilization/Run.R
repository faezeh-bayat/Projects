setwd("Projects/2018-05-Mean Variance Stabilization/Experiments/2018-06-21-Gene expression analysis/mean-variance stabilization/Run.R")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(affy)
library(R.basic)
source("Projects/2018-05-Mean Variance Stabilization/Source/Source.R")


load("sample1_rep1.Rdata")
load("sample1_rep2.Rdata")
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_expression <- as.matrix(gene_expression)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)


Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
Sample1_rep2_results <- gene_expression_analysis(rep2_scores,gene_expression,gene_coordinates)  #sample 1 replicate 2
save(Sample1_rep1_results,file="Sample1_rep1_results.Rdata")
save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
plot(Sample1_rep1_results[,11],Sample1_rep1_results[,12])
plot(Sample1_rep2_results[,11],Sample1_rep2_results[,12])
