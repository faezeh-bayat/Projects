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

#==============mean-variance replicate1============

#sample1_rep1_signals <- read.table("sample1_rep1.bed")
#sample2_rep1_signals <- read.table("sample2_rep1.bed")

load("replicate1_scores.Rdata")
rep1_scores <- replicate1_scores
rep1_scores[,4] <- NULL
#load("sample1_rep2.Rdata")
#rep2_scores[,4] <- NULL
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)

Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
#Sample1_rep2_results <- gene_expression_analysis(rep2_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
save(Sample1_rep1_results,file="Sample1_rep1_results.Rdata")
#save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
#save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
#plot(Sample1_rep2_results[,11],Sample1_rep2_results[,12])

load("Sample1_rep1_results.Rdata")
temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
ordered_x<- temp[order(temp[,11]),]
zoom_in_x <- ordered_x[1:(length(temp[,1])-10),]
ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-10),]
plot(total_zoom_in[,11],total_zoom_in[,12],main = "H3K4me3 GM12878 gene expression analysis(mean-variance)",xlab = "Signal at TSS",ylab = "Gene expression")
abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")

cor(total_zoom_in[,11],total_zoom_in[,12],method = c("spearman"))
cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))



#===========mean-variance replicate 2==============
load("replicate2_scores.Rdata")
rep2_scores <- replicate2_scores
rep2_scores[,4] <- NULL
#load("sample1_rep2.Rdata")
#rep2_scores[,4] <- NULL
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)

Sample1_rep2_results <- gene_expression_analysis(rep2_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
#Sample1_rep2_results <- gene_expression_analysis(rep2_scores,gene_expression,gene_coordinates)  #sample 1 replicate 1
save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
#save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
#save(Sample1_rep2_results,file="Sample1_rep2_results.Rdata")
#plot(Sample1_rep2_results[,11],Sample1_rep2_results[,12])

load("Sample1_rep2_results.Rdata")
temp <-subset(Sample1_rep2_results,Sample1_rep2_results[,12]!="NA")
ordered_x<- temp[order(temp[,11]),]
zoom_in_x <- ordered_x[1:(length(temp[,1])-10),]
ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-10),]
plot(total_zoom_in[,11],total_zoom_in[,12],main = "H3K4me3 GM12878 gene expression analysis(mean-variance)",xlab = "Signal at TSS",ylab = "Gene expression")
abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")

cor(total_zoom_in[,11],total_zoom_in[,12],method = c("spearman"))
cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))

