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

rep1<- read.table("rep1.txt")
rep2<- read.table("rep2.txt")
load("rep1_score.Rdata")
load("rep2_score.Rdata")
mean_variance <- Weighted_Mean_Variance(rep1,rep2,rep1_score,rep2_score,10,1000,1.000069)
replicate1_scores <- linear_function_score(mean_variance,rep1)
rep1_scores <- replicate1_scores
rep1_scores[,4] <- NULL
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)
Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates)
temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
ordered_x<- temp[order(temp[,11]),]
zoom_in_x <- ordered_x[1:(length(temp[,1])),]
ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
total_zoom_in <- ordered_y[1:(length(ordered_y[,1])),]
cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))



############testing more combinations of parameters##########
