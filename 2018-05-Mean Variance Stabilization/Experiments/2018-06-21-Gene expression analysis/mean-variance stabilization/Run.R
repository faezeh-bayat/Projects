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
start_time <- Sys.time()
results <-read.table("results.bed",header=FALSE)
results[,]<- NA
results <- results[1:8,]
rownames(results)[1] = "Dataset"
rownames(results)[2] = "Evalation type"
rownames(results)[3]="Method"
rownames(results)[4]="Alpha"
rownames(results)[5]="Width"
rownames(results)[6]="Bin size"
rownames(results)[7]="Curve"
rownames(results)[8]="Evaluation Value"



rep1<- read.table("rep1.txt")
rep2<- read.table("rep2.txt")
load("rep1_score.Rdata")
L=length(rep1_score)
load("rep2_score.Rdata")
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)


###############
scores=matrix(nrow=length(rep1_score),ncol=2)
scores=matrix(nrow=L,ncol=2)
scores[,1]  <- rep1_score[,1]
scores[,2]  <- rep2_score[,1]
x=scores[,1]
y=scores[,2]

ordered_scores <- scores[order(x,y),]
###############
curve_functions=c("linear","polynomial","step")
alpha_values=c(1,2,2^(1/300),2^(1/1000),2^(1/3000))
width_values=c(0,1,2,3)
bin_sizes=c(1000,10000,500)
counter=1
pdf("PLOTS.pdf")
for(m in 2:(length(bin_sizes)))
  for(i in 1:length(curve_functions))
    for(j in 1:length(alpha_values))
      for(k in 1:length(width_values))
        
          {
          print(counter)
          mean_variance <- Weighted_Mean_Variance(rep1,rep2,rep1_score,rep2_score,width_values[k],bin_sizes[m],alpha_values[j],L,ordered_scores)
          replicate1_scores <- score_calculation(mean_variance,rep1,curve_functions[i])
          rep1_scores <- replicate1_scores
          rep1_scores[,4] <- NULL 
          Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates)
          temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
          ordered_x<- temp[order(temp[,11]),]
          zoom_in_x <- ordered_x[1:(length(temp[,1])-0),]
          ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
          total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-0),]
          pearson_value=cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))
          q=plot(total_zoom_in[,11],total_zoom_in[,12],main=paste("Curve function",curve_functions[i],"alpha",alpha_values[j],"width",width_values[k],"bin size",bin_sizes[m], sep=" "),xlab = "Signal at TSS",ylab = "Gene expression")
          print(q)
          results[1,counter] = "HUVEC"
          results[2,counter] = "Gene expression analysis"
          results[3,counter]="Variance stabilization"
          results[4,counter]=alpha_values[j]
          results[5,counter]=width_values[k]
          results[6,counter]=bin_sizes[m]
          results[7,counter]=curve_functions[i]
          results[8,counter]=pearson_value
          counter=counter+1
          #abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")
          write.table(results,"results.xls",sep="\t",quote=FALSE,row.names=TRUE)
          }

dev.off()
end_time <- Sys.time()
taken=end_time-(start_time)
taken
