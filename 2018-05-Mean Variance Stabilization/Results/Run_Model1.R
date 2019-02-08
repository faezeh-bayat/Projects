setwd("/scratch/faezeh/Gene_expression_analysis/H3K4me3 K562")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(Weighted.Desc.Stat)
library(stats)
source("/scratch/faezeh/Gene_expression_analysis/Source.R")

start_time <- Sys.time()
results <-read.table("results.bed",header=FALSE)
results[,]<- NA
results <- results[1:9,]
rownames(results)[1] = "Dataset"
rownames(results)[2] = "Evalation type"
rownames(results)[3]="Method"
rownames(results)[4]="Alpha"
rownames(results)[5]="Beta"
rownames(results)[6]="Width"
rownames(results)[7]="Bin size"
rownames(results)[8]="Curve"
rownames(results)[9]="Evaluation Value"

rep1<- read.table("rep1.txt")
rep2<- read.table("rep2.txt")
load("rep1_score.Rdata")
load("rep2_score.Rdata")
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)
a=c(rep1_score,rep2_score)
b=c(rep2_score,rep1_score)
scores=matrix(nrow=length(a),ncol=2)
scores[,1]  <- a
scores[,2]  <- b
x=scores[,1]
y=scores[,2]
ordered_scores <- scores[order(x,y),]
L=length(a)
#curve_functions=c("linear","polynomial","step")
curve_functions=c("step")
#alpha_values=c(1,2,2^(1/300),2^(1/1000),2^(1/3000))
beta_values=c(Inf,1,300,1000,3000,1e4,1e5,1e6)
width_values=c()
alpha_values=c()
bin_sizes=c(10000,100000,1000)
for(ll in 1:length(beta_values))
{
  alpha_values[ll]=2^(1/beta_values[ll])
  width_values[ll]=ceiling(log(1/0.01)/log(alpha_values[ll]))
}
width_values[!is.finite(width_values)] <- 0
#a=length(width_values)
#width_values[a+1]=1
#width_values[a+2]=3
counter=1
#pdf("PLOTS.pdf")
for(m in 1:(length(bin_sizes)))
  for(i in 1:length(curve_functions))
    for(j in 1:length(alpha_values))
      #for(k in 1:length(width_values))
      {
        print(counter)
        #mean_variance <- Model3(rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
        mean_variance<-Weighted_Mean_Variance(rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
        without_INF<- subset(mean_variance,mean_variance[,2]!=Inf & mean_variance[,2]!=-Inf)
        a=max(without_INF[,2])
        mean_variance[!is.finite(mean_variance[,2])] <- a
        replicate1_scores <- score_calculation(mean_variance,rep1,curve_functions[i])
        save(replicate1_scores, file=paste(counter, '.Rdata') )
        rep1_scores <- replicate1_scores
        rep1_scores[,4] <- NULL 
        Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates,"K562")
        temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
        ordered_x<- temp[order(temp[,11]),]
        zoom_in_x <- ordered_x[1:(length(temp[,1])-0),]
        ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
        total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-0),]
        pearson_value=cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))
        #q=plot(total_zoom_in[,11],total_zoom_in[,12],main=paste("Curve function",curve_functions[i],"alpha",alpha_values[j],"width",width_values[k],"bin size",bin_sizes[m], sep=" "),xlab = "Signal at TSS",ylab = "Gene expression")
        #print(q)
        results[1,counter] = "K562"
        results[2,counter] = "Gene expression analysis"
        results[3,counter]="Variance stabilization"
        results[4,counter]=alpha_values[j]
        results[5,counter]=beta_values[j]
        results[6,counter]=width_values[j]
        results[7,counter]=bin_sizes[m]
        results[8,counter]=curve_functions[i]
        results[9,counter]=pearson_value
        counter=counter+1
        #abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")
        write.table(results,"results.xls",sep="\t",quote=FALSE,row.names=TRUE)
      }

#dev.off()
end_time <- Sys.time()
taken=end_time-(start_time)
taken



########################################################################

setwd("/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(Weighted.Desc.Stat)
library(stats)
source("/scratch/faezeh/Gene_expression_analysis/Source.R")

start_time <- Sys.time()
results <-read.table("results.bed",header=FALSE)
results[,]<- NA
results <- results[1:9,]
rownames(results)[1] = "Dataset"
rownames(results)[2] = "Evalation type"
rownames(results)[3]="Method"
rownames(results)[4]="Alpha"
rownames(results)[5]="Beta"
rownames(results)[6]="Width"
rownames(results)[7]="Bin size"
rownames(results)[8]="Curve"
rownames(results)[9]="Evaluation Value"

rep1<- read.table("rep1.txt")
rep2<- read.table("rep2.txt")
load("rep1_score.Rdata")
load("rep2_score.Rdata")
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)
a=c(rep1_score,rep2_score)
b=c(rep2_score,rep1_score)
scores=matrix(nrow=length(a),ncol=2)
scores[,1]  <- a
scores[,2]  <- b
x=scores[,1]
y=scores[,2]
ordered_scores <- scores[order(x,y),]
L=length(a)
#curve_functions=c("linear","polynomial","step")
curve_functions=c("step")
#alpha_values=c(1,2,2^(1/300),2^(1/1000),2^(1/3000))
beta_values=c(Inf,1,300,1000,3000,1e4,1e5,1e6)
width_values=c()
alpha_values=c()
bin_sizes=c(10000,100000,1000)
for(ll in 1:length(beta_values))
{
  alpha_values[ll]=2^(1/beta_values[ll])
  width_values[ll]=ceiling(log(1/0.01)/log(alpha_values[ll]))
}
width_values[!is.finite(width_values)] <- 0
#a=length(width_values)
#width_values[a+1]=1
#width_values[a+2]=3
counter=1
#pdf("PLOTS.pdf")
for(m in 1:(length(bin_sizes)))
  for(i in 1:length(curve_functions))
    for(j in 1:length(alpha_values))
      #for(k in 1:length(width_values))
    {
      print(counter)
      #mean_variance <- Model3(rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      mean_variance <-Weighted_Mean_Variance (rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      without_INF<- subset(mean_variance,mean_variance[,2]!=Inf & mean_variance[,2]!=-Inf)
      a=max(without_INF[,2])
      mean_variance[!is.finite(mean_variance[,2])] <- a
      replicate1_scores <- score_calculation(mean_variance,rep1,curve_functions[i])
      save(replicate1_scores, file=paste(counter, '.Rdata') )
      rep1_scores <- replicate1_scores
      rep1_scores[,4] <- NULL 
      Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates,"GM12878")
      temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
      ordered_x<- temp[order(temp[,11]),]
      zoom_in_x <- ordered_x[1:(length(temp[,1])-0),]
      ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
      total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-0),]
      pearson_value=cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))
      #q=plot(total_zoom_in[,11],total_zoom_in[,12],main=paste("Curve function",curve_functions[i],"alpha",alpha_values[j],"width",width_values[k],"bin size",bin_sizes[m], sep=" "),xlab = "Signal at TSS",ylab = "Gene expression")
      #print(q)
      results[1,counter] = "GM12878"
      results[2,counter] = "Gene expression analysis"
      results[3,counter]="Variance stabilization"
      results[4,counter]=alpha_values[j]
      results[5,counter]=beta_values[j]
      results[6,counter]=width_values[j]
      results[7,counter]=bin_sizes[m]
      results[8,counter]=curve_functions[i]
      results[9,counter]=pearson_value
      counter=counter+1
      #abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")
      write.table(results,"results.xls",sep="\t",quote=FALSE,row.names=TRUE)
    }

#dev.off()
end_time <- Sys.time()
taken=end_time-(start_time)
taken



##############################################

setwd("/scratch/faezeh/Gene_expression_analysis/H3K4me3 H1-hESC")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(Weighted.Desc.Stat)
library(stats)
source("/scratch/faezeh/Gene_expression_analysis/Source.R")

start_time <- Sys.time()
results <-read.table("results.bed",header=FALSE)
results[,]<- NA
results <- results[1:9,]
rownames(results)[1] = "Dataset"
rownames(results)[2] = "Evalation type"
rownames(results)[3]="Method"
rownames(results)[4]="Alpha"
rownames(results)[5]="Beta"
rownames(results)[6]="Width"
rownames(results)[7]="Bin size"
rownames(results)[8]="Curve"
rownames(results)[9]="Evaluation Value"

rep1<- read.table("rep1.txt")
rep2<- read.table("rep2.txt")
load("rep1_score.Rdata")
load("rep2_score.Rdata")
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)
a=c(rep1_score,rep2_score)
b=c(rep2_score,rep1_score)
scores=matrix(nrow=length(a),ncol=2)
scores[,1]  <- a
scores[,2]  <- b
x=scores[,1]
y=scores[,2]
ordered_scores <- scores[order(x,y),]
L=length(a)
#curve_functions=c("linear","polynomial","step")
curve_functions=c("step")
#alpha_values=c(1,2,2^(1/300),2^(1/1000),2^(1/3000))
beta_values=c(Inf,1,300,1000,3000,1e4,1e5,1e6)
width_values=c()
alpha_values=c()
bin_sizes=c(10000,100000,1000)
for(ll in 1:length(beta_values))
{
  alpha_values[ll]=2^(1/beta_values[ll])
  width_values[ll]=ceiling(log(1/0.01)/log(alpha_values[ll]))
}
width_values[!is.finite(width_values)] <- 0
#a=length(width_values)
#width_values[a+1]=1
#width_values[a+2]=3
counter=1
#pdf("PLOTS.pdf")
for(m in 1:(length(bin_sizes)))
  for(i in 1:length(curve_functions))
    for(j in 1:length(alpha_values))
      #for(k in 1:length(width_values))
    {
      print(counter)
      #mean_variance <- Model3(rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      mean_variance <-Weighted_Mean_Variance (rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      without_INF<- subset(mean_variance,mean_variance[,2]!=Inf & mean_variance[,2]!=-Inf)
      a=max(without_INF[,2])
      mean_variance[!is.finite(mean_variance[,2])] <- a
      replicate1_scores <- score_calculation(mean_variance,rep1,curve_functions[i])
      save(replicate1_scores, file=paste(counter, '.Rdata') )
      rep1_scores <- replicate1_scores
      rep1_scores[,4] <- NULL 
      Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates,"H1")
      temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
      ordered_x<- temp[order(temp[,11]),]
      zoom_in_x <- ordered_x[1:(length(temp[,1])-0),]
      ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
      total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-0),]
      pearson_value=cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))
      #q=plot(total_zoom_in[,11],total_zoom_in[,12],main=paste("Curve function",curve_functions[i],"alpha",alpha_values[j],"width",width_values[k],"bin size",bin_sizes[m], sep=" "),xlab = "Signal at TSS",ylab = "Gene expression")
      #print(q)
      results[1,counter] = "H1-hESC"
      results[2,counter] = "Gene expression analysis"
      results[3,counter]="Variance stabilization"
      results[4,counter]=alpha_values[j]
      results[5,counter]=beta_values[j]
      results[6,counter]=width_values[j]
      results[7,counter]=bin_sizes[m]
      results[8,counter]=curve_functions[i]
      results[9,counter]=pearson_value
      counter=counter+1
      #abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")
      write.table(results,"results.xls",sep="\t",quote=FALSE,row.names=TRUE)
    }

#dev.off()
end_time <- Sys.time()
taken=end_time-(start_time)
taken



################################

setwd("/scratch/faezeh/Gene_expression_analysis/H3K4me3 HUVEC")
rm(list=ls())
library(pracma)
library(readxl)
library(MASS)
library(data.table)
library(gdata)
library(Weighted.Desc.Stat)
library(stats)
source("/scratch/faezeh/Gene_expression_analysis/Source.R")

start_time <- Sys.time()
results <-read.table("results.bed",header=FALSE)
results[,]<- NA
results <- results[1:9,]
rownames(results)[1] = "Dataset"
rownames(results)[2] = "Evalation type"
rownames(results)[3]="Method"
rownames(results)[4]="Alpha"
rownames(results)[5]="Beta"
rownames(results)[6]="Width"
rownames(results)[7]="Bin size"
rownames(results)[8]="Curve"
rownames(results)[9]="Evaluation Value"

rep1<- read.table("rep1.txt")
rep2<- read.table("rep2.txt")
load("rep1_score.Rdata")
load("rep2_score.Rdata")
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)
a=c(rep1_score,rep2_score)
b=c(rep2_score,rep1_score)
scores=matrix(nrow=length(a),ncol=2)
scores[,1]  <- a
scores[,2]  <- b
x=scores[,1]
y=scores[,2]
ordered_scores <- scores[order(x,y),]
L=length(a)
#curve_functions=c("linear","polynomial","step")
curve_functions=c("step")
#alpha_values=c(1,2,2^(1/300),2^(1/1000),2^(1/3000))
beta_values=c(Inf,1,300,1000,3000,1e4,1e5,1e6)
width_values=c()
alpha_values=c()
bin_sizes=c(10000,1000,100000)
for(ll in 1:length(beta_values))
{
  alpha_values[ll]=2^(1/beta_values[ll])
  width_values[ll]=ceiling(log(1/0.01)/log(alpha_values[ll]))
}
width_values[!is.finite(width_values)] <- 0
#a=length(width_values)
#width_values[a+1]=1
#width_values[a+2]=3
counter=1
#pdf("PLOTS.pdf")
for(m in 1:(length(bin_sizes)))
  for(i in 1:length(curve_functions))
    for(j in 1:length(alpha_values))
      #for(k in 1:length(width_values))
    {
      print(counter)
      #mean_variance <- Model3(rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      mean_variance <-Weighted_Mean_Variance (rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      without_INF<- subset(mean_variance,mean_variance[,2]!=Inf & mean_variance[,2]!=-Inf)
      a=max(without_INF[,2])
      mean_variance[!is.finite(mean_variance[,2])] <- a
      replicate1_scores <- score_calculation(mean_variance,rep1,curve_functions[i])
      save(replicate1_scores, file=paste(counter, '.Rdata') )
      rep1_scores <- replicate1_scores
      rep1_scores[,4] <- NULL 
      Sample1_rep1_results <- gene_expression_analysis(rep1_scores,gene_expression,gene_coordinates,"HUVEC")
      temp <-subset(Sample1_rep1_results,Sample1_rep1_results[,12]!="NA")
      ordered_x<- temp[order(temp[,11]),]
      zoom_in_x <- ordered_x[1:(length(temp[,1])-0),]
      ordered_y <- zoom_in_x[order(zoom_in_x[,12]),]
      total_zoom_in <- ordered_y[1:(length(ordered_y[,1])-0),]
      pearson_value=cor(total_zoom_in[,11],total_zoom_in[,12],method = c("pearson"))
      #q=plot(total_zoom_in[,11],total_zoom_in[,12],main=paste("Curve function",curve_functions[i],"alpha",alpha_values[j],"width",width_values[k],"bin size",bin_sizes[m], sep=" "),xlab = "Signal at TSS",ylab = "Gene expression")
      #print(q)
      results[1,counter] = "HUVEC"
      results[2,counter] = "Gene expression analysis"
      results[3,counter]="Variance stabilization"
      results[4,counter]=alpha_values[j]
      results[5,counter]=beta_values[j]
      results[6,counter]=width_values[j]
      results[7,counter]=bin_sizes[m]
      results[8,counter]=curve_functions[i]
      results[9,counter]=pearson_value
      counter=counter+1
      #abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")
      write.table(results,"results.xls",sep="\t",quote=FALSE,row.names=TRUE)
    }

#dev.off()
end_time <- Sys.time()
taken=end_time-(start_time)
taken





