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
log_results <-read.table("results.bed",header=FALSE)
log_results[,]<- NA
log_results <- log_results[1:9,]
rownames(log_results)[1] = "Dataset"
rownames(log_results)[2] = "Evalation type"
rownames(log_results)[3]="Method"
rownames(log_results)[4]="Alpha"
rownames(log_results)[5]="Beta"
rownames(log_results)[6]="Width"
rownames(log_results)[7]="Bin size"
rownames(log_results)[8]="Curve"
rownames(log_results)[9]="Evaluation Value"

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
beta_values=c(Inf,1,300,1000,3000,1e4,1e5,1e6,1e7)
width_values=c()
alpha_values=c()
bin_sizes=c(100000,10000,1000)
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
diff_counter=1
#pdf("PLOTS.pdf")
ll=length(bin_sizes)*length(curve_functions)*length(alpha_values)
for(m in 1:(length(bin_sizes)))
  for(i in 1:length(curve_functions))
    for(j in 1:length(alpha_values))
      #for(k in 1:length(width_values))
    {
      print(counter)
      mean_variance<-Weighted_Mean_Variance(rep1,rep2,width_values[j],bin_sizes[m],alpha_values[j],ordered_scores)
      without_INF<- subset(mean_variance,mean_variance[,2]!=Inf & mean_variance[,2]!=-Inf)
      a=max(without_INF[,2])
      mean_variance[!is.finite(mean_variance[,2])] <- a
      log_likelihood=0
      log_log_likelihood=0
      asinh_log_likelihood=0
      fold_log_likelihood=0
      for (k in 1:L)
      {
        nearest=which.min(abs(mean_variance[,1] - scores[k,1]))
        
        sigma=1/mean_variance[nearest,2]
        log_sigma=scores[k,2]+1
        asinh_sigma=sqrt((scores[k,2]+1)^2+1)
        fold_sigma=1
        
        log_like=dnorm(scores[k,2],scores[k,1],sigma,log=TRUE)
        log_likelihood=log_likelihood+log_like
        
        log_log_like=dnorm(scores[k,2],scores[k,1],log_sigma,log=TRUE)
        log_log_likelihood=log_log_likelihood+log_log_like
        
        asinh_log_like=dnorm(scores[k,2],scores[k,1],asinh_sigma,log=TRUE)
        asinh_log_likelihood=asinh_log_likelihood+asinh_log_like
        
        fold_log_like=dnorm(scores[k,2],scores[k,1],fold_sigma,log=TRUE)
        fold_log_likelihood=fold_log_likelihood+fold_log_like
      }
      log_likelihood=log_likelihood/L
      log_log_likelihood=log_log_likelihood/L
      asinh_log_likelihood=asinh_log_likelihood/L
      fold_log_likelihood=fold_log_likelihood/L
      
      
      log_results[1,counter] = "GM12878"
      log_results[2,counter] = "Variance Stabilization Log Likelihood"
      log_results[3,counter]="Variance stabilization"
      log_results[4,counter]=alpha_values[j]
      log_results[5,counter]=beta_values[j]
      log_results[6,counter]=width_values[j]
      log_results[7,counter]=bin_sizes[m]
      log_results[8,counter]=curve_functions[i]
      log_results[9,counter]=log_likelihood
      
      ########
      log_results[1,counter+(1*ll)] = "GM12878"
      log_results[2,counter+(1*ll)] = "log(x+1) Log Likelihood"
      log_results[3,counter+(1*ll)]="Variance stabilization"
      log_results[4,counter+(1*ll)]=alpha_values[j]
      log_results[5,counter+(1*ll)]=beta_values[j]
      log_results[6,counter+(1*ll)]=width_values[j]
      log_results[7,counter+(1*ll)]=bin_sizes[m]
      log_results[8,counter+(1*ll)]=curve_functions[i]
      log_results[9,counter+(1*ll)]=log_log_likelihood
      ##########
      log_results[1,counter+(2*ll)] = "GM12878"
      log_results[2,counter+(2*ll)] = "asinh Log Likelihood"
      log_results[3,counter+(2*ll)]="Variance stabilization"
      log_results[4,counter+(2*ll)]=alpha_values[j]
      log_results[5,counter+(2*ll)]=beta_values[j]
      log_results[6,counter+(2*ll)]=width_values[j]
      log_results[7,counter+(2*ll)]=bin_sizes[m]
      log_results[8,counter+(2*ll)]=curve_functions[i]
      log_results[9,counter+(2*ll)]=asinh_log_likelihood
      ###########
      log_results[1,counter+(3*ll)] = "GM12878"
      log_results[2,counter+(3*ll)] = "Fold enrichment Log Likelihood"
      log_results[3,counter+(3*ll)]="Variance stabilization"
      log_results[4,counter+(3*ll)]=alpha_values[j]
      log_results[5,counter+(3*ll)]=beta_values[j]
      log_results[6,counter+(3*ll)]=width_values[j]
      log_results[7,counter+(3*ll)]=bin_sizes[m]
      log_results[8,counter+(3*ll)]=curve_functions[i]
      log_results[9,counter+(3*ll)]=fold_log_likelihood
      counter=counter+1
      #abline(rlm(total_zoom_in[,12]~total_zoom_in[,11]),col="red")
      write.table(log_results,"log_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
    }

#dev.off()
end_time <- Sys.time()
taken=end_time-(start_time)
taken



