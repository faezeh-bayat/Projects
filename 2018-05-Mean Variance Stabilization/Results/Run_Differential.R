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

#setwd("/Users/faezehbayat/Documents")
#counter=1
#path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/Model1/",counter,' .Rdata',sep="")
#path=paste("/Users/faezehbayat/Documents/",counter,' .Rdata',sep="")
#sample1_rep1 <- load(path)

differential_results <-read.table("results.bed",header=FALSE)
gene_expression=read.table("57epigenomes.RPKM.pc",header = TRUE)
gene_coordinates=read.table("gene_coordinates.txt",header = TRUE)
differential_results[,]<- NA
differential_results <- differential_results[1:7,]
rownames(differential_results)[1] = "Model"
rownames(differential_results)[2] = "GM12878-H1hESC"
rownames(differential_results)[3]="GM12878-HUVEC"
rownames(differential_results)[4]="GM12878-K562"
rownames(differential_results)[5]="H1hESC-HUVEC"
rownames(differential_results)[6]="H1hESC-K562"
rownames(differential_results)[7]="HUVEC-K562"

##################GM12878-H1--------------
for(counter in 1:24)
{
      print(counter)
      path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/Model1/",counter,' .Rdata',sep="")
      sample1_rep1 <- load(path)
      write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample1_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
      path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 H1-hESC/Model1/",counter,' .Rdata',sep="")
      sample2_rep1<- load(path)
      write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample2_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
      system("cat sample1_rep1.txt > sample1_rep1.bed")
      system("cat sample2_rep1.txt > sample2_rep1.bed")
      system("bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed")
      system("bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed")
      system("awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed")

      sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
      sample1_rep1_signals[,4] <- NULL
      sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
      sample2_rep1_signals[,4] <- NULL
      table_MA <-read.table("MAnorm.bed",header=FALSE)
      MAnorm_result <- differential_expression_analysis(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates,"GM12878","H1")
      temp <-subset(MAnorm_result,MAnorm_result[,17]!="NA")
      pearson_value=cor(temp[,14],temp[,17],method = c("pearson"))
      differential_results[1,counter] = "Model 1"
      differential_results[2,counter] = pearson_value
      write.table(differential_results,"differential_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
}
#############GM12878-HUVEC----------------
for(counter in 1:24)
{
  print(counter)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/Model1/",counter,' .Rdata',sep="")
  sample1_rep1 <- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample1_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 HUVEC/Model1/",counter,' .Rdata',sep="")
  sample2_rep1<- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample2_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  system("cat sample1_rep1.txt > sample1_rep1.bed")
  system("cat sample2_rep1.txt > sample2_rep1.bed")
  system("bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed")
  system("bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed")
  system("awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed")

  sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
  sample1_rep1_signals[,4] <- NULL
  sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
  sample2_rep1_signals[,4] <- NULL
  table_MA <-read.table("MAnorm.bed",header=FALSE)
  MAnorm_result <- differential_expression_analysis(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates,"GM12878","HUVEC")
  temp <-subset(MAnorm_result,MAnorm_result[,17]!="NA")
  pearson_value=cor(temp[,14],temp[,17],method = c("pearson"))
  differential_results[3,counter] = pearson_value
  write.table(differential_results,"differential_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
}
#############GM12878-K562
for(counter in 1:24)
{
  print(counter)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/Model1/",counter,' .Rdata',sep="")
  sample1_rep1 <- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample1_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 K562/Model1/",counter,' .Rdata',sep="")
  sample2_rep1<- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample2_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  system("cat sample1_rep1.txt > sample1_rep1.bed")
  system("cat sample2_rep1.txt > sample2_rep1.bed")
  system("bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed")
  system("bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed")
  system("awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed")

  sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
  sample1_rep1_signals[,4] <- NULL
  sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
  sample2_rep1_signals[,4] <- NULL
  table_MA <-read.table("MAnorm.bed",header=FALSE)
  MAnorm_result <- differential_expression_analysis(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates,"GM12878","K562")
  temp <-subset(MAnorm_result,MAnorm_result[,17]!="NA")
  pearson_value=cor(temp[,14],temp[,17],method = c("pearson"))
  differential_results[4,counter] = pearson_value
  write.table(differential_results,"differential_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
}
#############H1-HUVEC----------
for(counter in 1:24)
{
  print(counter)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 H1-hESC/Model1/",counter,' .Rdata',sep="")
  sample1_rep1 <- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample1_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 HUVEC/Model1/",counter,' .Rdata',sep="")
  sample2_rep1<- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample2_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  system("cat sample1_rep1.txt > sample1_rep1.bed")
  system("cat sample2_rep1.txt > sample2_rep1.bed")
  system("bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed")
  system("bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed")
  system("awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed")

  sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
  sample1_rep1_signals[,4] <- NULL
  sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
  sample2_rep1_signals[,4] <- NULL
  table_MA <-read.table("MAnorm.bed",header=FALSE)
  MAnorm_result <- differential_expression_analysis(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates,"H1","HUVEC")
  temp <-subset(MAnorm_result,MAnorm_result[,17]!="NA")
  pearson_value=cor(temp[,14],temp[,17],method = c("pearson"))
  differential_results[5,counter] = pearson_value
  write.table(differential_results,"differential_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
}
#############H1-K562-------------
for(counter in 1:24)
{
  print(counter)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 H1-hESC/Model1/",counter,' .Rdata',sep="")
  sample1_rep1 <- load()
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample1_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 K562/Model1/",counter,' .Rdata',sep="")
  sample2_rep1<- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample2_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  system("cat sample1_rep1.txt > sample1_rep1.bed")
  system("cat sample2_rep1.txt > sample2_rep1.bed")
  system("bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed")
  system("bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed")
  system("awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed")

  sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
  sample1_rep1_signals[,4] <- NULL
  sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
  sample2_rep1_signals[,4] <- NULL
  table_MA <-read.table("MAnorm.bed",header=FALSE)
  MAnorm_result <- differential_expression_analysis(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates,"H1","K562")
  temp <-subset(MAnorm_result,MAnorm_result[,17]!="NA")
  pearson_value=cor(temp[,14],temp[,17],method = c("pearson"))
  differential_results[6,counter] = pearson_value
  write.table(differential_results,"differential_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
}
#############HUVEC-K562------------
for(counter in 1:24)
{
  print(counter)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 HUVEC/Model1/",counter,' .Rdata',sep="")
  sample1_rep1 <- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample1_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  path=paste("/scratch/faezeh/Gene_expression_analysis/H3K4me3 K562/Model1/",counter,' .Rdata',sep="")
  sample2_rep1<- load(path)
  write.table(replicate1_scores,file="/scratch/faezeh/Gene_expression_analysis/H3K4me3 GM12878/sample2_rep1.txt",row.names=F, sep="\t",quote = FALSE,col.names = FALSE)
  system("cat sample1_rep1.txt > sample1_rep1.bed")
  system("cat sample2_rep1.txt > sample2_rep1.bed")
  system("bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed")
  system("bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed")
  system("awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed")

  sample1_rep1_signals <- read.table("sample1_rep1_signals.bed")
  sample1_rep1_signals[,4] <- NULL
  sample2_rep1_signals <- read.table("sample2_rep1_signals.bed")
  sample2_rep1_signals[,4] <- NULL
  table_MA <-read.table("MAnorm.bed",header=FALSE)
  MAnorm_result <- differential_expression_analysis(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates,"HUVEC","K562")
  temp <-subset(MAnorm_result,MAnorm_result[,17]!="NA")
  pearson_value=cor(temp[,14],temp[,17],method = c("pearson"))
  differential_results[7,counter] = pearson_value
  write.table(differential_results,"differential_results.xls",sep="\t",quote=FALSE,row.names=TRUE)
}


