Mean_Variance<-function(rep1,rep2)
{
  rep1_score=matrix(nrow=rep1[length(rep1[,1]),3],ncol=1) #Scores for all genomic positions
  rep2_score=matrix(nrow=rep1[length(rep1[,1]),3],ncol=1)
  
  for(i in 1:length(rep1[,4]))
  {
    print(i)
    rep1_score[rep1[i,2]:rep1[i,3],1]<-rep1[i,4]
    
  }
  
  
  for(i in 1:length(rep2[,4]))
  {
    print(i)
    rep2_score[rep2[i,2]:rep2[i,3],1]<-rep2[i,4]
    
  }
  
  
  scores=matrix(nrow=length(rep1_score),ncol=2)
  scores[,1]  <- rep1_score[,1]
  scores[,2]  <- rep2_score[,1]
  
  x=scores[,1]
  y=scores[,2]
  
  ordered_scores <- scores[order(x,y),]
  
  segment_points=matrix(nrow=1000,ncol=1)
  l=floor(length(ordered_scores[,1])/1000)
  mean_var=matrix(nrow=l,ncol=2)
  
  
  start=1
  end=1000
  
  
  for (i in 1:(l-1))
  {
    print(i)
    segment_points[1:1000,1]<-ordered_scores[start:end,2]
    #segment_points[1001:2000,1]<-ordered_scores[start:end,2]
    mean_var[i,1]<-mean(segment_points)
    mean_var[i,2]<-var(segment_points)
    start<- end+1
    end <- end+1000
  }
  
  endpoint=length(ordered_scores[,1])
  
  mean <- mean(ordered_scores[start:endpoint,1])
  var <- var(ordered_scores[start:endpoint,1])
  mean_var[l,1]<-mean
  mean_var[l,2]<-var
  
  
  l=length(mean_var[,1])
  #plot(mean_var[,1],mean_var[,2])
  
  mean_variance=matrix(nrow=(l-1),ncol=2)
  
  mean_variance[,1] <- mean_var[1:(l-1),1]
  mean_variance[,2] <- mean_var[1:(l-1),2]
  
  save(mean_var,file="mean_var.Rdata")
  save(mean_variance,file="mean_variance.Rdata")
  #return(mean_variance)
  meanvar_scores <- mean_variance
  rep1_scores <- rep1_meanvar_score(meanvar_scores,rep1)
  rep2_scores <- rep2_meanvar_score(meanvar_scores,rep2)
  rep1_occurance <- rep1_evaluation(motifs,rep1_scores)
  rep1_occurance
}

#########################################################
rep1_meanvar_score<-function(mean_variance,rep1)
{
  data<- mean_variance
  
  x=data[,1]
  y=data[,2]
  
  ordered_data <- data[order(x,y),]
  
  base_integeal<-trapz(ordered_data[,1],ordered_data[,2])
  
  replicate1_scores <- rep1
  min_mean=min(x)
  max_mean=max(x)
  min_var=min(y)
  max_var=max(y)
  
  l=length(rep1[,4])
  for (i in 1:l)
  {
    print(i)
    
    if(rep1[i,4]==0)
      replicate1_scores[i,5]<- 0
    
    else if(rep1[i,4]<=max_mean)
    {
      sub <- subset(ordered_data, ordered_data[,1]<=rep1[i,4])
      replicate1_scores[i,5] <- trapz(sub[,1],sub[,2])
    }
    else if(rep1[i,4]>max_mean)
    {
      replicate1_scores[i,5] <- base_integeal+(max_var*(rep1[i,4]-max_mean))
    }
    
  }
  replicate1_scores
}
#############################################################
rep2_meanvar_score<-function(mean_variance,rep2)
{
  data<- mean_variance
  
  x=data[,1]
  y=data[,2]
  
  ordered_data <- data[order(x,y),]
  
  base_integeal<-trapz(ordered_data[,1],ordered_data[,2])
  
  replicate1_scores <- rep1
  replicate2_scores <- rep2
  
  min_mean=min(x)
  max_mean=max(x)
  min_var=min(y)
  max_var=max(y)
  
  ll=length(rep2[,4])
  for (i in 1:ll)
  {
    print(i)
    
    if(rep2[i,4]==0)
      replicate2_scores[i,5]<- 0
    
    else if(rep2[i,4]<=max_mean)
    {
      sub <- subset(ordered_data, ordered_data[,1]<=rep2[i,4])
      replicate2_scores[i,5] <- trapz(sub[,1],sub[,2])
    }
    else if(rep2[i,4]>max_mean)
    {
      replicate2_scores[i,5] <- base_integeal+(max_var*(rep2[i,4]-max_mean))
    }
  }
  
  #save(replicate1_scores,file="replicate1_scores.Rdata")
  #trapz(ordered_data[,1],ordered_data[,2])
  replicate2_scores
}
############################################################################
rep1_evaluation <- function(motifs,replicate1_scores)
{
  
  whole_motifs <- motifs
  chr21_motifs <- subset(whole_motifs,whole_motifs[,1]=="chr21")# && whole_motifs[,4]=="NRSF"))
  chr21_NRSF_motifs <- subset(chr21_motifs,chr21_motifs[,4]=="REST")# && whole_motifs[,4]=="NRSF"))
  
  l=length(chr21_NRSF_motifs[,1])
  
  for (i in 1:l)
  {
    print(i)
    chr21_NRSF_motifs[i,7] <-chr21_NRSF_motifs[i,2]+(ceiling((chr21_NRSF_motifs[i,3]-chr21_NRSF_motifs[i,2])/2))
  }
  
  ordered_replicate1_scores <- replicate1_scores[order(replicate1_scores[,5],decreasing = TRUE),]
  
  
  l=length(chr21_NRSF_motifs[,1])
  replicate1_motifs<-data.frame(NULL,NULL,NULL)
  
  counter=0
  for (i in 1:3000)
  {
    print(i)
    rep1_seq=seq(from=ordered_replicate1_scores[i,2],to=ordered_replicate1_scores[i,3])
    for(j in 1:l)
    {
      motif_seq=seq(from=(chr21_NRSF_motifs[j,7]-50),to=(chr21_NRSF_motifs[j,7]+50))
      intersection=intersect(rep1_seq,motif_seq)
      if (length(intersection)!=0)
      {
        counter<-counter+1
        replicate1_motifs[counter,1]<-ordered_replicate1_scores[i,2]
        replicate1_motifs[counter,2]<-ordered_replicate1_scores[i,3]
        replicate1_motifs[counter,3]<-1
        replicate1_motifs[counter,4]<-ordered_replicate1_scores[i,5]
        break
      }
      
    }
    
  }
  counter
}

#===============functions for MAnorm evaluation==========
nearset_peak<- function(tss,gene_coordinates,gene_expression)
{
  upstream=tss-8000
  downstream=tss+2000
  nearest_genes<-subset(gene_coordinates,(gene_coordinates[,3]>=upstream && gene_coordinates[,3]<=downstream))
  if (length(nearest_genes[,1])!=0)
  {
    nearest_genes[,9]<- abs(tss-(nearest_genes[,3]))
    min_index <- which.min(nearest_genes[,9])
    nearest_gene <- nearest_genes[min_index,]
    nearest_gene<- as.matrix(nearest_gene)
    gene_expression_index <- which(gene_expression[,1]==nearest_gene[1,1])
    return (gene_expression_index)
  }
  else
    return (0)
}


#================================================
MAnorm <- function(sample1_rep1_signals,sample2_rep1_signals,table_MA)
{
  #--------Calculating M and A values for peaks-----------
  M<-log2((sample1_rep1_signals[,4]+1)/(sample2_rep1_signals[,4]+1))
  A<-0.5*log2((sample1_rep1_signals[,4]+1)*(sample2_rep1_signals[,4]+1))
  M <- as.matrix(M)
  A <- as.matrix(A)
  #--------linear regression-------------
  linear<-lm(M~A)$coefficients
  b<-rlm(M~A)$coefficients
  #--------------------------------------
  cat("M = b[1] + b[2] * A\n")
  log2_sample1_rep1_signals <- log2(sample1_rep1_signals[,4] + 1)
  log2_sample2_rep1_signals <- log2(sample2_rep1_signals[,4] + 1)
  log2_sample1_rep1_signals_rescaled <- (2-b[2])*log2_sample1_rep1_signals/(2+b[2]) - 2*b[1]/(2+b[2]);
  M_rescaled <- (log2_sample1_rep1_signals_rescaled - log2_sample2_rep1_signals);
  A_rescaled <- (log2_sample1_rep1_signals_rescaled + log2_sample2_rep1_signals)/2;
  #----------------------------------------------
  
  
  
  table_MA[,4] <- sample1_rep1_signals[,4]
  table_MA[,5] <- sample2_rep1_signals[,4]
  table_MA[,6] <- M_rescaled
  table_MA[,7] <- A_rescaled
  table_MA <-as.data.frame(table_MA)
  table_MA[,8] <- 0
  log2_sample1_rep1_signals_rescaled <- as.matrix(log2_sample1_rep1_signals_rescaled)
  #sample2_rep1_signals <- as.matrix(sample2_rep1_signals )
  for (n in c(1:nrow(table_MA))) 
  {
    print(n)
    #        cat(n,'\t',round(2^log2_peak_count_read1_rescaled[n]),'\t',peak_count_read2[n],'\n')
    table_MA[n,8]<--log10(pval(round(2^(log2_sample1_rep1_signals_rescaled[n])),round(sample2_rep1_signals[n,4])))
  }
  
  
  colnames(table_MA)[1] = "chr"
  colnames(table_MA)[2] = "start"
  colnames(table_MA)[3] = "end"
  colnames(table_MA)[4] = "#raw_read_1"
  colnames(table_MA)[5] = "#raw_read_2"
  colnames(table_MA)[6] = "M_value_rescaled"
  colnames(table_MA)[7] = "A_value_rescaled"
  colnames(table_MA)[8] = "-log10(p-value)"
  
  write.table(table_MA,"MAnorm_result.xls",sep="\t",quote=FALSE,row.names=FALSE)
  return(MAnorm_result)
}

#========================================
classify <- function(MAnorm_result,gene_expression,gene_coordinates)
{
  Results <- MAnorm_result
  gene_coordinates_21 <-subset(gene_coordinates,gene_coordinates[,2]==21)
  l=length(gene_coordinates_21[,1])
  for (i in 1:l) 
  {
    print(i)
    tss=gene_coordinates_21[i,3]
    if(gene_coordinates_21[i,5]==1)
    {
      upstream <- tss-8000
      downstream <- tss+2000
      gene_coordinates_21[i,9] <- upstream
      gene_coordinates_21[i,10] <- downstream
      gene_coordinates_21[i,11] <- seq(gene_coordinates_21[i,9],gene_coordinates_21[i,10],1)
    }
    else
    {
      upstream <- tss+8000
      downstream <- tss-2000
      gene_coordinates_21[i,9] <- upstream
      gene_coordinates_21[i,10] <- downstream
      gene_coordinates_21[i,11] <- seq(gene_coordinates_21[i,9],gene_coordinates_21[i,10],1)
    }
  }
  min_M=floor(min(MAnorm_result[,6]))
  max_M=ceil(max(MAnorm_result[,6]))
  num_groups=max_M-min_M
  #E116	GM12878 50
  #E123	K562 56
  l=length(Results[,1])
  
  for (i in 1:l) 
  {
    print(i)
    nearest_peak_index <- nearset_peak(MAnorm_result[i,2],gene_coordinates_21,gene_expression)
    if(nearest_peak_index!=0)
    {
      Results[i,9]<- gene_expression[nearest_peak_index,1]
      sample1_gene_expression <- gene_expression[nearest_peak_index,50]
      sample2_gene_expression <- gene_expression[nearest_peak_index,56]
      Results[i,10] <- sample1_gene_expression
      Results[i,11] <- sample2_gene_expression
    }
    else
    {
      Results[i,9] <- NA
      Results[i,10] <- NA
      Results[i,11] <- NA
    }
  }
  return(Results)
}
#=================================
gene_expression_analysis <- function(sample1_rep1_signals,gene_expression,gene_coordinates)
{
  rep1 <- sample1_rep1_signals
  rep1_score=matrix(nrow=rep1[length(rep1[,1]),3],ncol=1) #Scores for all genomic positions
  
  for(i in 1:length(rep1[,4]))
  {
    print(i)
    rep1_score[rep1[i,2]:rep1[i,3],1]<-rep1[i,4]
    
  }
  gene_coordinates_21 <-subset(gene_coordinates,gene_coordinates[,2]==21)
  
  #Sample 1 E116	GM12878 50  
  #Sample 2 E123	K562 56
  l=length(gene_coordinates_21[,1])
  
  for (i in 1:l) 
  {
    print(i)
    tss=gene_coordinates_21[i,3]
    if(gene_coordinates_21[i,5]==1)
    {
      upstream <- tss-8000
      downstream <- tss+2000
      gene_coordinates_21[i,9] <- upstream
      gene_coordinates_21[i,10] <- downstream
      peak_interval <- seq(gene_coordinates_21[i,9],gene_coordinates_21[i,10],1)
    }
    else
    {
      upstream <- tss+8000
      downstream <- tss-2000
      gene_coordinates_21[i,9] <- upstream
      gene_coordinates_21[i,10] <- downstream
      peak_interval <- seq(gene_coordinates_21[i,10],gene_coordinates_21[i,9],1)
    }
    
    if(min(peak_interval)>length(rep1_score))
    {
      gene_coordinates_21[i,11] <- NA
    }
    else
    {
      gene_coordinates_21[i,11] <- sum(rep1_score[peak_interval,1])#sum of the signals in the window size of 10000(from upstream to downstream)
    }
    
    gene_expression_index <- which(gene_expression[,1]==gene_coordinates_21[i,1])
    if(length(gene_expression_index)!=0)
    {
      gene_coordinates_21[i,12] <- gene_expression[gene_expression_index,50] # amount of the gene expression for the corresponding gene
    }
}
  
 Results <- gene_coordinates_21
 return(Results)
}

#==================================

#===============================
pval <- function(x, y)
{
  if (x+y<20) { # x + y is small
    p1<- nchoosek(x+y,x) * 2^-(x+y+1);
    p2<- nchoosek(x+y,y) * 2^-(x+y+1);
  }
  else { # if x+y is large, use approximation
    log_p1 <- (x+y)*log(x+y) - x*log(x) - y*log(y) - (x+y+1)*log(2);
    p1<-exp(log_p1);
    log_p2 <- (x+y)*log(x+y) - x*log(x) - y*log(y) - (x+y+1)*log(2);
    p2<- exp(log_p2);
  }
  pvalue=max(p1,p2)
  return(pvalue)
}
#===================================
MAnorm_raw_read <- function(rep1,rep2)
{
  
  save(mean_var,file="mean_var.Rdata")
}
MAnorm_var_stab <- function(rep1,rep2)
{
  
}

MAnorm_fold_enrich <- function(rep1,rep2)
{
  M<-log2((rep1[,4]+1)/(rep2[,4]+1))
  A<-0.5*log2((rep1[,4]+1)*(rep2[,4]+1))
  M <- as.matrix(M)
  A <- as.matrix(A)
}

Diff_var_stab <- function(rep1,rep2)
{
  l=length(rep1[,1])
  diff_signals <- matrix(nrow=l,ncol=4)
  for(i in 1:l)
  {
    print(i)
    diff_signals[i,1] <- rep1[i,1]
    diff_signals[i,2] <- rep1[i,2]
    diff_signals[i,3] <- rep1[i,3]
    diff_signals[i,4] <- rep2[i,4]-rep1[i,4]
  }
  diff_signals
}




#--------------------------------
target_gene<- function(start,end,gene_coordinates,gene_expression)
{
  summit <- ((end-start)/2)+start
  nearest_genes<-subset(gene_coordinates,(summit>=gene_coordinates[,9] && summit<=gene_coordinates[,10]))
  if (length(nearest_genes[,1])!=0)
  {
    tss <- gene_coordinates_21[,3]
    nearest_genes[,11]<- abs(tss-(nearest_genes[,3]))
    min_index <- which.min(nearest_genes[,9])
    nearest_gene <- nearest_genes[min_index,]
    nearest_gene<- as.matrix(nearest_gene)
    gene_expression_index <- which(gene_expression[,1]==nearest_gene[1,1])
    return (gene_expression_index)
  }
  else
    return (0)
}
