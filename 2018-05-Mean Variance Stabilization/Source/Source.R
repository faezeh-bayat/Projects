Weighted_Mean_Variance<-function(rep1,rep2,rep1_score,rep2_score,distance,bin,alpha)
{
  scores=matrix(nrow=length(rep1_score),ncol=2)
  scores[,1]  <- rep1_score[,1]
  scores[,2]  <- rep2_score[,1]
  x=scores[,1]
  y=scores[,2]
  
  ordered_scores <- scores[order(x,y),]
  L=ceiling(length(ordered_scores[,1])/bin) # number of the bins
  l=bin
  Y=matrix(nrow=L,ncol = 1)
  a=matrix(nrow=L,ncol = 1)
  z=matrix(nrow=((2*distance)+1),1)
  n=(2*distance)+1
  mean_var=matrix(nrow=L,ncol=2)
  #########################
  for(i in 1:(L-1))
  {
    print (i)
    Y[i,1]=sum(ordered_scores[(((i-1)*l)+1):(i*l),2]) # sum of points in each bin
    a[i,1]=sum((ordered_scores[(((i-1)*l)+1):(i*l),2])^2) # sum-square of points in each bin
  }
  Y[L,1]=sum(ordered_scores[(((L-1)*l)+1):length(ordered_scores[,1]),2])
  a[L,1]=sum((ordered_scores[(((L-1)*l)+1):length(ordered_scores[,1]),2])^2)
  ##########################
  for (j in 1:n)
  {
    z[j,1] <- l*(1/(alpha^(abs(distance+1-j)*l)))
  }
  ##########################
  for(i in 1:L)
  {
    print(i)
    output=Weighted_Mean_Var(Y,z,a,i,distance)
    mean_var[i,1]=output[1]
    mean_var[i,2]=output[2]
  }
  mean_var
}
############################################
Weighted_Mean_Var<- function(weighted_sum,sum_of_weights,sum_square,center_bin,width)
{
  y=weighted_sum
  W=weighted_sum
  L=length(y)
  z=sum_of_weights
  k=center_bin
  w=width
  a=sum_square
  A=sum_square
  
  lower=max(1,(k-w))
  upper=min(L,(k+w))
  xx=upper-lower+1
  center_z=ceiling(length(z)/2)
  lower_z=(center_z-(k-lower))
  counter=lower_z
  upper_z=(center_z+(upper-k))
  for(i in lower:upper)
  {
    y[i,1]=y[i,1]*z[counter,1]
    counter=counter+1
  }
  weighted_mean=sum(y[lower:upper])/sum(z[lower_z:upper_z])
  
  counter=lower_z
  for(i in lower:upper)
  {
    a[i,1]=a[i,1]*z[counter,1]
    counter=counter+1
  }
  weighted_var=(sum(a[lower:upper])/sum(z[lower_z:upper_z]))-(weighted_mean^2) #some negative points
  
  meanvar_values <- c(weighted_mean,weighted_var)
  meanvar_values
}

#########################################################
polynomial_function_score <- function(mean_variance,rep1)
{
  replicate1_scores <- rep1
  ds <- data.frame(mean_variance)
  ds <-subset(ds,ds[,1]!="NA")
  x <- ds[,1]
  y <- ds[,2]
  model <- lm(y ~ poly(x,2))
  ss <- summary(model)$coefficients
  a0 <- ss[1,1]
  a1 <- ss[2,1]
  #a2 <- ss[3,1]
  #a3 <- ss[4,1]
  xx <- apply(rep1[,4,drop=F],1,polynomial_integralfunction,a0,a1)
  #xx <- apply(rep1[,4,drop=F],1,polynomial_integralfunction,a0,a1)
  xx <- data.frame(xx)
  replicate1_scores[,5] <- xx[,1]
  save(replicate1_scores,file="replicate1_scores.Rdata")
  replicate1_scores
}
################################

step_function_score <- function(mean_variance,rep1)
{
  data<- mean_variance
  
  x=data[,1]
  y=data[,2]
  
  ordered_data <- data[order(x,y),]
  ordered_data <-subset(ordered_data,ordered_data[,1]!="NA")
  x=ordered_data[,1]
  y=ordered_data[,2]
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
#########################################################
linear_function_score <- function(mean_variance,rep1)
{
  replicate1_scores <- rep1
  ds <- data.frame(mean_variance)
  x <- ds[,1]
  y <- ds[,2]
  #z <- nls(y ~ a * x +b, data = ds, start = list(a=1, b=1))
  z <- nls(y ~ a*x+b, data = ds, start = list(a=1, b=1))
  #z <- nls(y ~ a*x, data = ds, start = list(a=1))
  ss <- summary(z)$parameters
  a <- ss[1,1]
  b <- ss[2,1]
  xx <- apply(rep1[,4,drop=F],1,linear_integralfunction,a,b)
  #xx <- apply(rep1[,4,drop=F],1,linear_integralfunction,a)
  xx <- data.frame(xx)
  replicate1_scores[,5] <- xx[,1]
  #save(replicate1_scores,file="replicate1_scores.Rdata")
  replicate1_scores
}
#########################################################
curve_function_score <- function(mean_variance,rep1)
{
  replicate1_scores <- rep1
  ds <- data.frame(mean_variance)
  x <- ds[,1]
  y <- ds[,2]
  z <- nls(y ~ a * x^b, data = ds, start = list(a=1, b=1))
  ss <- summary(z)$parameters
  a <- ss[1,1]
  b <- ss[2,1]
  xx <- apply(rep1[,4,drop=F],1,curve_integralfunction,a,b)
  xx <- data.frame(xx)
  replicate1_scores[,5] <- xx[,1]
  save(replicate1_scores,file="replicate1_scores.Rdata")
  replicate1_scores
}
#########################################################
#########################################################
rep1_meanvar_score<-function(mean_variance,rep1)
{
  
  replicate1_scores <- rep1
  ds <- data.frame(mean_variance)
  x <- ds[,1]
  y <- ds[,2]
  z <- nls(y ~ a * x^b, data = ds, start = list(a=1, b=1))
  ss <- summary(z)$parameters
  a <- ss[1,1]
  b <- ss[2,1]
  xx <- apply(rep1[,4,drop=F],1,integralfunction,a,b)
  xx <- data.frame(xx)
  replicate1_scores[,5] <- xx[,1]
  save(replicate1_scores,file="replicate1_scores.Rdata")
  replicate1_scores
}
#############################################################
rep2_meanvar_score<-function(mean_variance,rep2)
{
  replicate2_scores <- rep2
  ds <- data.frame(mean_variance)
  x <- ds[,1]
  y <- ds[,2]
  z <- nls(y ~ a * x^b, data = ds, start = list(a=1, b=1))
  ss <- summary(z)$parameters
  a <- ss[1,1]
  b <- ss[2,1]
  xx <- apply(rep2[,4,drop=F],1,integralfunction,a,b)
  xx <- data.frame(xx)
  replicate2_scores[,5] <- xx[,1]
  save(replicate2_scores,file="replicate2_scores.Rdata")
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
#================================================
#========================================
classify_MAnorm <- function(MAnorm_result)
{
  
  #E116	GM12878 50
  #E123	K562 56
  Results <-subset(MAnorm_result,MAnorm_result[,12]!="NA")
  l=length(Results[,1])
  
  
  data=matrix(nrow=l,ncol = 2)
  
  
  min_M=(min(Results[,11]))
  max_M=(max(Results[,11]))
  ll=max_M-min_M
  a=min_M+(ll/4)
  b=a+ll/4
  c=b+ll/4
  
  for (i in 1:l)
  {
    print(i)
    if(Results[i,11]<=a)
    {
      Results[i,14] <- 1 
      data[i,1] <- "group1"
      if (Results[i,12]>Results[i,13])
      {
        data[i,2] <- "GN12878 highly expressed"
        Results[i,15] <- 1
        Results[i,16] <- -1
      }
      else
      {
        data[i,2] <- "K562 highly expressed"
        Results[i,15] <- -1
        Results[i,16] <- 1
      }
    }
    else if(Results[i,11]>a & Results[i,11]<=b)
    {
      Results[i,14] <- 2
      data[i,1] <- "group2"
      if (Results[i,12]>Results[i,13])
      {
        data[i,2] <- "GN12878 highly expressed"
        Results[i,15] <- 1
        Results[i,16] <- -1
      }
      else
      {
        data[i,2] <- "K562 highly expressed"
        Results[i,15] <- -1
        Results[i,16] <- 1
      }
    }
    else if(Results[i,11]>b & Results[i,11]<=c)
    {
      Results[i,14] <- 3
      data[i,1] <- "group3"
      if (Results[i,12]>Results[i,13])
      {
        data[i,2] <- "GN12878 highly expressed"
        Results[i,15] <- 1
        Results[i,16] <- -1
      }
      else
      {
        data[i,2] <- "K562 highly expressed"
        Results[i,15] <- -1
        Results[i,16] <- 1
      }
    }
    else if(Results[i,11]>c)
    {
      Results[i,14] <- 4
      data[i,1] <- "group4"
      if (Results[i,12]>Results[i,13])
      {
        data[i,2] <- "GN12878 highly expressed"
        Results[i,15] <- 1
        Results[i,16] <- -1
      }
      else
      {
        data[i,2] <- "K562 highly expressed"
        Results[i,15] <- -1
        Results[i,16] <- 1
      }
    }
  }
  
  
  
  
  data <- data.frame(groups=data[,1],celltype=data[,2])
  colnames(Results)[14] <- c("group")
  colnames(Results)[15] <- c("sample1 highly expressed")
  colnames(Results)[16] <- c("sample2 highly expressed")
  
  chart=matrix(nrow=8,ncol = 3)
  chart[1,1] <- "group1"
  chart[2,1] <- "group1"
  chart[3,1] <- "group2"
  chart[4,1] <- "group2"
  chart[5,1] <- "group3"
  chart[6,1] <- "group3"
  chart[7,1] <- "group4"
  chart[8,1] <- "group4"
  
  for (i in 1:4)
  {
    group <- subset(Results,Results[,14]==i)
    l=length(group[,1])
    sample1_high=subset(Results,(Results[,14]==i & Results[,15]==1))
    sample2_high=subset(Results,(Results[,14]==i & Results[,16]==1))
    chart[2*i-1,2] <- "GN12878 highly expressed"
    chart[2*i,2] <- "K562 highly expressed"
    
    chart[2*i-1,3] <- length(sample1_high[,1])
    chart[2*i,3] <- length(sample2_high[,1])
    group <- NULL
  }
  ggplot(data,aes(groups,fill=celltype))+geom_bar(position = "dodge2")
  write.table(Results,"Classify_Results.xls",sep="\t",quote=FALSE,row.names=FALSE)
  #return(chart)
  return(data)
}
#=================================
gene_expression_analysis <- function(sample1_rep1_signals,gene_expression,gene_coordinates)
{
  rep1 <- sample1_rep1_signals
  rep1_score=matrix(nrow=rep1[length(rep1[,1]),3],ncol=1) #Scores for all genomic positions
  gene_names <- as.matrix(gene_expression)
  
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
      gene_coordinates_21[i,11] <- mean(rep1_score[peak_interval,1])#sum of the signals in the window size of 10000(from upstream to downstream)
    }
    
    gene_expression_index <- which(gene_names[,1]==gene_coordinates_21[i,1])
    if(length(gene_expression_index)!=0)
    {
      gene_coordinates_21[i,12] <- asinh(gene_expression[gene_expression_index,50]) # amount of the gene expression for the corresponding gene
    }
}
  
 Results <- gene_coordinates_21
 return(Results)
}

#==================================
differential_expression_analysis <- function(sample1_rep1_signals,sample2_rep1_signals,table_MA,gene_expression,gene_coordinates)
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
  table_MA <-as.data.frame(table_MA)
  
  M_values=matrix(nrow=table_MA[length(table_MA[,1]),3],ncol=1) #Scores for all genomic positions
  gene_names <- as.matrix(gene_expression)
  
  for(i in 1:length(table_MA[,6]))
  {
    print(i)
    M_values[table_MA[i,2]:table_MA[i,3],1]<-table_MA[i,6]
    
  }
  gene_coordinates_21 <-subset(gene_coordinates,gene_coordinates[,2]==21)
  
  #Sample 1 E116	GM12878 50  
  #Sample 2 E123	K562 56
  l=length(gene_coordinates_21[,1])
  sample1_counter <- 0
  sample2_counter <- 0
  
  
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
    
    if(min(peak_interval)>length(M_values))
    {
      gene_coordinates_21[i,11] <- NA
    }
    else
    {
      gene_coordinates_21[i,11] <- mean(M_values[peak_interval,1])#Average of the signals in the window size of 10000(from upstream to downstream)
    }
    
    gene_expression_index <- which(gene_names[,1]==gene_coordinates_21[i,1])
    if(length(gene_expression_index)!=0)
    {
      gene_coordinates_21[i,12] <- gene_expression[gene_expression_index,50]
      gene_coordinates_21[i,13] <- gene_expression[gene_expression_index,56]
    }
  }
  
  colnames(gene_coordinates_21)[9] <- c("upstream")
  colnames(gene_coordinates_21)[10] <- c("downstream")
  colnames(gene_coordinates_21)[11] <- c("mean M_values")
  colnames(gene_coordinates_21)[12] <- c("sample1 gene expression")
  colnames(gene_coordinates_21)[13] <- c("sample2 gene expression")
  
  Results <- gene_coordinates_21
  write.table(Results,"Results.xls",sep="\t",quote=FALSE,row.names=FALSE)
  return(Results)
}
#===============================
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
ROC_AUC <- function(MAnorm_result)
{
  temp <-subset(MAnorm_result,MAnorm_result[,12]!="NA")
  l=length(temp[,1])
  data <- temp
  colnames(data)[10] <- c("Class")
  for(i in 1:l)
  {
    if(temp[i,12]>temp[i,13])
      data[i,10] <- "high"
    else
      data[i,10] <- "low"
  }
  data[,1:9] <- NULL
  data[,3:4] <- NULL
  data$Class <- as.factor(data$Class)
  smp_size <- floor(0.75 * nrow(data))
  train_ind <- sample(seq_len(nrow(data)), size = smp_size)
  trainData <- data[train_ind, ]
  testData <- data[-train_ind, ]
  trainX <- trainData
  trainX[,1] <- NULL
  testX <- testData
  testX[,1] <- NULL
  # ctrl <- trainControl(method = "repeatedcv", number = 5,summaryFunction=twoClassSummary,classProbs=TRUE,	allowParallel = TRUE)					
  # #grid <- expand.grid(interaction.depth=c(1,2), n.trees=c(10,20),	shrinkage=c(0.01,0.1),n.minobsinnode = 20)		
  # gbm.tune <- train(x=trainX,y=trainData$Class,method = "gbm",metric = "ROC",trControl = ctrl,verbose=FALSE)
  # gbm.pred <- predict(gbm.tune,testX)
  # confusionMatrix(gbm.pred,testData$Class)   
  #gbm.probs <- predict(gbm.tune,testX,type="prob")
  #head(gbm.probs)
  #gbm.ROC <- roc(predictor=gbm.probs$high,response=testData$Class,levels=rev(levels(testData$Class)))
  gbm.ROC <- roc(predictor=data[,1],response=data[,2],levels=rev(levels(testData$Class)))
  gbm.ROC$auc
  plot(gbm.ROC,main="GBM ROC")
  
  # Plot the propability of poor segmentation
  #histogram(~gbm.probs$PS|testData$Class,xlab="Probability of Poor Segmentation")
  return(gbm.ROC$auc)
}
#-------------------------------------


trapzfunction <- function(x,a,b)
{
  mean_var_curve_formula <- function(x) a * x^b
  s <- trapzfun(mean_var_curve_formula, 0,x)
  s$value
}
########################################

########################################
curve_integralfunction <- function(x,a,b)
{
  curve_function <- function(x) a * x^b
  s <- integral(curve_function, 0,x)
  return(s)
}
########################################
polynomial_integralfunction <- function(x,a0,a1)
{
  polynomial_function <- function(x) a0 + a1*x
  #polynomial_function <- function(x) a0 + a1*x + a2 *(x^2)
  #polynomial_function <- function(x) a0 + a1*x + a2 *(x^2) + a3 *(x^3)
  s <- integral(polynomial_function, 0,x)
  return(s)
}
##########################################
linear_integralfunction <- function(x,a,b)
{
  linear_function <- function(x) a * x +b
  s <- integral(linear_function, 0,x)
  return(s)
}

# linear_integralfunction <- function(x,a)
# {
#   linear_function <- function(x) a * x
#   s <- integral(linear_function, 0,x)
#   return(s)
# }

#####################################
# poly_integralfunction <- function(x,a0,a1)
# {
#   poly_function <- function(x) a0 + a1*x
#   s <- integral(poly_function, x1,x2)
#   return(s)
# }
# poly_integral<- function(mean_variance,rep1)
# {
#   replicate1_scores <- rep1
#   ds <- data.frame(mean_variance)
#   ds <-subset(ds,ds[,1]!="NA")
#   x <- ds[,1]
#   y <- ds[,2]
#   model <- lm(y ~ poly(x,2))
#   ss <- summary(model)$coefficients
#   a0 <- ss[1,1]
#   a1 <- ss[2,1]
#   ordered_scores <- replicate1_scores[order(replicate1_scores[,4]),]
#   l=length(ordered_scores[,1])
#   
#   replicate1_scores[1,5] <- poly_integralfunction(replicate1_scores[1,4],a0,a1)
#   s <- replicate1_scores[1,5]
#   for (i in 2:l)
#   {
#     
#     replicate1_scores[i,5] <- s+poly_integralfunction(replicate1_scores[i,4],a0,a1)
#   }
#   xx <- apply(rep1[,4,drop=F],1,polynomial_integralfunction,a0,a1)
#   xx <- data.frame(xx)
#   replicate1_scores[,5] <- xx[,1]
#   save(replicate1_scores,file="replicate1_scores.Rdata")
#   replicate1_scores
# }

step2_function_score <- function(mean_variance,rep1)
{
  data<- mean_variance
  x=data[,1]
  y=data[,2]
  ordered_data <- data[order(x,y),]
  replicate1_scores <- rep1
  ordered_scores <- replicate1_scores[order(replicate1_scores[,4]),]
  #base_integeal<-trapz(ordered_data[,1],ordered_data[,2])
  replicate1_scores <- ordered_scores
  
  # min_mean=min(x)
  # max_mean=max(x)
  # min_var=min(y)
  # max_var=max(y)
  a <- c(0,ordered_data[1,1])
  b <- c(0,ordered_data[1,2])
  replicate1_scores[1,5] <- trapz(a,b)
  s <- replicate1_scores[1,5]
  l=length(rep1[,4])
  for (i in 2:l)
  {
    print(i)
    a <- c(ordered_data[i-1,1],ordered_data[i,1])
    b <- c(ordered_data[i-1,2],ordered_data[i,2])
    replicate1_scores[i,5] <- s + trapz(a,b)
    s <- replicate1_scores[i,5] 
    
    # if(rep1[i,4]==0)
    #   replicate1_scores[i,5]<- 0
    # 
    # else if(rep1[i,4]<=max_mean)
    # {
    #   sub <- subset(ordered_data, ordered_data[,1]<=rep1[i,4])
    #   replicate1_scores[i,5] <- trapz(sub[,1],sub[,2])
    # }
    # else if(rep1[i,4]>max_mean)
    # {
    #   replicate1_scores[i,5] <- base_integeal+(max_var*(rep1[i,4]-max_mean))
    # }
    
  }
  replicate1_scores
}


###############################
step<- function(mean_variance,rep1)
{
  data<- mean_variance
  x=data[,1]
  y=data[,2]
  
  ordered_data <- data[order(x,y),]
  ordered_data <-subset(ordered_data,ordered_data[,1]!="NA")
  x=ordered_data[,1]
  y=ordered_data[,2]
  #base_integeal<-trapz(ordered_data[,1],ordered_data[,2])
  
  replicate1_scores <- rep1
  #min_mean=min(x)
  #max_mean=max(x)
  #min_var=min(y)
  #max_var=max(y)
  
  l=length(rep1[,4])
  for (i in 1:l)
  {
    print(i)
    
    if(rep1[i,4]==0)
      replicate1_scores[i,5]<- 0
    
    else
    {
      sub <- subset(ordered_data, ordered_data[,1]<=rep1[i,4] & ordered_data[,1]>=rep1[i-1,4])
      replicate1_scores[i,5] <- replicate1_scores[i-1,5]+trapz(sub[,1],sub[,2])
    }
    
  }
  replicate1_scores
}
