setwd("/Users/faezehbayat/Documents/FINAL experiments(cedar)")
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/Total results.csv")
library(ggplot2)

##################Gene expression analysis visualization#######################

step.data <- data[data$Method == "SDPM",]
step.data <- step.data[step.data$Evaluation.type == "Gene expression analysis",]
#step.data <- step.data[step.data$Curve == "step",]
step.data <- step.data[step.data$Bin.size == 100000,]
#step.data <- step.data[step.data$Alpha == 1,]
step.data <- step.data[step.data$Beta == 1000000,]
#step.data <- step.data[step.data$Width == 3,]
step.data <- step.data[step.data$Model == "Model 1",]
other.data <- data[data$Method != "SDPM",]
other.data <- other.data[other.data$Evaluation.type == "Gene expression analysis",]
best.data <- rbind(other.data, step.data)
x=best.data$Dataset
y=best.data$Evaluation.Value
grp=best.data$Method
Methods=best.data$Method
p=ggplot(best.data, aes(x, y, group = Methods)) + geom_col(aes(fill = Methods), position = "dodge2") +theme_bw()+
  scale_fill_brewer(palette="Dark2")+labs(x = "\ncell types",y="\nLinearity score")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("newgeneexpression.pdf", p, width=15, height=8, units="cm",scale=1)

p
################   Mean of cell types in gene expression #########
step.data <- data[data$Method == "SDPM",]
step.data <- step.data[step.data$Evaluation.type == "Gene expression analysis",]
#step.data <- step.data[step.data$Curve == "step",]
step.data <- step.data[step.data$Bin.size == 100000,]
#step.data <- step.data[step.data$Alpha == 1,]
step.data <- step.data[step.data$Beta == 1000000,]
#step.data <- step.data[step.data$Width == 3,]
step.data <- step.data[step.data$Model == "Model 1",]
other.data <- data[data$Method != "SDPM",]
other.data <- other.data[other.data$Evaluation.type == "Gene expression analysis",]
best.data <- rbind(other.data, step.data)
mean.data=best.data
a1<- best.data[best.data$Method=="asinh(Fold enrichment)",]
mean_a1=mean(a1[,9])
a1[1,9]=mean_a1
a2<- best.data[best.data$Method=="asinh(P-value)",]
mean_a2=mean(a2[,9])
a2[1,9]=mean_a2
a3<- best.data[best.data$Method=="Fold enrichment",]
mean_a3=mean(a3[,9])
a3[1,9]=mean_a3
a4<- best.data[best.data$Method=="log(Fold enrichment +1)",]
mean_a4=mean(a4[,9])
a4[1,9]=mean_a4
a5<- best.data[best.data$Method=="P-value",]
mean_a5=mean(a5[,9])
a5[1,9]=mean_a5
a6<- best.data[best.data$Method=="log(P-value+1)",]
mean_a6=mean(a6[,9])
a6[1,9]=mean_a6
a7<- best.data[best.data$Method=="SDPM",]
mean_a7=mean(a7[,9])
a7[1,9]=mean_a7
new.data=rbind(a1[1,],a2[1,],a3[1,],a4[1,],a5[1,],a6[1,],a7[1,])
x=new.data$Method
y=new.data$Evaluation.Value

p=ggplot(new.data, aes(x, y)) + geom_bar(stat="identity",fill="steelblue")+theme_bw()+
  labs(x = "\nSignal transformation methods",y="\nAverage linearity score")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_text(aes(label=y), vjust=0,label=format(y, digits=2))

ggsave("meangeneexpression.pdf", p, width=15, height=11, units="cm",scale=1)
p
###################Differential expression analysis##########################

data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/Total results.csv")
step.data <- data[data$Method == "SDPM",]
step.data <- step.data[step.data$Evaluation.type == "Differential expression",]
#step.data <- step.data[step.data$Curve == "step",]
step.data <- step.data[step.data$Bin.size == 100000,]
#step.data <- step.data[step.data$Alpha == 1,]
step.data <- step.data[step.data$Beta == 1000000,]
#step.data <- step.data[step.data$Width == 3,]
step.data <- step.data[step.data$Model == "Model 1",]
other.data <- data[data$Method != "SDPM",]
other.data <- other.data[other.data$Evaluation.type == "Differential expression",]
best.data <- rbind(other.data, step.data)
x=best.data$Dataset
y=best.data$Evaluation.Value
grp=best.data$Method
Methods=best.data$Method
p=ggplot(best.data, aes(x, y, group = Methods)) + geom_col(aes(fill = Methods), position = "dodge2") +theme_bw()+
  scale_fill_brewer(palette="Dark2")+labs(x = "\nCell type pairs",y="\nDifferential expression score")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave("differentialexpression.pdf", p, width=15, height=10, units="cm",scale=1)

######################  Mean of cell types in differential expression##############################
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/Total results.csv")
step.data <- data[data$Method == "SDPM",]
step.data <- step.data[step.data$Evaluation.type == "Differential expression",]
#step.data <- step.data[step.data$Curve == "step",]
step.data <- step.data[step.data$Bin.size == 100000,]
#step.data <- step.data[step.data$Alpha == 1,]
step.data <- step.data[step.data$Beta == 1000000,]
#step.data <- step.data[step.data$Width == 3,]
step.data <- step.data[step.data$Model == "Model 1",]
other.data <- data[data$Method != "SDPM",]
other.data <- other.data[other.data$Evaluation.type == "Differential expression",]
best.data <- rbind(other.data, step.data)
mean.data=best.data
a1<- best.data[best.data$Method=="asinh(Fold enrichment)",]
mean_a1=mean(a1[,9])
a1[1,9]=mean_a1
a2<- best.data[best.data$Method=="asinh(P-value)",]
mean_a2=mean(a2[,9])
a2[1,9]=mean_a2
a3<- best.data[best.data$Method=="Fold enrichment",]
mean_a3=mean(a3[,9])
a3[1,9]=mean_a3
a4<- best.data[best.data$Method=="log(Fold enrichment +1)",]
mean_a4=mean(a4[,9])
a4[1,9]=mean_a4
a5<- best.data[best.data$Method=="P-value",]
mean_a5=mean(a5[,9])
a5[1,9]=mean_a5
a6<- best.data[best.data$Method=="log(P-value+1)",]
mean_a6=mean(a6[,9])
a6[1,9]=mean_a6
a7<- best.data[best.data$Method=="SDPM",]
mean_a7=mean(a7[,9])
a7[1,9]=mean_a7
new.data=rbind(a1[1,],a2[1,],a3[1,],a4[1,],a5[1,],a6[1,],a7[1,])
x=new.data$Method
y=new.data$Evaluation.Value

p=ggplot(new.data, aes(x, y)) + geom_bar(stat="identity",fill="steelblue")+theme_bw()+
  labs(x = "\nSignal transformation methods",y="\nAverage differential expression score")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_text(aes(label=y), vjust=0,label=format(y, digits=2))
p
ggsave("meandifferentialexpression.pdf", p, width=15, height=11, units="cm",scale=1)



##################Log likelihood#######################
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/log likelihood analysis.csv")
data[,9] <- -1*data[,9]
x=data$Dataset
y=data$Evaluation.Value
Methods=data$Evalation.type
p=ggplot(data, aes(x, y, group = Methods)) + geom_col(aes(fill = Methods), position = "dodge2") +theme_bw()+
  scale_fill_brewer(palette="Dark2")+labs(x = "\ncell types",y="\n-Log likelihood")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
#ggplot(data, aes(x, y, group = Methods)) + geom_col(aes(fill = Methods), position = "dodge2") +theme_bw()
ggsave("loglikelihood.pdf", p, width=15, height=7, units="cm",scale=1)

################## Mean of Log likelihood #############
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/log likelihood analysis.csv")
data[,9] <- -1*data[,9]
best.data <- data
a1<- best.data[best.data$Evalation.type=="SDPM",]
mean_a1=mean(a1[,9])
a1[1,9]=mean_a1
a2<- best.data[best.data$Evalation.type=="log(x+1)",]
mean_a2=mean(a2[,9])
a2[1,9]=mean_a2
a3<- best.data[best.data$Evalation.type=="asinh(x)",]
mean_a3=mean(a3[,9])
a3[1,9]=mean_a3
a4<- best.data[best.data$Evalation.type=="No transformation",]
mean_a4=mean(a4[,9])
a4[1,9]=mean_a4

new.data=rbind(a1[1,],a2[1,],a3[1,],a4[1,])
x=new.data$Evalation.type
y=new.data$Evaluation.Value

p=ggplot(new.data, aes(x, y)) + geom_bar(stat="identity",fill="steelblue")+theme_bw()+
  labs(x = "\nSignal transformation methods",y="\nAverage -log likelihood")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_text(aes(label=y), vjust=0,label=format(y, digits=2))
p
ggsave("meanloglikelihood.pdf", p, width=15, height=11, units="cm",scale=1)

######################################################parameter setting##################
setwd("/Users/faezehbayat/Documents/FINAL experiments(cedar)")
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/Total results.csv")
step.data <- data[data$Method == "SDPM",]
step.data <- step.data[step.data$Evaluation.type == "Gene expression analysis",]
data <- step.data
mean_data <- summarize(group_by(data, Beta, Width, Bin.size, Model), avg_eval=mean(Evaluation.Value))
step.data <- mean_data

#data$alpha.factor <- as.factor(data$Alpha)
#data$bin.factor <- as.factor(data$Bin.size)
#data$width.factor <- as.factor(data$Width)
#step.data <- data[data$Method == "SDPM",]
#y = Evaluation.Value
#p=ggplot(step.data) + aes(x = as.factor(Bin.size),y=avg_eval ) + geom_point(stat="identity") + theme_bw()+labs(x = "\nBin size",y="\nLinearity score")
#ggsave("newnewbinpearson.png", p, width=12, height=10, units="cm",scale=1)


#p=ggplot(step.data) + aes(x = as.factor(Beta), y = avg_eval) + geom_point(stat="identity") + theme_bw()+labs(x = "\nBand width",y="\nLinearity score")
#ggsave("newnewbandwidthpearson.png", p, width=12, height=10, units="cm",scale=1)


p=ggplot(step.data) + aes(x = as.factor(Beta), shape = as.factor(Bin.size), y = avg_eval) + geom_point(stat="identity",aes( color=as.factor(Bin.size))) + theme_bw()+labs(x = "\nBand width",y="\nLinearity score")
ggsave("expressionnewBinbandpearson.png", p, width=12, height=10, units="cm",scale=1)
##############
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/Total results.csv")
step.data <- data[data$Method == "SDPM",]
step.data <- step.data[step.data$Evaluation.type == "Differential expression",]
data <- step.data
mean_data <- summarize(group_by(data, Beta, Width, Bin.size, Model), avg_eval=mean(Evaluation.Value))
step.data <- mean_data
p=ggplot(step.data) + aes(x = as.factor(Beta), shape = as.factor(Bin.size), y = avg_eval) + geom_point(stat="identity",aes( color=as.factor(Bin.size))) + theme_bw()+labs(x = "\nBand width",y="\nLinearity score")
ggsave("diffnewBinbandpearson.png", p, width=12, height=10, units="cm",scale=1)
##################################
setwd("/Users/faezehbayat/Documents/FINAL experiments(cedar)")
load("mean_variance.Rdata")
xx=data.frame(mean_variance)
p=ggplot(xx, aes(xx[,1],sqrt(xx[,2]))) + geom_point(stat="identity")+theme_bw()+labs(x = "Mean",y=expression(sigma))
p
ggsave("newnewmeansigma2.png", p, width=12, height=10, units="cm",scale=1)
plot(xx[,1],1/sqrt(xx[,2]))
dd <- subset(xx,xx[,1]<=10 & (1/sqrt(xx[,2]))<=400000 )
p=ggplot(dd, aes(dd[,1],1/sqrt(dd[,2]))) + geom_point(stat="identity")+theme_bw()+labs(x = "Mean",y=expression(1/sigma))
p
ggsave("newnewmeansigma.png", p, width=12, height=10, units="cm",scale=1)


setwd("/Users/faezehbayat/Desktop")
load("1.Rdata")
replicate1_scores[,6] <- log(replicate1_scores[,4]+1)
replicate1_scores[,7] <- asinh(replicate1_scores[,4])
xx=data.frame(replicate1_scores)

xx=subset(xx,xx[,4]<=10)
p=ggplot(xx, aes(xx[,4],xx[,7])) + geom_line(stat="identity")+theme_bw()+labs(x = "x",y="f(x)")
p
ggsave("newk562xfx.png", p, width=12, height=10, units="cm",scale=1)

library(ggplot2)
library(ggrepel)
library(dplyr)
cyl <- c("SDPM","log(x+1)","asinh(x)")

p1=ggplot(xx, aes(xx[,4],xx[,7])) + geom_line(stat="identity")+theme_bw()+labs(x = "x",y="f(x)")
p1
p = ggplot()+ 
  geom_line(data = xx, aes(x = xx[,4], y = xx[,5]), color = "blue")+
  geom_line(data = xx, aes(x = xx[,4], y = xx[,6]), color = "red") +
  geom_line(data = xx, aes(x = xx[,4], y = xx[,7]), color = "purple") +
  xlab('x') +
  ylab('f(x)')+theme_bw()

p

p + scale_color_manual(name="Transformation methods", 
                        labels = c("SDPM", 
                                   "log(x+1)", 
                                   "asinh(x)"), 
                        values = c("SD"="blue", 
                                   "lo"="red", 
                                   "as"="purple"))
                                   
ggsave("xfx.png", p, width=12, height=10, units="cm",scale=1)
################################
setwd("/Users/faezehbayat/Documents/FINAL experiments(cedar)")
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/total_assay_log_results_analysis.csv")
library(ggplot2)
data[,9] <- -1*data[,9]
best.data <- data
a1<- best.data[best.data$Evalation.type=="SDPM",]
mean_a1=mean(a1[,9])
a1[1,9]=mean_a1
a2<- best.data[best.data$Evalation.type=="log(x+1)",]
mean_a2=mean(a2[,9])
a2[1,9]=mean_a2
a3<- best.data[best.data$Evalation.type=="asinh(x)",]
mean_a3=mean(a3[,9])
a3[1,9]=mean_a3
a4<- best.data[best.data$Evalation.type=="No transformation",]
mean_a4=mean(a4[,9])
a4[1,9]=mean_a4

new.data=rbind(a1[1,],a2[1,],a3[1,],a4[1,])
x=new.data$Evalation.type
y=new.data$Evaluation.Value

p=ggplot(new.data, aes(x, y)) + geom_bar(stat="identity",fill="steelblue")+theme_bw()+
  labs(x = "\nSignal transformation methods",y="\nAverage -log likelihood")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_text(aes(label=y), vjust=0,label=format(y, digits=2))
p
ggsave("all_assays_loglikelihood.pdf", p, width=15, height=11, units="cm",scale=1)


#########################
setwd("/Users/faezehbayat/Documents/FINAL experiments(cedar)")
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
data <- read.csv("/Users/faezehbayat/Documents/FINAL experiments(cedar)/assay_likelihood_analysis.csv")
step.data <- data[data$Evalation.type == "SDPM",]
step.data <- data[step.data$Dataset == "H3K36me3 H1-hESC",]
data <- step.data
mean_data <- summarize(group_by(data, Beta, Width, Bin.size), avg_eval=mean(Evaluation.Value))
step.data <- mean_data

#data$alpha.factor <- as.factor(data$Alpha)
#data$bin.factor <- as.factor(data$Bin.size)
#data$width.factor <- as.factor(data$Width)
#step.data <- data[data$Method == "SDPM",]
#y = Evaluation.Value
#p=ggplot(step.data) + aes(x = as.factor(Bin.size),y=avg_eval ) + geom_point(stat="identity") + theme_bw()+labs(x = "\nBin size",y="\nLinearity score")
#ggsave("newnewbinpearson.png", p, width=12, height=10, units="cm",scale=1)


#p=ggplot(step.data) + aes(x = as.factor(Beta), y = avg_eval) + geom_point(stat="identity") + theme_bw()+labs(x = "\nBand width",y="\nLinearity score")
#ggsave("newnewbandwidthpearson.png", p, width=12, height=10, units="cm",scale=1)


p=ggplot(step.data) + aes(x = as.factor(Beta), shape = as.factor(Bin.size), y = avg_eval) + geom_point(stat="identity",aes( color=as.factor(Bin.size))) + theme_bw()+labs(x = "\nBand width",y="\nLinearity score")
p
ggsave("parameter_settin_likelihood.pdf", p, width=12, height=10, units="cm",scale=1)
