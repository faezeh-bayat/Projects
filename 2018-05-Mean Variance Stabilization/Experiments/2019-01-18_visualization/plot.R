read.tcsv <- function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}


library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

data <- read.csv("../Differential expression analysis.csv")
#data <- read.csv("../Gene expression analysis.csv")
other_methods_data <- read.csv("../Other_methods_differential_expression.csv")
#data <- read.csv("../Other_methods_gene_expression.csv")
data$Model <- revalue(data$Model, c("Model 1"="Model_1", "Model 2"="Model_2", "Model 3"="Model_3"))
data <- filter(data, Model=="Model_1")

data$alpha.factor <- as.factor(data$Alpha)
data$bin.factor <- as.factor(data$Bin.size)
data$width.factor <- as.factor(data$Width)

# step.data <- data[data$Method == "Variance stabilization",]
# step.data <- step.data[step.data$Curve == "step",]
# step.data <- step.data[step.data$Bin.size == 10000,]
# step.data <- step.data[step.data$Alpha == 1,]
# step.data <- step.data[step.data$Width == 3,]
# other.data <- data[data$Method != "Variance stabilization",]
# best.data <- rbind(other.data, step.data)

mean_data <- summarize(group_by(data, Beta, Width, Bin.size, Model), avg_eval=mean(Evaluation.Value))


ot#ggplot(best.data) + aes(x = as.factor(Bin.size), y = Evaluation.Value) + geom_bar(stat="identity") + theme_bw()

#ggplot(x) + aes(x = as.factor(Bin.size), y = Evaluation.Value) + geom_bar(stat="identity") + theme_bw()

ggplot(mean_data) + 
  aes(x=as.factor(Beta), shape=as.factor(Bin.size), group=as.factor(Bin.size), y = avg_eval) + 
  geom_point() + geom_line() +
  facet_wrap(~Model) + 
  theme_bw()

#############
y <- dcast(x, Width + Bin.size + Beta ~ Model)
ggplot(y) + 
  aes(x=Model_1, y=Model_2) +
  geom_abline() +
  geom_point() + 
  geom_line() +
  theme_bw()

###################
other_methods_mean_data <- summarize(group_by(other_methods_data, Method), avg_eval=mean(Pearson.Correlation))
