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

data <- read.tcsv("HUVEC gene expression results.csv")

library(ggplot2)

data$alpha.factor <- as.factor(data$Alpha)
data$bin.factor <- as.factor(data$Bin.size)
data$width.factor <- as.factor(data$Width)

step.data <- data[data$Method == "Variance stabilization",]
step.data <- step.data[step.data$Curve == "step",]
step.data <- step.data[step.data$Bin.size == 10000,]
step.data <- step.data[step.data$Alpha == 1,]
step.data <- step.data[step.data$Width == 3,]


other.data <- data[data$Method != "Variance stabilization",]

best.data <- rbind(other.data, step.data)

ggplot(best.data) + aes(x = as.factor(Method), y = Evaluation.Value) + geom_bar(stat="identity") + theme_bw()