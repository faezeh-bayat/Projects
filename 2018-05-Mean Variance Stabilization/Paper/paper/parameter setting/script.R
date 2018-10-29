
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
setwd("/Users/faezehbayat/Documents/new experiments/gene expression analysis/mean-variance/H3K4me3 K562")
library(ggplot2)
data <- read.tcsv("results.csv")
ggplot(data) + aes(x = bin.factor, y = Pearson.Correlation) + geom_point(stat="identity") + theme_bw()+labs(x = "Bin size",y="Pearson correlation")
ggplot(data) + aes(x = bandwidth.factor, y = Pearson.Correlation) + geom_point(stat="identity") + theme_bw()+labs(x = "Band width",y="Pearson correlation")
ggplot(data) + aes(x = bandwidth.factor, shape = bin.factor, y = Pearson.Correlation) + geom_point(stat="identity",aes(shape=bin.factor, color=bin.factor)) + theme_bw()+labs(x = "Band width",y="Pearson correlation")
