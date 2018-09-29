library("airway")
data("airway")
library("DESeq2")
library("ggplot2")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library("qvalue")


data <- DESeqDataSet(airway, design=~0 + cell + dex)
clean_data <- data[rowSums(counts(data)) > 1, ]
analyzed_data <- DESeq(clean_data)

pval <- results(analyzed_data)$pvalue
ggplot(data.frame(pval = pval), aes(x = pval))+
  geom_histogram(fill = "coral", col = "blue")+
  labs(title = "P-Values")
  theme_bw()
  
#qvalue package vio bionconductor
qval <- qvalue(pval)
plot(qval)


# Exercise 2 

f0 <- function(x) 1/sqrt(2*pi)*exp(-(x^2)/2) 
f1 <- function(x)  1/2 * 1/sqrt(2*pi)*exp(-((x-2.5)^2)/2) + 1/2 * 1/sqrt(2*pi)*exp(-((x+2.5)^2)/2)
p0 <- 0.95
p1 <- 1-p0

x <- seq(from = -5, to = 5, by = .05)
plot(x, p0*f0(x), xlim = c(-5,5), ylim = c(0,.4),type = "l", col = "blue", main = "Non-Zero Assumptions")
points(x, p1*f1(x), xlim = c(-5,5), ylim = c(0,.4),type = "l", col = "red")
points(x, p0*f0(x) + p1*f1(x), xlim = c(-5,5), ylim = c(0,.4),type = "l", col = "green")

x <- seq(from = -1, to = 1, by = .05)
plot(x, p0*f0(x), xlim = c(-1,1), ylim = c(0,.4),type = "l", col = "blue", main = "Non-Zero Assumptions")
points(x, p1*f1(x), xlim = c(-1,1), ylim = c(0,.4),type = "l", col = "red")
points(x,p0*(f0(x) + p1/p0*f1(x)), xlim = c(-1,1), ylim = c(0,.4),type = "l", col = "green")

