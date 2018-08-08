#Mouse data visualization----------------------------------------
setwd("Desktop/GLM_Final_Project/Mouse_Data")

#Define useful functions-------------------------------------------------
format_global = function(x){
  n = length(x)
  ind = combn(1:n, 2)
  X = matrix(NA, ncol = 2, nrow = n*(n-1)/2)
  for(i in 1:(n*(n-1)/2)){
    X[i,] = c(x[ind[1,i]], x[ind[2,i]])
  }
  return(X)
}
#Read in/Format data-------------------------------------------------

#response
brain = as.matrix(read.csv("N54859.csv", header = F))
n = nrow(brain)
colnames(brain) = rownames(brain) = 1:n
head(brain)

#covariates
covar = read.csv("first_page.csv")[,-c(1:2,9:11)]
covar = covar[-nrow(covar),]
for(i in 1:ncol(covar)){
  covar[,i] = as.character(droplevels(covar[,i]))
}
covar = as.matrix(covar)

#Left versus Right Brain network----------------------------------------
side = format_global(covar[,1])
indL = which(side[,1] == side[,2] &  side[,2] == "Left")
indr = which(side[,1] == side[,2] &  side[,2]== "Right")
indd = which(side[,1] != side[,2] )




