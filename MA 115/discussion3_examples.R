#----------------------------------
#
#   Summary Statistics Lab
#
#----------------------------------

#sampling functions
samp_mix <- function(n, p, c1=-5, c2=5, sd1 = 1, sd2 = 1){
  samps <- numeric(n)
  for(i in 1:n) samps[i] <- ifelse(sample(1:2, 1, prob = c(p, 1-p)) == 1, rnorm(1, c1, sd1), rnorm(1, c2,sd2))
  return(samps)
} 
samp_box <- function(n,r) runif(n,-r,r)

#set parameters
n <- 20

#get samples
samp1 <- rnorm(n, 0, 1)
samp2 <- rnorm(n, 0, 3)
samp3 <- samp_mix(n, .5)
samp4 <- samp_box(n, 5)

#----------------------------------
#
#   Exloratory Data Analysis
#
#----------------------------------

plot(density(samp1, kernel = "epanechnikov"), col = 1, 
     xlim = c(-10, 10), main = "")
lines(density(samp2, kernel = "epanechnikov"), col = 2)
lines(density(samp3, kernel = "epanechnikov"), col = 3)
lines(density(samp4, kernel = "epanechnikov"), col = 4)

#----------------------------------
#
#   Summary Statistics
#
#----------------------------------

#variance functions
var <- function(x) 1/(length(x) - 1) * sum((x - mean(x))^2)
sd <- function(x) sqrt(var(x))


df <- cbind(samp1,samp2,samp3, samp4)

colMeans(df) #means 
apply(df, 2, var) #variances
apply(df, 2, sd) #sd


#----------------------------------
#
#   Summary Statistics
#
#----------------------------------


tuition <- c(50980, 52500, 48560, 44990, 49580, 44032) # BU, BC, NE, Harvard, MIT, Emmerson

#summary
summary(tuition)

#standardize
ave <- mean(tuition)
s <- sd(tuition)
z_scores <- (tuition - ave)/s
























