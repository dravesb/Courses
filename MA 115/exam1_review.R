#-----------------------------------
#
#       Discussion Example
#
#-----------------------------------

#-------------------------
#   Set up data
#-------------------------

n <- 10
f <- function(x) rnorm(1, mean = 3*x, sd = 2)
x <- round(rnorm(n), 1)
y <- round(sapply(x,f),1)

#-----------------------------
#   Univariate Data Analysis
#-----------------------------

#print data
print(x)

#measures of center
mean(x) 
median(x)

#measures of dispersion
var(x) 
sd(x)
sd(x)/mean(x)

#5 number summary
summary(x)

#IQR
Q1 <- unname(summary(x)[2])
Q3 <- unname(summary(x)[5])
IQR <-Q3 - Q1

#Check for outliers
LIF <- Q1 - 1.5 * IQR
LOF <- Q1 - 3 * IQR

UIF <- Q3 + 1.5 * IQR
UOF <- Q3 + 3 * IQR

#plot these
plot(x, rep(0,n), xlim = c(LOF, UOF))
abline(v = Q1, lty = 2); abline(v = Q3, lty = 2)
abline(v = LIF, col = "red", lty = 2); abline(v = UIF, col = "red", lty = 2)
abline(v = LOF, col = "blue", lty = 2); abline(v = UOF, col = "blue", lty = 2)

#boxplots
boxplot(x)
abline(h = min(x[which(x > LIF)]))
abline(h = max(x[which(x < UIF)]))
abline(h = median(x))

#correlation
plot(x, y)

#correlation
cor(x, y)
cor(-x, y)
cor(-x, -y)

#regression coefficent
rho <- cor(x,y)
sdx <- sd(x)
sdy <- sd(y)


beta1 <- rho * sdy/sdx
beta0 <- mean(y) - beta1 * mean(x)

plot(x, y)
abline(a = beta0, b = beta1)




