#---------------------------------
#
#     Bayesian HW
#
#---------------------------------

#2.10 a 

summand <- function(n) 1/n * (.99)^n
C <- 1/sum(sapply(203:2000, summand))
C


#2.10 b 

summand <- function(n) (n - 279.5)^2 * 1/n * (.99)^n
sd <- sqrt(C * sum(sapply(203:2000, summand)))
sd


#2.13 a 
set.seed(1985)

alpha <- 3
beta <- 0.5
sy <- 238 
n <- 10 

quantile(rnbinom(1000, sy + alpha, 1 - 1/(1 + beta + n)), prob = c(.025, .975))


#2.13 b

#get sum of xi
passenger_deaths <- c(734, 516,754,877,814,362,764,809,223,1066)
death_rate <- c(.19,0.12,0.15,0.16,0.14,0.06,0.13,0.13,0.03,0.15)
x <- passenger_deaths/death_rate
sx <- sum(x)


set.seed(1985)

alpha <- 3
beta <- .5
sy <- 238 
n <- 10 
xtil <- 8000

quantile(rnbinom(1000, sy + alpha, 1 - xtil/(xtil + beta + sx)), prob = c(.025, .975))



#2.13 c
set.seed(1985)

alpha <- 3
beta <- 0.5
sy <- sum(passenger_deaths)
n <- length(passenger_deaths)

quantile(rnbinom(1000, sy + alpha, 1 - 1/(1 + beta + n)), prob = c(.025, .975))


#2.13 d

#get sum of xi

sx <- sum(x)
set.seed(1985)
alpha <- 3
beta <- .5 
xtil <- 8000

quantile(rnbinom(1000, sy + alpha, 1 - xtil/(xtil + beta + sx)), prob = c(.025, .975))


fatal_accidents <- c(24 , 25, 31, 31, 22, 21, 26, 20, 16, 22)





