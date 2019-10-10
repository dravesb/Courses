#-------------------------------------
#
#       Bayesian HW3 
#
#-------------------------------------

setwd("~/Documents/Work/github/Courses/MA 578/hw3/")
set.seed(1985)

#------------------
# BDA 3.8 
#------------------

#  c
bikes_y <- c(16, 9, 10, 13, 19,20, 18, 17, 35, 55)
cars_y <- c(58, 90, 48, 57, 103, 57, 86, 112, 273, 64)
n <- bikes_y + cars_y

bikes_z <- c(112, 1, 2, 4, 9,7, 9, 8)
cars_z <- c(113, 18, 14, 44, 208, 67, 29, 154)
m <- bikes_z + cars_z

sim_post_y <- rbeta(1000, 1 + sum(bikes_y), 1 + sum(n - bikes_y))
sim_post_z <- rbeta(1000, 1 + sum(bikes_z), 1 + sum(m - bikes_z))

jpeg('post_dens.jpg')
plot(density(sim_post_y), col = "red", main = "3.8: Posterior Densities")
points(density(sim_post_z), col = "blue", type = "l")
dev.off()

diff_means <- numeric(1000)
for(i in 1:1000){
  
  #sample binomial
  
  sim_post_y <- rbeta(1000, 1 + sum(bikes_y), 1 + sum(n - bikes_y))
  sim_post_z <- rbeta(1000, 1 + sum(bikes_z), 1 + sum(m - bikes_z))
  
  diff_means[i] <- mean(sim_post_y) - mean(sim_post_z)
}

jpeg('diff_means.jpg')
hist(sim_post_y - sim_post_z,  main = "3.8: Difference in Means")
dev.off()

#------------------
# BDA 3.12 
#------------------

# e
#passenger_deaths <- c(734, 516,754,877,814,362,764,809,223,1066)
fatal_accidents <- c(24 , 25, 31, 31, 22, 21, 26, 20, 16, 22)
years <- 1976:1985

#crude linear models estimats 
lm_coef <- coef(lm(fatal_accidents ~ years))
alpha_lm<- unname(lm_coef[1])
beta_lm<- unname(lm_coef[2])

#crude Poisson regression model estimates
pois_coef <- coef(glm(fatal_accidents ~ years, family = "poisson"))
alpha_pois<- unname(pois_coef[1])
beta_pois <- unname(pois_coef[2])


















