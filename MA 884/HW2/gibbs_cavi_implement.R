#helper functions 
library(limma); library(MASS)
LOGEPS <- log(.Machine$double.eps / 2)
log1pe <- function (x) { # vectorized version: `x` can be a vector
  l <- ifelse(x > 0, x, 0) # shift
  x <- ifelse(x > 0, -x, x) # range reduction: `x = -abs(x)`
  ifelse(x < LOGEPS, l, l + log(1 + exp(x)))
}

#-----------------------------------------------
#
#             Gibbs sampler 
#
#-----------------------------------------------
#Inputs
#   y - response vector
#   X - design matrix
#   hyper = c(q, r, gamma, sigma2) hyper parameter
#   theta - initial theta vector
#   delta - initial delta vector
#   max.iter - max iteration
gibbs <- function(y, X, hyper, theta, delta, max.iter = 10000){

  #set up initial parameters
  p <- ncol(X)
  n <- nrow(y)
  
  #unpack hyper parameters
  q <- hyper[1]
  r <- hyper[2]
  gamma <- hyper[3]
  sigma2 <- hyper[4]
  
  #set up frequently used values
  XtX <- crossprod(X)
  
  #set up parameter storage
  Delta <- matrix(NA, nrow = max.iter , ncol = p)
  Delta[1,] <- delta
  Theta <- matrix(NA, nrow = max.iter , ncol = p)
  Theta[1,] <- theta
  
  
  #begin Markov Chain
  for(t in 2:max.iter){
    if(t %% 10 == 0){print(t)}
    #sample deltas--------------
    for(j in 1:p){
      #set log[p] = log[a/(a+b)] = loga - log a p b 
      loga <- log(q*sqrt(r)) -.5 *log(2*pi) - r*.5 * Theta[t-1, j]^2
      logb <- log(1-q) -.5 *log(2*pi*gamma) - 1/(2*gamma) * Theta[t-1, j]^2
      logapb <- logsumexp(loga, logb)
      
      #get log p
      logp <- loga - logapb
      
      #sample delta
      Delta[t, j] <- ifelse(log(runif(1))<logp, 1, 0)
    }
    
    #sample thetas---------------
    
    #set up Dinverse
    Dinv <- diag(ifelse(Delta[t,] == 1, r, 1/gamma))
    
    #get (mu, Sigma_til) Sigma_til = Sigma / sigma^2 
    Sigma_til <- sigma2 * solve(XtX + sigma2 *Dinv) 
    mu <- Sigma_til %*% crossprod(X,y)
    
    #get sample
    Theta[t, ] <- mvrnorm(n=1, mu, Sigma = sigma2 * Sigma_til)
    
  }

  return(list(Delta, Theta))  
}

#-----------------------------------------------
#
#             Variational
#             Approximation
#
#-----------------------------------------------
#Inputs
#   y - response vector
#   X - design matrix
#   hyper = c(q, r, gamma, sigma2) hyper parameter
#   q - initial q vector
#   mu - initial mu vector
#   v2 - initial v vector
#   max.iter - max iteration

VA <- function(y, X, hyper, q, mu, v2, max.iter = 100){
  
  #set up initial parameters
  p <- ncol(X)
  n <- nrow(y)
  
  #set up parameter storage
  logQ <- matrix(NA, nrow = max.iter, ncol = p)
  logQ[1,] <- log(q)
  M <- matrix(NA, nrow = max.iter, ncol = p)
  M[1,] <- mu
  logV <- matrix(NA, nrow = max.iter , ncol = p)
  logV[1,] <- log(v2)
  
  #unpack hyper parameters
  q <- hyper[1]
  r <- hyper[2]
  gamma <- hyper[3]
  sigma2 <- hyper[4]
  
  #begin CAVI updates
  for(t in 2:max.iter){
    
    if(t %% 100 == 0){print(t)}
    
    #sample q--------------
    for(j in 1:p){
      #set log[p] = log[a/(a+b)] = loga - log a p b 
      loga <- log(q*sqrt(r)) -.5 *log(2*pi) - r*.5 * (M[t-1, j]^2 + exp(logV[t-1, j]))
      logb <- log(1-q) -.5 *log(2*pi*gamma) - 1/(2*gamma) * (M[t-1, j]^2 + exp(logV[t-1, j]))
      
      #Update Q
      logQ[t, j] <- loga - logsumexp(loga, logb)
    }
    
    #update (mu, v2) ---------------
    #copy over precious mean values for ease of comupation
    M[t, ] <- M[t-1,]
    for(j in 1:p){
      
      #helper variable
      den <- -2 * crossprod(X[,j]) + exp(logQ[t, j])*r  + (1 - exp(logQ[t, j]))/gamma
      #update log variance
      logV[t, j] <- -log(den)
      
      #update log mean
      M[t, j] <- (2* crossprod(X[,j], y - X[,-j] %*% M[t,-j]) )/den
    }
    
    
  }
  
  return(list(logQ, M, logV))
}

#-----------------------------------------------
#
#             Simulation Study
#
#-----------------------------------------------

set.seed(330)
p <- 500
n <- p/10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
coef <- matrix(c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15), nrow = 10)
Y <- X[,1:10] %*% coef + rnorm(50)
hyper <- c(1/p^2, 1/sqrt(n), .1/max(eigen(crossprod(X))$values), sigma2 = 1)
  
gibbs_results1 <- gibbs(Y, X, hyper, rnorm(p), rbinom(p, 1, .25))  
variational_results1 <- VA(Y, X, hyper, rbeta(p, 2,3), rnorm(p), rgamma(p, 3))  
  
set.seed(331)
p <- 1000
n <- p/10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
coef <- matrix(c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15), nrow = 10)
Y <- X[,1:10] %*% coef + rnorm(50)
hyper <- c(1/p^2, 1/sqrt(n), .1/max(eigen(crossprod(X))$values), sigma2 = 1)

gibbs_results2 <- gibbs(Y, X, hyper, rnorm(p), rbinom(p, 1, .25))  
variational_results2 <- VA(Y, X, hyper, rbeta(p, 2,3), rnorm(p), rgamma(p, 3))  

set.seed(332)
p <- 2000
n <- p/10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
coef <- matrix(c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15), nrow = 10)
Y <- X[,1:10] %*% coef + rnorm(50)
hyper <- c(1/p^2, 1/sqrt(n), .1/max(eigen(crossprod(X))$values), sigma2 = 1)

gibbs_results3 <- gibbs(Y, X, hyper, rnorm(p), rbinom(p, 1, .25))  
variational_results3 <- VA(Y, X, hyper, rbeta(p, 2,3), rnorm(p), rgamma(p, 3))  



#write out Gibbs data
Delta1 <- gibbs_results1[[1]]
Theta1 <- gibbs_results1[[2]]

Delta2 <- gibbs_results2[[1]]
Theta2 <- gibbs_results2[[2]]

Delta3 <- gibbs_results3[[1]]
Theta3 <- gibbs_results3[[2]]

write.csv(Delta1, "delta1_gibbs")
write.csv(Theta1, "theta1_gibbs")
write.csv(Delta2, "delta2_gibbs")
write.csv(Theta2, "theta2_gibbs")
write.csv(Delta3, "delta3_gibbs")
write.csv(Theta3, "theta3_gibbs")


logQ1 <- variational_results1[[1]] 
M1 <- variational_results1[[2]]
logV1 <- variational_results1[[3]]

logQ2 <- variational_results2[[1]] 
M2 <- variational_results2[[2]]
logV2 <- variational_results2[[3]]

logQ3 <- variational_results3[[1]] 
M3 <- variational_results3[[2]]
logV3 <- variational_results3[[3]]


write.csv(logQ1, "logQ1_VA")
write.csv(M1, "M1_VA")
write.csv(logV1, "logV1_VA")
write.csv(logQ2, "logQ2_VA")
write.csv(M2, "M2_VA")
write.csv(logV2, "logV2_VA")
write.csv(logQ3, "logQ3_VA")
write.csv(M3, "M3_VA")
write.csv(logV3, "logV3_VA")


#-----------------------------------------------
#
#             Visualizations Gibbs
#
#-----------------------------------------------

d1_hat <- apply(Delta1, 2, mean)
t1_hat <- apply(Theta1, 2, mean)

df <- data.frame(cbind(d1_hat, t1_hat, c(rep(1, 10), rep(0, 500 - 10))))
colnames(df) <- c("Delta", "Theta", "Active")
df$Active <- as.factor(df$Active)
p1 <- ggplot(df, aes(x = Delta, y = Theta, col = Active))+
        geom_point(alpha = .8)+
        theme_minimal()+
        labs(x = "Estimated Delta Value", 
             y = "Estimate Theta Value",
            title = "Gibbs Sampler: p = 500")
p1

theta_true <- c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15, rep(0, 500 - 10))
p2 <- ggplot(data.frame(that = t1_hat, theta =  theta_true), 
             aes(x = theta, y = that, col = "coral"))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "True Theta ", 
       y = "Estimated Theta",
       title = "Gibbs Sampler: p = 500")
p2


d2_hat <- apply(Delta2, 2, mean)
t2_hat <- apply(Theta2, 2, mean)
theta_true <- c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15, rep(0, 1000 - 10))

df <- data.frame(cbind(d2_hat, t2_hat, c(rep(1, 10), rep(0, 1000- 10))))
colnames(df) <- c("Delta", "Theta", "Active")
df$Active <- as.factor(df$Active)
p3 <- ggplot(df, aes(x = Delta, y = Theta, col = Active))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "Estimated Delta Value", 
       y = "Estimate Theta Value",
       title = "Gibbs Sampler: p = 1000")
p3


p4 <- ggplot(data.frame(that = t2_hat, theta =  theta_true), 
             aes(x = theta, y = that, col = "coral"))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "True Theta ", 
       y = "Estimated Theta",
       title = "Gibbs Sampler: p = 1000")
p4


d3_hat <- apply(Delta3, 2, mean)
t3_hat <- apply(Theta3, 2, mean)
theta_true <- c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15, rep(0, 2000 - 10))

df <- data.frame(cbind(d3_hat, t3_hat, c(rep(1, 10), rep(0, 2000 - 10))))
colnames(df) <- c("Delta", "Theta", "Active")
df$Active <- as.factor(df$Active)
p5 <- ggplot(df, aes(x = Delta, y = Theta, col = Active))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "Estimated Delta Value", 
       y = "Estimate Theta Value",
       title = "Gibbs Sampler: p = 2000")
p5

p6 <- ggplot(data.frame(that = t3_hat, theta =  theta_true), 
             aes(x = theta, y = that, col = "coral"))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "True Theta ", 
       y = "Estimated Theta",
       title = "Gibbs Sampler: p = 2000")
p6

ggsave("gibbs1_estimates.pdf", plot = p1, width = 5, height = 5)
ggsave("gibbs1_theta.pdf", plot = p2, width = 5, height = 5)
ggsave("gibbs2_estimates.pdf", plot = p3, width = 5, height = 5)
ggsave("gibbs2_theta.pdf", plot = p4, width = 5, height = 5)
ggsave("gibbs3_estimates.pdf", plot = p5, width = 5, height = 5)
ggsave("gibbs3_theta.pdf", plot = p6, width = 5, height = 5)


#-----------------------------------------------
#
#             Visualizations CAVI
#
#-----------------------------------------------

q1hat <- exp(logQ1[nrow(logQ1),])
m1hat <- M1[nrow(M1),]
v1hat <- exp(logV1[nrow(logV1),])
theta_true <- c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15, rep(0, 500 - 10))

df <- data.frame(q = q1hat, mu = m1hat, 
                 v = v1hat, Active = as.factor(c(rep(1, 10), rep(0, 500 - 10))),
                 theta = theta_true)

p7 <- ggplot(df, aes(x = theta, y = mu, col = Active))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "True Coefficient Value",
       y = "Estimated Coefficient Value", 
       title = "CAVI: p = 500")

p7


q2hat <- exp(logQ2[nrow(logQ2),])
m2hat <- M2[nrow(M2),]
v2hat <- exp(logV2[nrow(logV2),])
theta_true <- c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15, rep(0, 1000 - 10))

df <- data.frame(q = q2hat, mu = m2hat, 
                 v = v2hat, Active = as.factor(c(rep(1, 10), rep(0, 1000 - 10))),
                 theta = theta_true)

p8 <- ggplot(df, aes(x = theta, y = mu, col = Active))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "True Coefficient Value",
       y = "Estimated Coefficient Value", 
       title = "CAVI: p = 1000")
p8


q3hat <- exp(logQ3[nrow(logQ3),])
m3hat <- M3[nrow(M3),]
v3hat <- exp(logV3[nrow(logV3),])
theta_true <- c(-15, -12, -10, -7, -3, 4, 6, 9 ,11, 15, rep(0, 2000 - 10))


df <- data.frame(q = q3hat, mu = m3hat, 
                 v = v3hat, Active = as.factor(c(rep(1, 10), rep(0, 2000 - 10))),
                 theta = theta_true)

p9 <- ggplot(df, aes(x = theta, y = mu, col = Active))+
  geom_point(alpha = .8)+
  theme_minimal()+
  labs(x = "True Coefficient Value",
       y = "Estimated Coefficient Value", 
       title = "CAVI: p = 2000")

p9


ggsave("CAVI1.pdf", plot = p7, width = 5, height = 5)
ggsave("CAVI2.pdf", plot = p8, width = 5, height = 5)
ggsave("CAVI3.pdf", plot = p9, width = 5, height = 5)


p10 <- ggplot(data.frame(v = v1hat, Active = as.factor(c(rep(1, 10), rep(0, 500 - 10)))), 
                  aes(x = Active, y = v)) + 
          geom_boxplot()+
          theme_minimal()+
          labs(x = "Active Status", 
                y = "Variance Estimate",
              title = "CAVI: p = 500")

p11 <- ggplot(data.frame(v = v2hat, Active = as.factor(c(rep(1, 10), rep(0, 1000 - 10)))), 
              aes(x = Active, y = v)) + 
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Active Status", 
       y = "Variance Estimate",
       title = "CAVI: p = 1000")

p12 <- ggplot(data.frame(v = v3hat, Active = as.factor(c(rep(1, 10), rep(0, 2000 - 10)))), 
              aes(x = Active, y = v)) + 
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Active Status", 
       y = "Variance Estimate",
       title = "CAVI: p = 2000")

ggsave("var1.pdf", plot = p10, width = 5, height = 5)
ggsave("var2.pdf", plot = p11, width = 5, height = 5)
ggsave("var3.pdf", plot = p12, width = 5, height = 5)



