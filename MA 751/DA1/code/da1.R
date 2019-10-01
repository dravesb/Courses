#----------------------------------
#
#         Data Assignment 1 
#
#----------------------------------


#read in data 
train <- read.csv("~/Documents/Work/github/Courses/MA 751/DA1/Data/vowel.train.txt")
test <- read.csv("~/Documents/Work/github/Courses/MA 751/DA1/Data/vowel.test.txt")

#summarize data
head(train)

#clean data
train <- as.matrix(train[,-1])
test <- as.matrix(test[,-1])

#----------------------------------
#
#         QDA Implemenation
#
#----------------------------------

#function: est - return mean/covariance estimates
#arguments: 
#       X - data matrix to take mean/covariance of

est <- function(X){
  mu <-colMeans(X)
  Sigma <- cov(scale(X, scale = FALSE))
  return(list(mu = mu, Sigma = Sigma))
}

#function: deta - return discriminant function
#arguments: 
#       X - data matrix to take mean/covariance of

delta <- function(x, mu, Sigma, lp){
  
    #numerically stable matrix inversion
    C <- chol(Sigma)
    y <- backsolve(C, matrix(x - mu, ncol = 1), trans = TRUE)
    b <- backsolve(C,y)
    
    #return discriminant function 
    -.5 * log(det(Sigma)) - .5 * crossprod(x - mu, b) + lp 
    
}


#function: QDA - return classification 
#arguments: 
#       X - training design matrix 
#       Y - training responses
#       x - test point
QDA <-local({
  #calculate prior class probabilities
  lpi <- log(as.numeric(prop.table(table(train[,1]))))
  
  #calculate class means & covariances
  K <- length(unique(train[,1]))
  mus <- list(K)
  Sigmas <- list(K)
  
  for(k in 1:K){
    #subset dataframe & get estimates 
    df <- train[,-1][which(train[,1] == k), ]
    tmp <- est(df)
    
    #store estimates
    mus[[k]] <- tmp$mu
    Sigmas[[k]] <- tmp$Sigma
    
  }
  
  function(x) which.max(sapply(1:K, function(y) delta(x, mus[[y]], Sigmas[[y]], lpi[y])))
  
})

#----------------------------------
#
#        Applicaiton to 
#           Test Set 
#
#----------------------------------

#format test data
test.X <- as.matrix(test[,-1]); colnames(test.X) = rep("", ncol(test.X))

#get predictions
yhat <- apply(test.X, 1, QDA)

#get get misclassification rate
MR <- sum(yhat != test[,1])/nrow(test)
MR

#make figure
df1 <- data.frame(X1 = test.X[,1], X2 = test.X[,2], class = as.factor(yhat))
df2 <- data.frame(X1 = test.X[,1], X2 = test.X[,2], class = as.factor(test[,1]))

p1 <- ggplot(df1, aes(X1, X2, col = class))+
  geom_point()+
  labs(title = "Predicted Classifications")+
  theme_bw()

p2 <- ggplot(df2, aes(X1, X2, col = class))+
  labs(title = "Truth")+
  geom_point()+
  theme_bw()

grid.arrange(p1, p2, ncol=2)

#make grid values
dummydf <- 




#----------------------------------
#
#         LDA Implemenation
#
#----------------------------------

#Set up new training/testing data
train.new <- matrix(NA, nrow = nrow(train) , ncol  = 1 + 10 + 10 + choose(10, 2))
train.new[,1:11] <- as.matrix(train)

test.new <- matrix(NA, nrow = nrow(test) , ncol  = 10 + 10 + choose(10, 2))
test.new[,1:10] <- as.matrix(test[,-1])


#fill in all cross terms and quadratic terms
count <- 1
for(i in 1:10){
  for(j in i:10){
    train.new[,11 + count] <- train[,i + 1] * train[,j + 1]
    test.new[,10 + count] <- test[,i + 1] * test[,j + 1]
    count <- count + 1
  }
}

#function: deta - return discriminant function
#arguments: 
#       X - data matrix to take mean/covariance of

deltaLDA <- local({
  
    #calculate sigma
    Sigma <- cov(train.new[,-1])
  
    function(x, mu, lp){
      #numerically stable matrix inversion
      C <- chol(Sigma)
      y <- backsolve(C, matrix(mu, ncol = 1), trans = TRUE)
      b <- backsolve(C,y)
  
      #return discriminant function 
      crossprod(x, b) - .5 * crossprod(mu, b) + lp 
    }
})

#function: QDA - return classification 
#arguments: 
#       x - test point
LDA <-local({
  #calculate prior class probabilities
  lpi <- log(as.numeric(prop.table(table(train.new[,1]))))
  
  #calculate class means & covariances
  K <- length(unique(train.new[,1]))
  mus <- list(K)
  for(k in 1:K){
    #subset dataframe & get estimates 
    df <- train.new[,-1][which(train.new[,1] == k), ]
    
    #store estimates
    mus[[k]] <- colMeans(df)
    
  }
  
  #return classifications
  function(x) which.max(sapply(1:K, function(y) deltaLDA(x, mus[[y]], lpi[y])))
  
})



#----------------------------------
#
#        Applicaiton to 
#           Test Set 
#
#----------------------------------

#get predictions
yhat.new <- apply(test.new, 1, LDA)

#get get misclassification rate
MR.new <- sum(yhat.new != test[,1])/nrow(test)
MR.new

#make figure 
df1 <- data.frame(X1 = test.new[,1], X2 = test.new[,2], class = as.factor(yhat.new))
df2 <- data.frame(X1 = test.new[,1], X2 = test.new[,2], class = as.factor(test[,1]))

p3 <- ggplot(df1, aes(X1, X2, col = class))+
  geom_point()+
  labs(title = "Predicted")+
  theme_bw()

p4 <- ggplot(df2, aes(X1, X2, col = class))+
  labs(title = "Truth")+
  geom_point()+
  theme_bw()

grid.arrange(p3, p4, ncol=2)



library(MASS)
df <- data.frame( Y = test[,1], X = test.new[,1], Z = test.new)
model <- lda(Y~ ., data=df)
decisionplot(model, df, class = "X2")



