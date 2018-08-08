#-----------------------------------------------------------------------
#
#                   Bernoulli Simulation
#
#-----------------------------------------------------------------------

library(igraph)
library(ggplot2)
library(reshape)
#set seed 
set.seed(576)

#Define useful functions-------------------------------------------------
format_global = function(x){
	tmp = outer(x, x, "paste")
	tmp = tmp[upper.tri(tmp)]
	split = function(a) as.numeric(unlist(strsplit(a, " ")))
	X = t(sapply(tmp, split))
	return(X)		
}
getp = function(x) exp(x)/(1+exp(x))

#set up interesting sample graph-----------------------------------------
n = 50
size.gp1 = 20
size.gp2 = n - size.gp1

#make covariates based on groups
x1 = rpois(n, 4)
x2 = c(rep(1, size.gp1),rep(2, size.gp2))

X1 = format_global(x1)
X2 = format_global(x2)

X = cbind(X1, X2)
colnames(X) = rep("", ncol(X))
rownames(X) = rep("", nrow(X))
#create probability matrix
combine = function(x){
  ret = ifelse(x[3] + x[4] == 2, -2, ifelse(x[3] +x[4] == 4 ,-1,0)) +   
    +1/10*(x[1] + x[2]) 
  return(ret)
}
eta = apply(X, 1,combine)
hist(eta)

P = matrix(0, nrow = n, ncol = n)
P[upper.tri(P)] = pnorm(eta)
P[lower.tri(P)] = t(P)[lower.tri(P)]

#make pretty P plot
Pmelt = melt(P)
colnames(Pmelt) = c("Index","Index.","P")
ggplot(data = Pmelt, aes(x=Index, y=Index., fill=P)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1)) 

#set up sample function
bern_sample = function(p){
	return(rbinom(1,1,p))	
}
A_sample = function(P){
	A = apply(P,c(1,2),bern_sample)
	A[lower.tri(A)] = t(A)[lower.tri(A)]
	diag(A) = 0
	return(A)
}

#get sample and make some nice visuals
A = A_sample(P)

Amelt = melt(A)
colnames(Amelt) = c("Index","Index.","A")
ggplot(data = Amelt, aes(x=Index, y=Index., fill=A)) + 
  geom_tile()+ 
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1)) 

par(mfrow = c(1,1))
net = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
V(net)$color = c(rep("red", size.gp1),rep("blue", size.gp2))
V(net)$label = ""
V(net)$size = 5
E(net)$width = .3
E(net)$color = "#55555555"
E(net)$arrow.mode = 0 
plot(net)

#set up logistic regression model-----------------------------
#set up response
A = as.matrix(as_adjacency_matrix(net))
Y = A[upper.tri(A, diag = F)]

#set up inter-nodal attributes
#group
group = ifelse(X[,3] == X[,4] & X[,4] == 2, 2, ifelse(X[,3] == X[,4] & X[,4] == 1, 1, 0))

#Construct model---------------------------------------------------

#homogenous group effect 
df = data.frame(Response = Y, SameGroup = ifelse(X[,3] == X[,4], 1,0))
m1 = glm(Response~.,data = df, family = binomial)
summary(m1)
par(mfrow = c(2,2)); plot(m1)

#get prob predictions
eta_diff = coef(m1)[1]
eta_same = sum(coef(m1))

getp(eta_same)
getp(eta_diff)

Phat = matrix(0, nrow = n, ncol = n)
Phat[upper.tri(Phat)] = m1$fitted.values
Phat[lower.tri(Phat)] = t(Phat)[lower.tri(Phat)]

#make pretty P plot
Pmelt = melt(Phat)
colnames(Pmelt) = c("Index","Index.","P")
ggplot(data = Pmelt, aes(x=Index, y=Index., fill=P)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1)) 

#construct test for overdispersion
sigma2 = sum(residuals(m1, type = "pearson")^2/m1$df.residual)
teststat = sigma2 * m1$df.residual
cuttoff = qchisq(.95, m1$df.residual)
reject = teststat > cuttoff


#Construct model---------------------------------------------------

#inhomogenous group effect 
df = data.frame(Response = Y, SameGroup = as.factor(group))
m2 = glm(Response~.,data = df, family = binomial)
summary(m2)
par(mfrow = c(2,2)); plot(m1)

Phat = matrix(0, nrow = n, ncol = n)
Phat[upper.tri(Phat)] = m1$fitted.values
Phat[lower.tri(Phat)] = t(Phat)[lower.tri(Phat)]

#construct test for overdispersion
sigma2 = sum(residuals(m2, type = "pearson")^2/m2$df.residual)
teststat = sigma2 * m2$df.residual
cuttoff = qchisq(.95, m2$df.residual)
reject = teststat > cuttoff



#make CIs
C = vcov(m2)
p0 = getp(coef(m2)[1])
p0l = getp(coef(m2)[1] - qnorm(0.975)*sqrt(C[1,1]))
p0r = getp(coef(m2)[1] + qnorm(0.975)*sqrt(C[1,1]))

p1 = getp(coef(m2)[1] + coef(m2)[2])
p1l = getp(coef(m2)[1] +coef(m2)[2] - qnorm(0.975)*sqrt(C[1,1] + C[2,2] + 2*C[1,2]))
p1r = getp(coef(m2)[1] +coef(m2)[2] + qnorm(0.975)*sqrt(C[1,1] + C[2,2] + 2*C[1,2]))

p2 = getp(coef(m2)[1] + coef(m2)[3])
p2l = getp(coef(m2)[1] +coef(m2)[3] - qnorm(0.975)*sqrt(C[1,1] + C[3,3] + 2*C[1,3]))
p2r = getp(coef(m2)[1] +coef(m2)[3] + qnorm(0.975)*sqrt(C[1,1] + C[3,3] + 2*C[1,3]))

c(p0l,p0r)
c(p1l,p1r)
c(p2l,p2r)

#include poisson covartiate-----------------------------------------------
df = data.frame(Response = Y, 
                SameGroup = as.factor(group),
                GlobalPois = as.numeric(X[,1]) +as.numeric(X[,2]))
m3 = glm(Response~.,data = df, family = binomial)
summary(m3)
par(mfrow = c(2,2)); plot(m3)

Phat = matrix(0, nrow = n, ncol = n)
Phat[upper.tri(Phat)] = m3$fitted.values
Phat[lower.tri(Phat)] = t(Phat)[lower.tri(Phat)]

Pmelt = melt(Phat)
colnames(Pmelt) = c("Index","Index.","P")
ggplot(data = Pmelt, aes(x=Index, y=Index., fill=P)) + 
  geom_tile()+ 
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = 0.5, limit = c(0,1)) 

#test for which model is best
anova(m1, m2,m3, test = "Chisq")




