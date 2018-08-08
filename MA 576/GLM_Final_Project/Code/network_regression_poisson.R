library(igraph)
library(ggplot2)
library(reshape)
library(MASS)
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

#set up interesting sample graph-----------------------------------------
n = 25
size.gp1 = 10
size.gp2 = n - size.gp1

#make covariates based on groups
x1 = rpois(n, 2)
x2 = c(rep(1, size.gp1),rep(2, size.gp2))

X1 = format_global(x1)
X2 = format_global(x2)

X = cbind(X1, X2)
colnames(X) = rep("", ncol(X))
rownames(X) = rep("", nrow(X))
#create probability matrix
combine = function(x){
  ret = ifelse(x[3] == x[4], 10, 0) + (x[1]+x[2])
  return(ret)
}
eta = apply(X, 1,combine)
hist(eta)

L = matrix(0, nrow = n, ncol = n)
L[upper.tri(L)] = eta
L[lower.tri(L)] = t(L)[lower.tri(L)]

#make pretty Lambda plot
Lmelt = melt(L)
colnames(Lmelt) = c("Index","Index.","L")
ggplot(data = Lmelt, aes(x=Index, y=Index., fill=L)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = median(L), limit = c(0,max(L)))

#set up sample function
pois_sample = function(l){
  return(rpois(1,l))	
}
A_sample = function(L){
  A = apply(L,c(1,2),pois_sample)
  A[lower.tri(A)] = t(A)[lower.tri(A)]
  diag(A) = 0
  return(A)
}

#get sample and make some nice visuals
A = A_sample(L)

Amelt = melt(A)
colnames(Amelt) = c("Index","Index.","W")
ggplot(data = Amelt, aes(x=Index, y=Index., fill=W)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = median(A), limit = c(0,max(A)))



net = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
V(net)$color = c(rep("red", size.gp1),rep("blue", size.gp2))
V(net)$label = ""
V(net)$size = 5
E(net)$width = E(net)$weight/median(E(net)$weight)
E(net)$color = "#55555555"
E(net)$arrow.mode = 0 
dev.off()
plot(net)

#set up poisson regression model-----------------------------
#set up response
A = as.matrix(as_adjacency_matrix(net, attr = 'weight'))
Y = A[upper.tri(A, diag = F)]

#set up inter-nodal attributes
#group
group = ifelse(X[,3] == X[,4], 1, 0)

#Construct model---------------------------------------------------
#group only model
df = data.frame(Response = Y, SameGroup = group)
m4 = glm(Response~.,data = df, family = poisson)
summary(m4)

#full model
df = data.frame(Response = Y, SameGroup = group, Cont = X[,1] + X[,2])
m5 = glm(Response~.,data = df, family = poisson)
summary(m5)

#construct test for overdispersion
sigma2 = sum(residuals(m5, type = "pearson")^2/m5$df.residual)
teststat = sigma2 * m5$df.residual
cuttoff = qchisq(.95, m5$df.residual)
reject = teststat > cuttoff
reject

#quassi poisson model
m6 = glm(Response~.,data = df, family = quasipoisson)
summary(m6)

#negative binomial model
m7 = glm.nb(Response~.,data = df, link = log)
summary(m7)

#construct test for overdispersion
sigma2 = sum(residuals(m7, type = "pearson")^2/m7$df.residual)
teststat = sigma2 * m7$df.residual
cuttoff = qchisq(.95, m7$df.residual)
reject = teststat > cuttoff
reject

#compare models
pchisq(-2*(logLik(m6) - logLik(m7)), 1, lower.tail = FALSE)

#get model estimates 
Lhat = matrix(0, nrow = n, ncol = n)
Lhat[upper.tri(Lhat)] = m7$fitted.values
Lhat[lower.tri(Lhat)] = t(Lhat)[lower.tri(Lhat)]

Lmelt = melt(Lhat)
colnames(Lmelt) = c("Index","Index.","L")
ggplot(data = Lmelt, aes(x=Index, y=Index., fill=L)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = median(L), limit = c(0,max(L)))


