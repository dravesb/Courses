#-----------------------------------------------------------------------
#
#                   Mouse Connectome
#
#-----------------------------------------------------------------------


setwd("~/Desktop/GLM_Final_Project/Mouse_Data")

#useful libraries-------------------------------------------------
library(igraph)
library(ggplot2)
library(reshape)
library(MASS)

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

#brain
brain = as.matrix(read.csv("N54859.csv", header = F))
n = nrow(brain)
colnames(brain) = rownames(brain) = 1:n
head(brain)

brain2 = matrix(0, nrow = nrow(brain), ncol = ncol(brain))
for(i in 1:nrow(brain2)){
  ind = as.numeric(which(brain[i,] == max(brain[i,])))
  brain2[i,ind] = brain[i, ind]
}

#covariates
covar = read.csv("first_page.csv")[,-c(1:2,9:11)]
covar = covar[-nrow(covar),]
for(i in 1:ncol(covar)){
  covar[,i] = as.character(droplevels(covar[,i]))
}
covar = as.matrix(covar)

#Make pretty visuals-------------------------------------------------
net = graph_from_adjacency_matrix(brain2,weighted = TRUE)
V(net)$color = c(rep("red", nrow(brain)/2), rep("blue", nrow(brain)/2))
V(net)$label = ""
V(net)$size = 2
E(net)$width = min(E(net)$weight, 100)/50
E(net)$color = "#55555555"
E(net)$arrow.mode = 0 
dev.off()
plot(net, layout = layout_on_sphere)


brain2 = log(brain + 1)
brainmelt = melt(brain2)
colnames(brainmelt) = c("Index","Index.","W")
ggplot(data = brainmelt, aes(x=Index, y=Index., fill=W)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = median(brain2), limit = c(0,max(brain2)))

#Build homogenous models-------------------------------------------------
#build response 
Y = brain[upper.tri(brain)]

#build covariates 
hemi = format_global(covar[,1])
level1 = format_global(covar[,2])
level2 = format_global(covar[,3])
level3 = format_global(covar[,4])
level4 = format_global(covar[,5])
subdivision = format_global(covar[,6])

hemi_group = ifelse(hemi[,1] == hemi[,2], 1, 0)
level1_group = ifelse(level1[,1] == level1[,2], 1, 0)
level2_group = ifelse(level2[,1] == level2[,2], 1, 0)
level3_group = ifelse(level3[,1] == level3[,2], 1, 0)
level4_group = ifelse(level4[,1] == level4[,2], 1, 0)
subdivision_group = ifelse(subdivision[,1] == subdivision[,2], 1, 0)

X = cbind(hemi, level1 , level2, level3, level4, subdivision) 
nested = numeric(nrow(X))
for(i in 1:nrow(X)){
  nested[i] = ifelse(X[i,1] != X[i,2], 0, 
              ifelse(X[i,3] != X[i,4], 1, 
              ifelse(X[i,5] != X[i,6], 2,
              ifelse(X[i,7] != X[i,8], 3,
              ifelse(X[i,9] != X[i,10], 4,
              ifelse(X[i,11] != X[i,12], 5,6))))))
}

df = data.frame(Response = Y,
                hemisphere = hemi_group, 
                level1 = level1_group,
                level2 = level2_group,
                level3 = level3_group,
                level4 = level4_group,
                subdivision = subdivision_group,
                nested = as.factor(nested))

#construct nested varible
m1 = glm(Response~hemisphere + level1 + level2 + level3, data = df, family = poisson)
m2 = glm(Response~hemisphere + level1 + level2 + level3 + level4, data = df, family = poisson)
m3 = glm(Response~hemisphere + level1 + level2 + level3 +level4 + subdivision, data = df, family = poisson)
#
anova(m1, m2, m3, m4, test = "Chisq")
sigma2 = sum(residuals(m3, "pearson")^2)/m3$df.residual
summary(m3, sigma2)

#fix overdispersion by using negative binomial with log link
#test difference between full Poisson and full NegBin
nb_full = glm.nb(Response ~ hemisphere +level1+ level2 + level3 +level4 +subdivision, data = df, link = log)
pchisq(2 * (logLik(nb_full) - logLik(m3)), df = 1, lower.tail = FALSE)
#p value is indistingishable from zero suggesting that the neg bin is a signficantly better fit
#we only consider nb models from here on

nb1 = glm.nb(Response ~ hemisphere, data = df, link = log)
nb2 = glm.nb(Response ~ hemisphere + level3 +level4, data = df, link = log)
nb3 = glm.nb(Response ~ hemisphere +level1+ level2 + level3 +level4, data = df, link = log)
nb4 = glm.nb(Response ~ hemisphere +level1+ level2 + level3 +level4 +subdivision, data = df, link = log)
anova(nb1, nb2, nb3,nb4, test = "Chisq")

x = c(AIC(nb1),AIC(nb1), AIC(nb3), AIC(nb4))
summary(nb4)
#subdivision appears to be too fine getting rid of it
summary(nb3) # best model by AIC and ANOVA standards

est = cbind(Estimates = coef(nb3), confint(nb3))
est
#going to consider adding nonhomogenous for groups 1 and 2 
#as their effects seem small

#Build inhomogenous models-------------------------------------------------
#build response 
Y = brain[upper.tri(brain)]

#make inhomogenous group labels
make_group_labels = function(x){
  X = format_global(x)
  names = unique(x)
  
  ret = numeric(nrow(X))
  for(i in 1:nrow(X)){
    if(X[i,1] == X[i,2]){
      ret[i] = names[which(names == X[i,1])]
    }else{
      ret[i] = 0
    }
  }
  return(ret)
}


hemi = make_group_labels(covar[,1])
level1 = make_group_labels(covar[,2])
level2 = make_group_labels(covar[,3])
level3 = make_group_labels(covar[,4])
level4 = make_group_labels(covar[,5])
subdivision = make_group_labels(covar[,6])


#g 
df = data.frame(Response = Y,
                hemisphere = as.factor(hemi), 
                level1 = as.factor(level1),
                level2 = level2_group,
                level3 = level3_group,
                level4 = level4_group 
                )

nb1 = glm.nb(Response ~(. - degree), data = df, link = log)
summary(nb1)
nb2 = glm.nb(Response ~., data = df, link = log)
summary(nb2)

anova(nb1, nb2, test = "Chisq")

#get model estimates 
Phat = matrix(0, nrow = n, ncol = n)
Phat[upper.tri(Phat)] = m2$fitted.values 
Phat[lower.tri(Phat)] = t(Phat)[lower.tri(Phat)]

Pmelt = melt(Phat)
colnames(Pmelt) = c("Index","Index.","Count")
ggplot(data = Pmelt, aes(x=Index, y=Index., fill=Count)) + 
  geom_tile()+
  scale_fill_gradient2(low = "white", 
                       high = "black", mid = "grey", 
                       midpoint = median(Phat), limit = c(0,max(Phat)))



#final model-------------------------------------------------------
df = data.frame(Response = Y,
                hemisphere = as.factor(hemi), 
                level1 = as.factor(level1),
                level2 = level2_group,
                level3 = level3_group,
                level4 = level4_group,
                level5 = subdivision_group
)

model = glm.nb(Response ~., data = df, link = log)
model$theta + qnorm(c(0.025, 0.975)) * model$SE.theta

fit = model$fitted.values
varest = fit + fit^2/model$theta
