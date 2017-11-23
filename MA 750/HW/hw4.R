#------------------------
#	1
#------------------------

#define a_m funcitons 
a0 = function(p){
	return(.75*(p-p^3/3 + 2/3))
}

a1 = function(p){
	return(.75*(p^2/2-p^4/4 - 1/4))
}

a2 = function(p){
	return(.75*(p^3/3-p^5/5 + 2/15))
}

#define Epanechnikov Kernel
K = function(t){
	return(.75*(1-t^2))
}

#Define indicator function
I = function(t,p){
	return(ifelse(t>p, 0,1))
}

#define Boundary kernel 
B = function(t,p){
	top = a2(p) - a1(p)*t
	bottom = a2(p)*a0(p) - a1(p)^2
	ret = top/bottom * K(t) * I(t,p)
	return(ret)
}

#plot boundary kernel
x = seq(-1, 1, by = .01)
y = list()
P = c(.25, .5, .75)
for(i in 1:3){
	p = P[i]
	y[[i]] = B(x, p)
}


plot(x, y[[1]], type = "l", lty = 1, ylab = "B(x)")
points(x, y[[2]], type = "l", lty = 2)
points(x, y[[3]], type = "l",  lty = 3)

legend("topleft", legend = c("p = .25", "p=.5", "p=.75"), lty = c(1,2,3))

#------------------------
#	2
#------------------------

#plot fourth order Epanechnikov kernel 
curve(175/32*(9/35 - 3*x^2/5)*(1-x^2), xlim = c(-1.1, 1.1), xlab = "", ylab = "", main = "Fourth Order Epanechnikov Kernel", col = "blue")


#------------------------
#	4
#------------------------

#get data 
pl = iris$Petal.Length

#define indicator and Epanechnikov kernel functions
I = function(x){
	return(ifelse(x> 1 || x< -1, 0, 1))
}

K = function(t){
	if(I(t) == 0){
		return(0)
	}else{
		return(0.75 *(1-t^2))
	}
}


#define Kappa
Kappa = function(upper){
	
	if(upper < -1){
		return(0)
	}
	
	if(upper > 1){
		return(1)
	}
	
	#Rough integral approximation 
	x = seq(-1, upper, by = 0.01)
	dx = 0.01
	y = K(x)
	return(sum(y*dx))
	
}

#define CDF 
CDF = function(x,h){
	arg = (x - pl)/h
	
	y = numeric(length(arg))
	for(i in 1:length(arg)){
		y[i] = Kappa(arg[i])
	}
	
	return(mean(y))
}

#plot results 
x = seq(1, 7, by = .1)
Y = list()
H = c(.25,.5,.75)
for(j in 1:length(H)){
	F = numeric(length(x))
	for(i in 1:length(x)){
		F[i] = CDF(x[i], H[j])
	}
	Y[[j]] = F
}

plot(x, Y[[1]], type = "l", main = "Epanechnikov Kernel CDF", lty = 2, col = "red", ylab = "F_n	(x)")
points(x, Y[[2]], type = "l", main = "Epanechnikov Kernel CDF", lty = 5,col = "green")
points(x, Y[[3]], type = "l", main = "Epanechnikov Kernel CDF", lty = 10, col = "blue")

legend("topleft", legend = c(".25", ".5", ".75"), lty = c(2,5,10))


#Build constant CDF estimator
# Compare this to the smooth kernel estimators 
histCDF = function(x){
	return(sum(x>pl)/length(pl))
}

histY = numeric(length(x))
for(i in 1:length(x)){
	histY[i] = histCDF(x[i])
}
plot(x, histY, type = "l", ylab = "F_hist(x)",main = "CDF Constant Estimator")









