#-------------------------------------------------------
#
#		Implementation 
#
#-------------------------------------------------------

library("ftnonpar")

#check out distribition
plot(dclaw(seq(-3, 3, len = 1000)), type = "l")

#get sample
X = rclaw(100)

#define kernel guassian kernel
K = function(y){
	return(1/sqrt(2*pi)*exp(-(y^2)/2))	
}
#define fhat 
fhat = function(x, h, est){
	n = length(X)
	return((1/(n*h))*sum(sqrt(est$y) * K((x-X)*sqrt(est$y)/h)))	
}

#define Kernel density estimate
f_var_bw = function(X, k, bw, h){
	
	#get initial bandwidth
	if(bw == "rot"){
		g = 1.06*sd(X)/length(X)^(1/5)
	}else{
		g = bw
	}
	
	#get initial estimate 
	ftil = density(X, kernel = k, bw = g)
	
	#approximate ftil(X) for our sample
	ftilX = approx(ftil$x, ftil$y, X)
	
	#get estimated kernel
	x = seq(-3,3,by = .05)
	f = matrix(NA, ncol = 2, nrow = length(x))
	for(i in 1:nrow(f)){
		f[i,1] = x[i]
		f[i,2] = fhat(x[i], h, ftilX)
	}

	#return estimate 
	return(f)
}

#-------------------------------------------------------
#
#		Find opt h 
#
#-------------------------------------------------------

H = seq(.05, .3, .05)
par(mfrow=c(3,2), mai = c(1, 0.1, 0.1, 0.1))
for(h in H){
	f = f_var_bw(X, "gaussian", "rot", h)
	plot(f[,1],f[,2], type = "l", main = paste("H: ", h))
}

H = seq(.1, .2, .02)
par(mfrow=c(3,2), mai = c(1, 0.1, 0.1, 0.1))
for(h in H){
	f = f_var_bw(X, "gaussian", "rot", h)
	plot(f[,1],f[,2], type = "l", main = paste("H: ", h))
}

hopt = .12 

#-------------------------------------------------------
#
#		Compare to fixed bandwidth 
#
#-------------------------------------------------------
 #find variable bandwidth (vbw) kernel
f_vbw = f_var_bw(X, "gaussian", "rot", hopt)

#find fixed bandwidth (fbw) kernel
rot_bw = 1.06 * sd(X) * (length(X))^(-1/5)
f_fbw = density(X, bw = rot_bw, kernel = "gaussian")

#Comparison to ftil
plot(f_vbw[,1],f_vbw[,2], col = "blue", type = "l", ylab = "", main = "Variable Bandwidth vs Fixed Bandwidth", xlab = "")
points(f_fbw, col = "orange", type = "l")
legend("topright",  c("Variable h", "Fixed h"), col = c("blue", "orange"), lty = 1)

#-----------------------------------------------------------
#
#		Compare starting kernels 
#
#-----------------------------------------------------------
Kern = c("gaussian","rectangular", "triangular")
BW = c(.1,.25,.5)

par(mfrow = c(1,3),  mai = c(1, 0.1, 0.1, 0.1))
for(k in 1:length(Kern)){
	plot(c(1,1), col = "white", xlim = c(-3,3), ylim = c(0,0.4), xlab = Kern[k], main = "")
	legend("topright",c("h = .1", "h = 0.25", "h = 0.5"), col = c("red", "green", "blue"), lty = 1)
	for(bw in 1:length(BW)){
		#get estimate 
		f = f_var_bw(X, Kern[k], BW[bw], hopt)
	
		#plot with different color
		points(f[,1], f[,2], col = bw+1, type = "l")
	}
}

par(mfrow = c(1,3),  mai = c(1, 0.1, 0.1, 0.1))
for(bw in 1:length(BW)){
	plot(c(1,1), col = "white", xlim = c(-3,3), ylim = c(0,0.4), xlab = paste("h = ",BW[bw]), main = "")
	legend("topright",c("Gaussian", "Uniform", "Triangular"), col = c("red", "green", "blue"), lty = 1)
	for(k in 1:length(Kern)){
		#get estimate 
		f = f_var_bw(X, Kern[k], BW[bw], hopt)
	
		#plot with different color
		points(f[,1], f[,2], col = k+1, type = "l")
	}
}

#--------------------------------------------
#
#	Opt H calculations
#
#--------------------------------------------

#package to estimate ~f''
library(kedd)

#define n 
n = length(X)

#Assuming Guassian Kernel 
Knorm = 1/(2*sqrt(pi)) 
mu2K = 1

#numerical integration for ||~f''||_2^2	
ftil2prime = dkde(X, kernel = "gaussian", deriv.order = 2, h = 1.06 * sd(X) * (length(X))^(-1/5))  

#appoximate area under curve
dx = diff(ftil2prime$eval.points)[1]
A = sum(dx * ftil2prime$est.fx^2)

#define ci
ftil = density(X, kernel = "gaussian", bw = 1.06*sd(X)/length(X)^(1/5))
ci = approx(ftil$x, ftil$y, X)$y

#fint opt h 
top = Knorm * mean(ci)
bottom = mu2K * A* (mean(ci^(-2)))^2
hopt = (top/bottom)^(1/5) * n^(-1/5)

#plot estimate 
est = f_var_bw(X, "gaussian", "rot", hopt)
plot(est[,1], est[,2], type = "l", ylab = "", xlab = "", main = "h = 0.077", col = "blue")




