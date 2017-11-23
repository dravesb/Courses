#-------------------------
#
#		Exercise 9.4  
#
#-------------------------

#initialize variables
n = 200
lambda = 5

#get uniform distribution sample
x = runif(n)

#Define F^{-1}(U) for Exp(\lambda)
Finverse = function(x, lambda){
	y = -log(1-x)/lambda
}

#get F^{-1}(U) and a sample from the Exp in R 
y = Finverse(x, lambda)
z = rexp(n, 5)

#plot kernel densities
plot(density(z), lty = 2, main = "Density Estimates",xlab = "x", ylab = "f(x)")
points(density(y), lty =3,type = "l")
legend("topright", legend = c("Exponential", "F-1(Unif)"), lty = c(2,3))

