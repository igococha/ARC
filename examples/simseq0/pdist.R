library(ggplot2)



x <- seq(-4, 4, length=100)
hx <- dnorm(x)



alpha <- 0.5
theta <- 2.0
beta <- 1/theta
x <-  seq(0, 5, length=100)


dx <- dgamma(x,alpha,beta)

plot(x, dx, type="l", lty=2, xlab="x value",
  ylab="CDF", main="Gamma Gistribution: DCF")


px <- pgamma(x,alpha,beta)


plot(x, px, type="l", lty=2, xlab="x value",
  ylab="CDF", main="Gamma Gistribution: DCF")


# this is giving problems 

kappa <- 0.0026641698338647143 
theta <- 1.370734568776759

alpha <- kappa
beta <- 1/theta

x <- seq(0, 0.005, length=1500)
dx <- dgamma(x,alpha,beta)

plot(x, dx, type="l", lty=2, xlab="x value",
  ylab="density", main="Gamma Gistribution: dgamma")


px <- pgamma(x,alpha,beta)


plot(x, px, type="l", lty=2, xlab="x value",
  ylab="distribution", main="Gamma Gistribution: pgamma")

p <- 0.125

qgamma(p,alpha,beta)

