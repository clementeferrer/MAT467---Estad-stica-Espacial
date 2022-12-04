#Librerias a utilizar

library(gstat)
library(sp)
library(latex2exp)
library(MASS)
library(spatstat)

#P1

#Matriz de combinaciones de valores de phi. 
phis <-function(n){
  phis_comb=c()
  h=1/n
  for (i in 1:(n-1)){
    phis_comb = c(phis_comb,h*i)}
    phis_comb = merge(phis_comb,phis_comb)
    return(phis_comb) 
}

phis = phis(5)

#Imagen para cada combinación (phi_1, phi_2)

for(k in (1:nrow(phis))){
  n = 101
  X = matrix(rnorm((n-1)^2,0,1), ncol=(n-1), nrow=(n-1))
  phi1=phis[k,1]
  phi2=phis[k,2]
    for (i in 2:(n-1)){
      for(j in 2:(n-1)){
        X[i,j]=phi1*X[i-1,j]+phi2*X[i,j-1]+X[i,j]
  }
}

image(X, col = hcl.colors(11, "purples", rev = TRUE), main = paste("phi1 = ", format(phi1), "phi2 =", format(phi2)))
}

#P3

#Plot dataset
plot(anemones, markscale=1, main = "Patrón puntual planar de Anemonas")

#Quadrat Count Test for CSR
QCI = as.im(anemones,dimyx=6)
plot(QCI,main="Anemones")
points(anemones, bg=0, pch=21, col="black")
qtest=quadrat.test(anemones,6,6)
qtest

#Clark Evans
clarkevans(anemones)
clarkevans.test(anemones)

#K function
K <- envelope(anemones, Kest, nsim=100, fix.n=TRUE)
plot(K,main="K function")
G <- envelope(anemones, Gest, nsim=100, fix.n=TRUE)
plot(G,main="G function")

#P4

#Poisson processes
lambda <- function(x,y) { exp(5*x+2*y) }
lambda_const <- 100
pp <- rpoispp(lambda, win=square(1))
pp_const <- rpoispp(lambda_const, win=square(1))

plot(density(pp), main = "Proceso de Poisson no homogéneo")
plot(pp, bg=0, pch=21, add = TRUE)

plot(density(pp_const), main = "Proceso de Poisson homogéneo")
plot(pp_const, bg=0, pch=21, add = TRUE)

#K function
K <- envelope(pp, Kest, nsim=40, fix.n=TRUE)
plot(K)
