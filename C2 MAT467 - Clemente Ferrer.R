# Clemente Ferrer
# 201910002-5


## Librerias
library(jpeg)
library(dichromat)
library(GeoModels)


## Lectura de las fotografías
setwd("C:/Users/ccfer/Documents/MAT467/Recortadas")

files <- list.files(path = "C:/Users/ccfer/Documents/RONY/Recortadas", pattern=".jpg", all.files = T)
harvard <- list()

for (i in 1:length(files))
{
  harvard[[i]] <- readJPEG(files[i])
}


## Función GCC
gcc <- function(y) {
  gcc_y <- array(numeric(), dim = c(nrow(y),ncol(y)))
  for (i in 1:nrow(y)){
    for (j in 1:ncol(y)){
    gcc_y[i,j] = (y[i,j,2])/(y[i,j,1]+y[i,j,2]+y[i,j,3])
    }
  }
  return(gcc_y)
}


## GCC fotografía recortadas
image(gcc(harvard[[1]]), axes=FALSE, col= colorRampPalette(c("white", "Green"))(64));
image(gcc(harvard[[8]]), axes=FALSE, col= colorRampPalette(c("white", "Green"))(64));
image(gcc(harvard[[11]]), axes=FALSE, col= colorRampPalette(c("white", "Green"))(64));
image(gcc(harvard[[15]]), axes=FALSE, col= colorRampPalette(c("white", "Green"))(64));


## Ajuste modelo autoregresivo espacial (código basado en los estimadores 
# hallados en la tarea 6, pregunta 1)

a=0
b=0
c=0
d=0
e=0
f=0

gcc_2018 = gcc(harvard[[11]])

for (i in 2:40){
  for (j in 2:40){
    a=a+(gcc_2018[i-1,j])^2
    b=b+(gcc_2018[i,j-1]*gcc_2018[i-1,j])
    c=c+(gcc_2018[i,j]*gcc_2018[i-1,j])
    d=d+(gcc_2018[i-1,j]*gcc_2018[i,j-1])
    e=e+(gcc_2018[i,j-1])^2
    f=f+(gcc_2018[i,j]*gcc_2018[i,j-1])
  }
}

phi1 = (b*f-c*e)/(b*d-a*e)
phi2 = (c*d-a*f)/(b*d-a*e)
print(c(phi1,phi2)) #Parámetros estimados

res=matrix(0,nrow = nrow(gcc_2018),ncol=ncol(gcc_2018))
for (i in 2:nrow(gcc_2018)){
  for (j in 2:ncol(gcc_2018)){
    res[i,j]=gcc_2018[i,j]-phi1*gcc_2018[i-1,j]-phi2*gcc_2018[i,j-1]
  }
}

image(res[2:40,2:40], col= colorRampPalette(c("white", "Green"))(64))


## Grilla de 40x40 ajustada a imagen GCC
grilla40x40 <- function(a, n){
  b <- a[1:n, 1:n]
  return(b)
}

harvard_adjusted_gcc <- list()

for (i in 1:length(harvard)){
  aux <- gcc(harvard[[i]])
  harvard_adjusted_gcc[[i]] <- grilla40x40(aux,40)
  image(harvard_adjusted_gcc[[i]], axes=FALSE, col= colorRampPalette(c("white", "Green"))(64))
}


## Modelo espacio-temporal con covarianza de Gneiting 
CorrParam("Gneiting")

model="Gaussian"
corrmodel="Gneiting"

coordt = 1:15 #Número de instantes temporales 
T = length(coordt)
NN = 1600 #Número de localizaciones espaciales
x = seq(1/40,1,1/40)
y = seq(1/40,1,1/40)
coords=merge(x,y)

lista = c()
for (i in 1:T){
  vec <- rep(i, NN)
  lista <- append(lista, vec) 
}

X=cbind(rep(1,NN*T),lista)

data = matrix(0,nrow=15,ncol=(40*40))
for (j in 1:15){
  data2 = c()
  for (i in 1:40){
    data2 = c(data2, harvard_adjusted_gcc[[j]][1:40,i])
  }
  data[j,] = data2
}

start=list(mean = 0.5, mean1 = 0.5, scale_s=0.2, scale_t = 2, sill=2, power_s = 0.2, power_t = 0.2, sep= 0.5)
fixed=list(nugget=0)

fit = GeoFit(data=data,coordx=coords, coordt=coordt, corrmodel=corrmodel,
             maxdist=0.05,maxtime=1,X=X,
             optimizer="BFGS",start=start,fixed=fixed)


### Análisis de residuos
res=GeoResiduals(fit)
mean(res$data) #approx 0

vario = GeoVariogram(data=res$data,coordx=coords, coordt=coordt, maxdist=20,maxtime=1)
GeoCovariogram(res,vario=vario,fix.lagt=1,fix.lags=1,show.vario=TRUE,pch=20)

GeoQQ(res)


## Predicción t=16
coords = as.matrix(coords)
n_loc = nrow(coords)
times = c(16)
n_tim = length(times)
Xloc = cbind(rep(1,n_loc*n_tim),rep(16, n_loc*n_tim))

param_est = as.list(c(fit$param,fixed))
data3 = data[14:15,]
coordt2 = c(14:15)
X2 = X[((40*40*13)+1):(40*40*15),]

pr = GeoKrig(data=data3,coordx=coords, coordt=coordt2, corrmodel=corrmodel,
             X=X2,Xloc=Xloc,sparse=FALSE,
             model=model,mse=TRUE,loc=coords,time=times,param=param_est)

par(mfrow=c(1,2))
colour <- colorRampPalette(c("white", "Green"))(64)
image.plot(x, y, matrix(pr$pred,ncol=length(x)),col=colour, main = paste("SP Kriging at Time=" , 16),ylab="")
image.plot(x, y, matrix(pr$mse,ncol=length(x)),col=colour,
           main = paste("MSE at Time=" , 16),ylab="")

#El último código está basado en los ejemplos de GeoModels R Package.

