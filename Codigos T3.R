library(ape)
library(scatterplot3d)
library(gstat)
library(MASS) #Box Cox

#Lectura del archivo
santiago <- read.csv(file = "C:/Users/ccfer/Downloads/Santiago.csv", sep = ';', dec = ',')  
head(santiago)
colnames(santiago)

#Grafico
with(santiago, 
     scatterplot3d(ï..Latitud,
                   Longitud, 
                   KPI, 
                   main="KPI",
                   xlab = "Latitud",
                   ylab = "Longitud",
                   zlab = "KPI",
                   pch = 20, color="red",type="h"))

plot(santiago$ï..Latitud,santiago$Longitud, xlab="Latitude", ylab="Longitude") 
plot(santiago$ï..Latitud, santiago$KPI, xlab="Latitude", ylab="KPI") 
plot(santiago$Longitud, santiago$KPI, xlab="Longitude", ylab="KPI") 

#Histograma: no normal
hist(santiago$KPI)

#Box Cox
b <- boxcox(lm(santiago$KPI ~ 1))
lambda <- b$x[which.max(b$y)]
nuevo_KPI <- (santiago$KPI ^ lambda - 1) / lambda
hist(nuevo_KPI)
shapiro.test(nuevo_KPI)
summary(santiago$KPI)

#Indice de Moran

santiago.dists <- as.matrix(dist(cbind(santiago$Longitud, santiago$ï..Latitud)))
santiago.dists.inv <- 1/santiago.dists
diag(santiago.dists.inv) <- 0

Moran.I(santiago$KPI, santiago.dists.inv)

#CLOUD VARIOGRAM
coordinates(santiago) = ~ï..Latitud+Longitud

# Variograma usando gstat
g <- gstat(id="KPI", formula = KPI~1, data = santiago)
stgo.vgm = variogram(g)
plot(stgo.vgm)
stgo.fit = fit.variogram(stgo.vgm, model = vgm("Mat"),fit.kappa=TRUE)
plot(stgo.vgm,stgo.fit)

# Updating the variogram
variog<-gstat(g, id="KPI",model=stgo.fit)

# Constructing the grid for prediction
grilla = expand.grid( x=seq(min(santiago$ï..Latitud), max(santiago$ï..Latitud), l=200), y=seq(min(santiago$Longitud), max(santiago$Longitud), l=200) )
gridded(grilla) = ~x+y

#Predictions
pred <- predict(variog, grilla, interval="confidence", level=0.95)

#plotting predictions
spplot(pred["KPI.pred"])
spplot(pred["KPI.var"])

pred
