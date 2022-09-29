library(ape)
library(scatterplot3d)
library(gstat)

santiago <- read.csv(file = "C:/Users/ccfer/Downloads/Santiago.csv", sep = ';', dec = ',')  
head(santiago)
colnames(santiago)
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
hist(santiago$KPI)
summary(santiago$KPI)

santiago.dists <- as.matrix(dist(cbind(santiago$Longitud, santiago$ï..Latitud)))
santiago.dists.inv <- 1/santiago.dists
diag(santiago.dists.inv) <- 0

Moran.I(santiago$KPI, santiago.dists.inv)

#CLOUD VARIOGRAM
coordinates(santiago) = ~ï..Latitud+Longitud
vario.cloud<-variogram(KPI ~1,santiago,cloud=TRUE)
plot(vario.cloud)

#EMPIRICAL VARIOGRAM
vario<-variogram(KPI ~1,santiago)
plot(vario)

#Cressie Hawkins VARIOGRAM
varioRob<-variogram(KPI ~1,santiago, cressie=TRUE)
plot(varioRob)

#BOXPLOT
boxplot(santiago$KPI)

# Fitting a Matern variogram to the variable KPI
fit1=fit.variogram(vario,vgm("Mat"), fit.kappa = TRUE)
fit1

#theoretical variogram model
V1<-vgm(psill=10287.33,model="Mat", range=0.03560912, nugget = 12000)

# Plotting bothe variogram in the same graph
plot(vario,V1)