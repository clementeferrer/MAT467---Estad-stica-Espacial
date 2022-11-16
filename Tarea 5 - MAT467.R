# Librerias a utilizar
library(fastmatrix)
library(SpatialPack)
library(akima)
library(fields)
library(geoR)
library(GeoModels)

#Lectura y manipulación de los datos
Santiago <- read.csv2("Santiago.csv")
coords <- Santiago[c("Longitud","ï..Latitud")]
x <- Santiago$Potencia.Nominal
y <- Santiago$KPI

#Test de hipótesis correlación espacial. Clifford, Richardson (1989)
hip.test <- modified.ttest(x,y,coords)
hip.test

#Gráfico de codispersión
codisp.graph <-codisp(x,y,coords)
plot(codisp.graph)

#Coeficiente de Tjøstheim's 
santiago.cor<- cor.spatial(x, y, coords)
santiago.cor

#Mapa de codispersión
codisp.map <- function(x, y, coords, nclass = 2, ncell = 300, plot.it = TRUE)
  {
    require("akima")
    require("fields")
    require("geoR")
    rho <- function(x, y, uvec, max.dist, angle)
    {
      z <- as.geodata(cbind(x$coords, x$data + y$data))
      nz <- variog(z, uvec = uvec, max.dist = max.dist, direction = angle, messages = FALSE)
      dx <- variog(x, uvec = uvec, max.dist = max.dist, direction = angle, messages = FALSE)
      dy <- variog(y, uvec = uvec, max.dist = max.dist, direction = angle, messages = FALSE)
      rho <- .5 * (nz$v - dx$v - dy$v) / sqrt(dx$v * dy$v)
    }
    
    x <- as.geodata(cbind(coords, x))
    y <- as.geodata(cbind(coords, y))
    dmax <- .5 * max(dist(coords))
    angles <- seq(from = 0, to = pi, by = 0.01)
    nangles <- length(angles)
    uvec <- seq(from = 0, to = dmax, length = nclass + 1)[-1]
    
    xcirc <- 0
    ycirc <- 0
    
    for (i in seq_len(nclass)) {
      xcirc[(nangles*(i-1)+1):(nangles*i)] <- seq(-uvec[i], uvec[i], length = nangles)
      ycirc[(nangles*(i-1)+1):(nangles*i)] <- sqrt(uvec[i]^2 - xcirc[(nangles*(i-1)+1):(nangles*i)]^2)
    }
    z <- matrix(0, nrow = nangles, ncol = nclass)
    for (i in seq_len(nangles))
      z[i,] <- rho(x, y, uvec = uvec, max.dist = dmax, angle = angles[i])
    z <- as.vector(z)
    xl <- seq(min(xcirc), max(ycirc), length=ncell)
    yl <- seq(min(ycirc), max(ycirc), length=ncell)
    
    if (plot.it) {
      par(pty = "s")
      image.plot(interp(xcirc, ycirc, as.vector(z), xo = xl,yo = yl), col = topo.colors(256),
                 xlab = "x", ylab = "y")
      title(main = "Codispersion Map")
    }
    
    invisible(list(xcirc = xcirc, ycirc = ycirc, z = z))
  }

codis=codisp.map(x,y,coords)

#Matern bivariado (Extraído de Vallejos, Osorio y Bevilacqua (2020))

data = data.frame(x,y)

corrmodel = "Bi_matern"

CorrParam("Bi_matern")

start <- list(sill_1 = var(data[,1]), sill_2 = var(data[,2]), scale_1 =0.5, scale_2 = 0.5, scale_12 = 0.5,pcol = cor(data[,1],data[,2]))

fixed <- list(nugget_1 = 0, nugget_2 = 0,smooth_1 = 1.5, smooth_2 = 1.5, smooth_12 = 1.5, mean_1 = mean(data[,1]),
              mean_2 = mean(data[,2]))

data=t(data)

fitml <- GeoFit(data , coordx = coords , corrmodel = corrmodel ,
                likelihood = " Full", type = "Standard", start = start ,
                fixed = fixed , varest = TRUE)

fitml
