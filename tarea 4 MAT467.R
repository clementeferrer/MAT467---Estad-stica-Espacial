library(readxl)
library(geoR)
library(nlme)
library(mgcv)
library(sp)
library(nortest)
library(aod)
library(car)


##Importar datos

forest_excel <- read_excel("C:/Users/ccfer/Downloads/Datos forestales.xls")
df_forest = data.frame(forest_excel)
df_Y <- df_forest[ , c("ALMEDTOT")]
df_X <- df_forest[ , c("ABTOT", "EJE_X", "EJE_Y")]
df_XY <- df_forest[ , c("ALMEDTOT", "ABTOT", "EJE_X", "EJE_Y")]
df_arbol <- df_forest[ , c("ALMEDTOT", "ABTOT")]
df_coord <- df_forest[ , c("EJE_X", "EJE_Y")]
X =  as.matrix(df_X)
Y = as.matrix(df_Y)
XY = as.matrix(df_XY)
coord = as.matrix(df_coord)

## Ploteo, variograma y REML

geodata_Y <- as.geodata(XY, coords.col = 3:4, data.col = 1)
plot(geodata_Y)
vario <- variog(geodata_Y)
plot(vario)

reml <- likfit(geodata_Y, ini=c(50,3000), cov.model = "spherical", lik.met = "REML")
summary(reml)
plot(vario)
lines(reml, lty =2)

## Construccion de matriz de cov

mhw.sp <- df_arbol

coordinates(mhw.sp) <- df_coord
summary(mhw.sp)

D <- spDists(mhw.sp)
dim(D)

exp.cov <- function(h, c, a) {
  c * (1-(3/2)*a*h+(1/2)*(a*h)**3)
}
V <- 0.4838 + exp.cov(D, 0.7907, 3000)
V[1:5, 1:5]

## Betas

ajuste_err_corr <- gls(
  ALMEDTOT ~ 1 + ABTOT,
  data = df_XY,
  correlation = corSpher(form =  ~ EJE_X + EJE_Y),
  method = "REML")

summary(ajuste_err_corr)
coef(ajuste_err_corr)
vcov(ajuste_err_corr)

## Test de Wald y KR, junto a chequeo de supuesto

wald.test(Sigma = vcov(ajuste_err_corr), b = coef(ajuste_err_corr), Terms = 1:2)

hist(Y)
qqPlot(Y)
ad.test(Y)