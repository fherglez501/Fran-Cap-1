#############
# Libraries #
#############
library(pscl)
library(MASS)
library(lmtest)
library(countreg)
#install.packages("countreg", repos="http://R-Forge.R-project.org") # Correlo desde aca porque no esta en CRAN

#########
# Datos #
#########
datos <- read.table("ParasiteCod.txt", h = T)
str(datos)
datos<-na.omit(datos) # sacando unos NA

datos$Area <- factor(datos$Area)
datos$Year <- factor(datos$Year)

############
# Analisis #
############
hist(datos$Intensity, ylab = "Frequencies", xlab = "Observed intensity values", breaks = 300) 
plot(table(datos$Intensity),ylab = "Frequencies", xlab = "Observed intensity values") 

# En los modelos inflados con 0s, lo que esta a la izquierda de la | muestra los predictores "tradicionales" (y se modela como conteo), y lo que está a la derecha muestra los predictores que inflan los 0s (se modela con binomial)

# Poisson #
zip1 <- zeroinfl(Intensity ~ Area * Year + Length | Area * Year + Length, dist = "poisson", link = "logit", data = datos)

# Overdispersion en ZIP? # Supuesto en conteos
E2 <- resid(zip1, type = "pearson")
N  <- nrow(datos)
p  <- length(coef(zip1))  
sum(E2^2) / (N - p) # Si, Mucho mayor que 1.2

# Binomial negativa #
Nb1 <- zeroinfl(Intensity ~ Area * Year + Length | Area * Year + Length, dist = "negbin", link = "logit", data = datos) 

# Overdispersion en Nb? #
E3 <- resid(Nb1, type = "pearson")
N  <- nrow(datos)
p  <- length(coef(Nb1)) + 1 # theta  
sum(E3^2) / (N - p) # 1.3

# Que modelo se ajusta mejor (i.e. básicamente para corregir por sobredispersion)
lrtest(zip1, Nb1) # Nos quedamos con Nb1
AIC(zip1, Nb1)

summary(Nb1)

###################
# Model selection #
###################
Nb1 <-    zeroinfl(Intensity ~ Area * Year + Length | Area * Year + Length, dist = "negbin", link = "logit", data = datos) # Modelo1 con binomial negativa
Nb1.LC <- zeroinfl(Intensity ~ Area * Year          | Area * Year + Length, dist = "negbin",link = "logit", data = datos) # Sin Length in count
lrtest(Nb1,Nb1.LC) # Sacar Length de los counts es significativo (i.e. NO lo sacamos)
AIC(Nb1,Nb1.LC)

Nb1.IntC<-zeroinfl(Intensity ~ Area + Year + Length | Area * Year + Length, dist = "negbin", link = "logit", data = datos) # Sin interaccion in count
lrtest(Nb1,Nb1.IntC) # No podemos sacar la interaccion

Nb1.LB <- zeroinfl(Intensity ~ Area * Year + Length | Area * Year         , dist = "negbin", link = "logit", data = datos) # Sin Length binomial
lrtest(Nb1,Nb1.LB) # No podemos sacar length de la parte binomial

Nb1.IntB<-zeroinfl(Intensity ~ Area * Year + Length | Area + Year + Length, dist = "negbin", link = "logit", data = datos) # Sin interaccion binomial
lrtest(Nb1,Nb1.IntB) # No podemos sacar la interaccion de la parte binomial

AIC(Nb1,Nb1.LC,Nb1.IntC,Nb1.LB,Nb1.IntB)

####################
# Model validation # Ya corregimos lo de la sobredispersion
####################
Nb1.residuals <- residuals(Nb1, type = "pearson")
Nb1.fitted <- fitted(Nb1)
plot(Nb1.residuals ~ Nb1.fitted) # No muy fantastico
plot(Intensity ~ Nb1.fitted, datos) # deberia ser linea de igualdad

rootogram(Nb1) # luce bien (explicacion: https://data.library.virginia.edu/simulating-data-for-count-models/)

##############
# Plot final # Hay que mejorarlo :)
##############

# COUNTS
new.data <- expand.grid(Area = levels(datos$Area), Year = levels(datos$Year), Length = seq(min(datos$Length), max(datos$Length), 1))
new.data$predicted <- predict(Nb1, newdata1, type ="count") # predecimos los conteos

# counts
ggplot(new.data, aes(x = Length, y = predicted, colour = factor(Year))) +
  geom_point() +
  geom_line() +
  facet_wrap(~Area, nrow=1) 
  

# ZERO
new.data2 <- expand.grid(Area = levels(datos$Area), Year = levels(datos$Year), Length = seq(min(datos$Length), max(datos$Length), 1))
new.data2$predicted <- predict(Nb1, newdata1, type ="zero") # predecimos los 0s

# counts
ggplot(new.data2, aes(x = Length, y = predicted, colour = factor(Year))) +
  geom_point() 

# Interpretacion
# https://math.usu.edu/jrstevens/biostat/projects2013/rep_ZIP.pdf
# https://timeseriesreasoning.com/contents/zero-inflated-poisson-regression-model/
# https://bookdown.org/roback/bookdown-BeyondMLR/ch-poissonreg.html#cs:drinking
# https://rstudio-pubs-static.s3.amazonaws.com/11785_e937a5f782864e1c9053d90b2b66c796.html
# https://stats.idre.ucla.edu/r/dae/zinb/