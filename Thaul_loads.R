# Para modelos mixtos. Chequear esta libreria!!!!
# https://github.com/nyiuab/NBZIMM

#############
# Libraries #
#############
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggfortify)
library(MASS)
library(pscl)
library(lmtest) #(install.packages("lmtest", repos="http://R-Forge.R-project.org"))
library(countreg) #(install.packages("countreg", repos="http://R-Forge.R-project.org"))
library(boot)

#########
# Datos #
#########
datos_2 <- read.csv("Datos.csv", header = TRUE, stringsAsFactors = TRUE) # Saqué inviernos
str(datos_2)
summary(datos_2)

ggplot(datos_2, aes(loads))+
  geom_histogram(bins=50)+
  xlab("Intensidad de valores observados")+ ylab("Frecuencia")+
  theme_bw()

###############
# PCA - Clima #
###############
loads_pca <- prcomp(datos_2[,c(17:35)], scale =TRUE)
summary(loads_pca)
Loadings_loads <- loads_pca$rotation
Loadings_loads[,1:4]

datos_2 <- data.frame(PC1=predict(loads_pca)[,1], PC2=predict(loads_pca)[,2], PC3=predict(loads_pca)[,3], datos_2)
str(datos_2)

summary(lm(PC1 ~ tmean, datos_2)) # ~ 87%
summary(lm(PC2 ~ trange, datos_2)) # ~ 94% 
summary(lm(PC3 ~ isotherm, datos_2)) # ~ 94% 

##################
# Preparar datos #
##################
datos_2$loads<-as.integer(datos_2$loads) 

datos_2$Year <- scale(datos_2$Year)
datos_2$lat <- scale(datos_2$lat)
datos_2$alt <- scale(datos_2$alt)
datos_2$tmean <- scale(datos_2$tmean)
datos_2$trange <- scale(datos_2$trange)
datos_2$isotherm <- scale(datos_2$isotherm)
datos_2$Prevalence <- scale(datos_2$Prevalence)
datos_2$Bd_total <- scale(datos_2$Bd_total)
datos_2$age <- as.factor(datos_2$age)
datos_2$Season <- as.factor(datos_2$Season)
#############################
# Modelos I: Poisson or ZIP #
#############################

# Poisson #
zip_loads <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					  Year + Bd_total + Prevalence + lat + alt, 
				      dist = "poisson", link = "logit", data = datos_2)
				
# Overdispersion en ZIP? # Supuesto en conteos
E2 <- resid(zip_loads, type = "pearson")
N  <- nrow(datos_2)
p  <- length(coef(zip_loads))  
sum(E2^2) / (N - p) # Si, Mucho mayor que 1.2

# Binomial negativa #
Nb_loads <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					 Year + Bd_total + Prevalence + lat + alt, 
					 dist = "negbin", link = "logit", data = datos_2) 

# Overdispersion en Nb? #
E3 <- resid(Nb_loads, type = "pearson")
N  <- nrow(datos_2)
p  <- length(coef(Nb_loads)) + 1 # theta  
sum(E3^2) / (N - p) # 0.85

# Que modelo se ajusta mejor (i.e. básicamente para corregir por sobredispersion)
lrtest(zip_loads, Nb_loads) # Nos quedamos con Nb
AIC(zip_loads, Nb_loads) # Nos quedamos con Nb

###############################
# Modelos II: Model selection #
###############################

# Binomial #
Nb_loads <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					 Year + Bd_total + Prevalence + lat + alt, 
					 dist = "negbin", link = "logit", data = datos_2) 

Nb_loads_Year <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					       Bd_total + Prevalence + lat + alt, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_Bd_total <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					 Year       + Prevalence + lat + alt, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_Prevalence <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					 Year + Bd_total +      lat + alt, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_lat <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					 Year + Bd_total + Prevalence +  alt, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_alt <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | 
					 Year + Bd_total + Prevalence + lat, 
					 dist = "negbin", link = "logit", data = datos_2) 

lrtest(Nb_loads, Nb_loads_Year) # NS
lrtest(Nb_loads, Nb_loads_Bd_total) # NS
lrtest(Nb_loads, Nb_loads_Prevalence) # Significativo
lrtest(Nb_loads, Nb_loads_lat) # NS
lrtest(Nb_loads, Nb_loads_alt) # NS

# Conteos #
Nb_loads_2 <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 

Nb_loads_2_age      <- zeroinfl(loads ~       Year + Season + lat + alt + tmean + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_Year     <- zeroinfl(loads ~ age +     Season + lat + alt + tmean + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_Season   <- zeroinfl(loads ~ age + Year +       lat + alt + tmean + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_lat      <- zeroinfl(loads ~ age + Year + Season +    alt + tmean + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_alt      <- zeroinfl(loads ~ age + Year + Season + lat +    tmean + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_tmean    <- zeroinfl(loads ~ age + Year + Season + lat + alt      + trange + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_trange   <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean       + isotherm | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 
Nb_loads_2_isotherm <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange         | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 


lrtest(Nb_loads_2, Nb_loads_2_age) # Significative
lrtest(Nb_loads_2, Nb_loads_2_Year) # Significative
lrtest(Nb_loads_2, Nb_loads_2_Season) # Significative
lrtest(Nb_loads_2, Nb_loads_2_lat) # Significative
lrtest(Nb_loads_2, Nb_loads_2_alt) # Significative
lrtest(Nb_loads_2, Nb_loads_2_tmean) # Significative
lrtest(Nb_loads_2, Nb_loads_2_trange) # Significative
lrtest(Nb_loads_2, Nb_loads_2_isotherm) # NS

# Modelo final #
Nb_loads_final <- zeroinfl(loads ~ age + Year + Season + lat + alt + tmean + trange | Prevalence, 
					 dist = "negbin", link = "logit", data = datos_2) 

#############
# Supuestos # 
#############
# Nb_loads_final.residuals <- residuals(Nb_loads_final, type = "pearson")
# Nb_loads_final.fitted <- fitted(Nb_loads_final)
# plot(Nb_loads_final.residuals ~ Nb_loads_final.fitted) # No muy fantastico
# plot(loads ~ Nb_loads_final.fitted, datos_2) # deberia ser linea de igualdad

countreg::rootogram(Nb_loads_final) # luce bien (explicacion: https://data.library.virginia.edu/simulating-data-for-count-models/)

summary(Nb_loads_final)

###################
# Figuras finales # 
###################

# COUNTS
new.data <- expand.grid(age = levels(datos_2$age), 
						Year = seq(min(datos_2$Year), max(datos_2$Year), 0.75),
						Season = levels(datos_2$Season), 
						lat = seq(min(datos_2$lat), max(datos_2$lat), 0.75),
						alt = seq(min(datos_2$alt), max(datos_2$alt), 0.75),
						tmean = seq(min(datos_2$tmean), max(datos_2$tmean), 0.75),
   					    trange = seq(min(datos_2$trange), max(datos_2$trange), 0.75),
   					    Prevalence = seq(min(datos_2$Prevalence), max(datos_2$Prevalence), 0.5))

new.data$predicted <- predict(Nb_loads_final, new.data, type ="count") # predecimos los conteos


plot.age <- ggplot(new.data, aes(x = age, y = predicted, colour = age)) +
  geom_boxplot() +
  xlab("Age")+
  ylab("Predicted counts")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))
  
plot.season <- ggplot(new.data, aes(x = Season, y = predicted, colour = Season)) +
  geom_boxplot() +
  xlab("Season")+
  ylab("Predicted counts")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))
        
plot.lat <- ggplot(new.data, aes(x = lat, y = predicted)) +
  geom_point(size = 2.5, alpha = 0.1, colour = "red") +
  xlab("Latitude")+
  ylab("Predicted counts")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))

plot.alt <- ggplot(new.data, aes(x = alt, y = predicted)) +
  geom_point(size = 2.5, alpha = 0.1, colour = "darkgreen") +
  xlab("Altitude (m)")+
  ylab("Predicted counts")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))

plot.tmean <- ggplot(new.data, aes(x = tmean, y = predicted)) +
  geom_point(size = 2.5, alpha = 0.1, colour = "orange") +
  xlab("Mean temperature (ºC)")+
  ylab("Predicted counts")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))

plot.trange <- ggplot(new.data, aes(x = trange, y = predicted)) +
  geom_point(size = 2.5, alpha = 0.1, colour = "blue") +
  xlab("Range temperature (ºC)")+
  ylab("Predicted counts")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))
  
grid.arrange(plot.age, plot.season, plot.lat, plot.alt, plot.tmean, plot.trange, nrow = 3)


# ZEROs
new.data$predicted_zero <- predict(Nb_loads_final, new.data, type ="zero") # predecimos los 0s

ggplot(new.data, aes(x = Prevalence, y = predicted_zero)) +
  geom_point(size = 2.5, alpha = 0.1, colour = "red") +
  xlab("Prevalence")+
  ylab("Predicted zeros")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, size = 1.5))

# BOOT ----
db_Nb_loads_final <- datos_2[,c("age", "Year", "Season", "lat",
                        "alt", "tmean", "trange", "Prevalence", "loads")]
Nb_loads_final <- zeroinfl(loads ~ age + Year + Season + lat + alt + 
          tmean + trange | Prevalence, dist = "negbin", link = "logit", 
          data = db_Nb_loads_final) 


stat <- function(df, inds, new.dat) {
  model <- formula(loads ~ age + Year + Season + lat + alt + tmean + trange | Prevalence)
  predict(
    zeroinfl(model, dist = "negbin", link = "logit", data = df[inds, ]),
    newdata = new.dat)
}


new.data <- expand.grid(age = levels(datos_2$age), 
                        Year = seq(min(datos_2$Year), max(datos_2$Year), 0.75),
                        Season = levels(datos_2$Season), 
                        lat = 0,
                        alt = 0,
                        tmean = seq(min(datos_2$tmean), max(datos_2$tmean), 0.75),
                        trange = seq(min(datos_2$trange), max(datos_2$trange), 0.75),
                        Prevalence = seq(min(datos_2$Prevalence), max(datos_2$Prevalence), 0.5))


set.seed(1234)
res <- boot(db_Nb_loads_final, stat, R = 200)

all.equal(res$t0, predict(Nb_loads_final)) ### TRUE


CI <- setNames(as.data.frame(t(sapply(1:nrow(db_Nb_loads_final), function(row)
  boot.ci(res, conf = 0.95, type = "basic", index = row)$basic[, 4:5]))),
  c("lower", "upper"))
CI <- data.frame(db_Nb_loads_final, Predicted=res$t0, CI)

