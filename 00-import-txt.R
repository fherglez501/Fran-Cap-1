# library ----
 pacman::p_load(geodata, terra, tidyverse, readr, pscl, MASS, lmtest, countreg)

 # data ----
 datos <- readxl::read_excel(here::here("Bd_tesis_FC_2025.xlsx"), sheet = 1)
 
 # functions ----
 dms_to_decimal <- function(coord, negativo = TRUE) {
   # Extrae grados, minutos y segundos con expresiones regulares
   partes <- stringr::str_match(coord, "([0-9]+)[°º]_?([0-9]+)'_?([0-9\\.]+)")
   grados <- as.numeric(partes[, 2])
   minutos <- as.numeric(partes[, 3])
   segundos <- as.numeric(partes[, 4])
   decimal <- grados + minutos / 60 + segundos / 3600
   # Aplica el signo negativo si corresponde
   decimal <- if (negativo) -decimal else decimal
   
   # Redondea a 6 decimales
   round(decimal, 6)
 }

# Procesamiento de datos
datos <- datos |> 
  # mutate(DATE = ymd(DATE)) |> 
  mutate(DATE = as.Date(DATE, origin = "1899-12-30")) |> 
  # Paso 1: convertir coordenadas y EG
  mutate(
    lat = dms_to_decimal(Coordinate_S, negativo = TRUE),
    lon = dms_to_decimal(Coordinate_W, negativo = TRUE),
    EG = round(EG)
  ) |> 
  # convertir caracteres a factor
  mutate(across(where(is.character), as.factor)) %>% 
  # extraer bioclimáticas de WorldClim
  bind_cols(
    terra::extract(
      geodata::worldclim_global(var = "bio", res = 2.5, path = tempdir()),
      terra::vect(.[, c("lon", "lat")], crs = "EPSG:4326")
    ) |> 
      tibble::as_tibble() |> 
      dplyr::select(-ID)
  ) |> 
  janitor::clean_names() |> 
  janitor::remove_empty(c("rows", "cols"))
        
# Cambiar fecha, peso y longitud
datos <- datos |>
  mutate(
    across(c(length_mm, weight_gr), as.numeric)
  ) |> 
  mutate(
    age = case_when(
      age %in% c("Adulto", "Adulto_") ~ "Adulto",
      age == "Juvenil"                ~ "Juvenil",
      age == "Pre-adulto"             ~ "Pre-adulto",
      TRUE                            ~ NA_character_  # para valores no esperados
    )
  )
  
summary(datos$eg)

# plot de fechas
plot(table(as.factor(datos$date)))

datos_so <- datos |> 
  filter(eg < 10000) |> 
  mutate(
    n = n(),
    pos = case_when(
      eg > 0 ~ 1,
      eg == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) |> 
  group_by(location) |> 
  mutate(prev = mean(pos, na.rm = TRUE)) |> 
  ungroup()

plot(
  table(datos_so$eg),
    type = "h",
    lwd = 5,
    col = "darkgreen",
    main = "Frecuencia de cada valor de EG",
    xlab = "EG (conteos)", ylab = "Frecuencia"
)


# Grafica de Prevalencia de Bd por localidad
ggplot(datos_so, aes(x = reorder(location, lat), y = prev)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Opcional para una visualización horizontal
  labs(
    title = "Prevalencia de Bd por localidad",
    x = "Localidad",
    y = "Prevalencia (proporción de positivos)"
  ) +
  theme_classic()


# Lambda: el promedio de EG para todos los animales
lambda_est <- mean(datos_so$eg, na.rm = TRUE)
lambda_est

p0_expected <- dpois(0, lambda_est)
p0_expected

p0_observed <- mean(datos_so$eg == 0, na.rm = TRUE)
p0_observed

cat("Proporción esperada de ceros bajo Poisson:", round(p0_expected, 3), "\n")
cat("Proporción observada de ceros en tus datos:", round(p0_observed, 3), "\n")

# test formal con DHARMa
library(DHARMa)

modelo_pois <- glm(eg ~ 1, data = datos, family = poisson)
res <- simulateResiduals(modelo_pois)
testZeroInflation(res)

# Regla de Vredemburg: ZEG 100, 1000, 10000

# Tus datos tienen demasiados ceros para ser explicados por un modelo Poisson.
# Necesitas considerar un modelo que permita más ceros, como:
# Modelo Poisson inflado en ceros (ZIP)
# Modelo binomial negativo (NB)
# Modelo binomial negativo inflado en ceros (ZINB) si también hay sobre-dispersión.

# PCA para reducir el numero de variables bioclimáticas
bio_vars <- datos_so[, grep("wc2_1_2_5m_bio_", names(datos_so))]

# Eliminar filas con NA (opcional según tu análisis)
bio_vars <- na.omit(bio_vars)

# PCA
pca_bio <- prcomp(bio_vars, scale. = TRUE)

# Calcular la varianza acumulada
var_exp <- summary(pca_bio)$importance["Cumulative Proportion", ]

# Determinar cuántos componentes explican al menos 80% de la varianza
num_pc <- which(var_exp >= 0.80)[1]

# imprimir en la consola
cat("Se conservarán", num_pc, "componentes principales que explican el",
    round(var_exp[num_pc] * 100, 2), "% de la varianza total.\n")

# Resumen PCA
summary(pca_bio)
# plot PCA
plot(pca_bio, type = "l")  # Scree plot

# Crear nuevo dataframe con los componentes seleccionados
pc_scores_chat <- as.data.frame(pca_bio$x[, 1:num_pc])

# PCA como nuevas variables
pc_scroes <- as.data.frame(pca_bio$x)

# marco de datos con los valores del PCA (cargas)
datos_so_pca <- cbind(datos_so, pc_scroes[, 1:3])

# cargas del PCA
loadings <- pca_bio$rotation[, 1:num_pc]  # Solo para los componentes retenidos

# ordenarlos para ver cuáles dominan cada componente:
apply(loadings, 2, function(x) sort(abs(x), decreasing = TRUE))

# Variables con carga absoluta > 0.6 en algún componente retenido
important_vars <- rownames(loadings)[apply(abs(loadings), 1, max) > 0.3]
print(important_vars)

# Modelos ZIP
# En los modelos inflados con 0s, lo que esta a la izquierda de la | muestra los predictores "tradicionales" (y se modela como conteo), y lo que está a la derecha muestra los predictores que inflan los 0s (se modela con binomial)

m0 <- zeroinfl(
  eg ~ age + lat + altitude +
    wc2_1_2_5m_bio_1 + wc2_1_2_5m_bio_2 + wc2_1_2_5m_bio_3 +
    wc2_1_2_5m_bio_4 + wc2_1_2_5m_bio_5 + wc2_1_2_5m_bio_6 +
    wc2_1_2_5m_bio_7 + wc2_1_2_5m_bio_8 + wc2_1_2_5m_bio_9 +
    wc2_1_2_5m_bio_10 + wc2_1_2_5m_bio_11 + wc2_1_2_5m_bio_12 +
    wc2_1_2_5m_bio_13 + wc2_1_2_5m_bio_14 + wc2_1_2_5m_bio_15 +
    wc2_1_2_5m_bio_16 + wc2_1_2_5m_bio_17 + wc2_1_2_5m_bio_18 +
    wc2_1_2_5m_bio_19,
  dist = "poisson", link = "logit", data = datos_so
)

m1 <- zeroinfl(
  eg ~ age + lat + altitude +
    wc2_1_2_5m_bio_3 + wc2_1_2_5m_bio_4 + wc2_1_2_5m_bio_6 +
    wc2_1_2_5m_bio_7 + wc2_1_2_5m_bio_9 + wc2_1_2_5m_bio_11,
  dist = "poisson", link = "logit", data = datos_so
)

m2 <- zeroinfl(
  eg ~ age + lat + altitude +
    wc2_1_2_5m_bio_3 + wc2_1_2_5m_bio_6 +
    wc2_1_2_5m_bio_7 + wc2_1_2_5m_bio_9 + wc2_1_2_5m_bio_11 | age + lat + altitude,
  dist = "poisson", link = "logit", data = datos_so
)

# 
library(glmmTMB)

m2_glmm <- glmmTMB(
  eg ~ age + lat + altitude +
    wc2_1_2_5m_bio_3 + wc2_1_2_5m_bio_6 +
    wc2_1_2_5m_bio_7 + wc2_1_2_5m_bio_9 + wc2_1_2_5m_bio_11 +
    (1 | location), # efecto aleatorio = localidad
  ziformula = ~ age + lat + altitude, # estos son los predictores que inflan los ceros
  family = poisson(),
  data = datos_so
)

performance::check_overdispersion(m2_glmm)

# Asegúrate de escalar peso y longitud para estabilidad numérica
datos_so$peso_z <- scale(datos_so$weight_gr)
datos_so$longitud_z <- scale(datos_so$length_mm)

# Modelo con covariables individuales
m3_glmm <- glmmTMB(
  eg ~ age + lat + altitude + peso_z + longitud_z +
    wc2_1_2_5m_bio_3 + 
    #wc2_1_2_5m_bio_6 +
    wc2_1_2_5m_bio_7 + 
    #wc2_1_2_5m_bio_9 + 
    wc2_1_2_5m_bio_11 +
    (1 | location),
  ziformula = ~ age + lat + altitude + prev + peso_z + longitud_z,
  family = poisson(),
  data = datos_so
)

performance::check_overdispersion(m3_glmm)