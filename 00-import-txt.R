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

datos_so <- datos |> filter(eg < 10000) 

plot(
  table(datos_so$eg),
    type = "h",
    lwd = 5,
    col = "darkgreen",
    main = "Frecuencia de cada valor de EG",
    xlab = "EG (conteos)", ylab = "Frecuencia"
)

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


# Modelos 
library(pscl)
library(MASS)
library(lmtest)
library(countreg)
