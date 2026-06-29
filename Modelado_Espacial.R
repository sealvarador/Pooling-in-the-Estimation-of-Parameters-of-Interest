################################################################################
# Modelado_Espacial.R
#
# Objetivo:
# Ajustar y comparar modelos de estimación para áreas pequeñas
# bajo un enfoque frecuentista espacial utilizando los estratos 
# construidos a partir de variables sociodemográficas.
#
# El script realiza:
#   1. Construcción de estimadores directos y varianzas de muestreo.
#   2. Ajuste del modelo SAR Frecuentista.
#   3. Ajuste del modelo SAR Bayesiano.
#   4. Diagnóstico y evaluación de supuestos.
#   5. Obtención de predicciones SAE por estrato.
#   6. Visualización y almacenamiento de resultados.
#   7. Exportación de resultados para análisis posteriores.
#
# Autor: Selene Alvarado Rosario
# Fecha: Mayo 2025
################################################################################


# Ajustar ruta de trabajo según la ubicación del proyecto


library(dplyr)
library(spdep)
library(ggplot2)
library(Matrix)
library(MASS)
library(fastDummies)
library(sae)
library(patchwork) 
library(saeHB.spatial)
library(coda)
library(mvtnorm)
library(forcats)


df_areas <- read.csv("complete_birth_rec.csv", stringsAsFactors = FALSE)        #Importo el dataset que se había preparado, ya con ethnic_rec

df_areas$type_area <- factor(
  df_areas$type_area,
  levels = c("rural", "urbana")
)
df_areas$high_edu <- factor(
  df_areas$high_edu,
  levels = c("sin educación", "primaria", "secundaria", "superior")
)
df_areas$ethnic_rec <- factor(
  df_areas$ethnic_rec,
  levels = c("otros", "indigena", "blanco/mestizo", "afro")
)
df_areas$cohort_quinquenal <- factor(
  df_areas$cohort_quinquenal,
  levels = c("1965–1969", "1970–1974", "1975–1979", "1980–1984",
             "1985–1989", "1990–1994", "1995–2002" )
)


strata_obs <- df_areas %>%
  group_by(type_area, high_edu, ethnic_rec, cohort_quinquenal) %>%
  summarise(
    n_h = n(),
    y_h = mean(age_1st_b, na.rm = TRUE),
    var_h = var(age_1st_b, na.rm = TRUE) / n_h,
    .groups = "drop"
  )                                                                             #resumen de los estratos observados


strata_obs <- strata_obs %>%
  mutate(stratum_id = row_number())                                             #Incluyo ID por estrato

strata_obs$strata <- interaction(
  strata_obs$type_area,
  strata_obs$high_edu,
  strata_obs$ethnic_rec,
  strata_obs$cohort_quinquenal,
  sep = "."
)


df_indiv <- df_areas %>%                                                        #incluyo la información a cada individuo del estrato al que pertenece
  left_join(
    strata_obs %>%
      dplyr::select(
        type_area,
        high_edu,
        ethnic_rec,
        cohort_quinquenal,
        stratum_id,
        strata
      ),
    by = c(
      "type_area",
      "high_edu",
      "ethnic_rec",
      "cohort_quinquenal"
    )
  )


strats <- df_indiv %>%
  dplyr::select(strata, stratum_id) %>%
  distinct(strata, .keep_all = TRUE) %>%
  mutate(stratum_id = as.character(stratum_id))


### Varianza global de respaldo ###
# Se utiliza para imputar la varianza de muestreo en estratos con n_i = 1, donde
# la varianza muestral no puede calcularse.
S_global <- as.numeric(var(log(df_indiv$age_1st_b), na.rm = TRUE))


# Se calcula el estimador directo de la media logarítmica, la varianza de 
# muestreo aproximada y el tamaño muestral.
df_small_area <- df_indiv %>%
  group_by(strata, cohort_quinquenal, high_edu, type_area, ethnic_rec) %>%
  summarize(
    y_dir = mean(log(age_1st_b), na.rm = TRUE),
    # Varianza de muestreo aproximada para la media
    psi_i = ifelse(n() > 1, 
                   var(log(age_1st_b), na.rm = TRUE) / n(), 
                   S_global / n()),
    n_i = n(),
    .groups = "drop"
  ) %>%
  mutate(
    y_dir = as.numeric(y_dir),
    psi_i = as.numeric(psi_i)
  )

df_small_area <- as.data.frame(df_small_area)

sum(is.na(df_small_area$psi_i))  #0
sum(is.na(df_small_area$y_dir))  #0


##### construcción de paletas #####

grunge_fun <- colorRampPalette(
  c("#2F3E2E",  # verde bosque oscuro
    "#3A5A40",  # verde musgo
    "#588157",  # sage
    "#A3B18A",  # verde grisáceo
    "#B7B7A4",  # gris oliva
    "#7F5539",  # café tierra
    "#5E503F")  # marrón oscuro
)


weights_fun <- colorRampPalette(
  c(#"#2F3E2E",  # verde bosque oscuro
    "#3A5A40",  # verde musgo
    #"#588157",  # sage
    #"#A3B18A",  # verde grisáceo
    "#B7B7A4",  # gris oliva
    "#7F5539",  # café tierra
    #"#5E503F",  # marrón oscuro
    "#1D0B14")  
)

weights_cols <- weights_fun(100)


weights_fun_hb <- colorRampPalette(c("#3A5A40", "#B7B7A4", "#7F5539", "#1D0B14"))
weights_cols_hb <- weights_fun_hb(100)


##### Matriz W de vecinos en el espacio de atributos #####
# función para obtener distancias de Mahalanobis entre estratos
pairwise_maha <- function(mat, S) {
  l <- nrow(mat); M <- matrix(0, l, l)
  for(i in 1:l){
    for(j in i:l){
      d2 <- as.numeric((mat[i,]-mat[j,]) %*% MASS::ginv(S) %*% (mat[i,]-mat[j,]))
      M[i,j] <- d2; M[j,i] <- d2
    }
  }
  sqrt(M)
}


# función para optimizar d0
find_best_d0 <- function(maha_dist, y_h, l) {
  y_centered <- y_h - mean(y_h)
  denom <- sum(y_centered^2)
  
  # distancias únicas
  d_vals <- sort(unique(maha_dist[lower.tri(maha_dist)]))
  
  best_I <- -Inf
  best_d0 <- NA
  
  for (d0 in d_vals) {
    # Matriz de pesos binaria
    W <- (maha_dist <= d0) * 1
    diag(W) <- 0  # asegúrate de que w_ii = 0
    
    # Numerador
    num <- sum(W * (y_centered %o% y_centered))
    # Suma total de pesos
    w_sum <- sum(W)
    
    # Moran’s I
    I <- (l / w_sum) * (num / denom)
    
    # Actualiza el máximo
    if (I > best_I) {
      best_I <- I
      best_d0 <- d0
    }
  }
  
  moran_vals <- sapply(sort(unique(maha_dist[lower.tri(maha_dist)])), function(d0) {
    W <- (maha_dist <= d0) * 1
    diag(W) <- 0
    num <- sum(W * (y_centered %o% y_centered))
    w_sum <- sum(W)
    (length(y_h) / w_sum) * (num / denom)
  })
  
  plot(sort(unique(maha_dist[lower.tri(maha_dist)])), moran_vals, type = "l",
       xlab = "d0 (umbral de distancia de Mahalanobis)",
       ylab = "Índice de Moran I",
       main = "Moran I vs d0")
  
  list(best_d0 = best_d0, best_I = best_I)
}


# función que crea la matriz W de pesos para el umbral d0
matriz_pesos_fun <- function(maha_dist, d0) {
  # Si la distancia es menor o igual a d0, es 1 (vecino), si no, 0.
  W_binaria <- (maha_dist <= d0)
  diag(W_binaria) <- 0 
  return(W_binaria)
}


### Construcción de la matriz de covariables ### 
X_preparada <- strata_obs %>%
  
  # Cohorte quinquenalm es ordinal (Punto medio para distancia real - según de Lockwood y Nandram, 2024)
  mutate(
    cohort_quinquenal_i = case_when(
      cohort_quinquenal == "1965–1969" ~ 1967,
      cohort_quinquenal == "1970–1974" ~ 1972,
      cohort_quinquenal == "1975–1979" ~ 1977,
      cohort_quinquenal == "1980–1984" ~ 1982,
      cohort_quinquenal == "1985–1989" ~ 1987,
      cohort_quinquenal == "1990–1994" ~ 1992,
      cohort_quinquenal == "1995–2002" ~ 1998.5
    ),
    
    # Nivel educativo es ordinal
    high_edu_i = as.integer(factor(high_edu, 
                                   levels = c("sin educación", "primaria", "secundaria", "superior"), 
                                   ordered = TRUE))
  ) %>%
  
  # Generar dummies para variables nominales
  dummy_cols(select_columns = c("type_area", "ethnic_rec"), 
             remove_first_dummy = FALSE) %>% 
  
  # se incluyen todas las categorías de etnia y área para que la distancia sea
  # consistente, mas adelante quitamos la categoría base de cada una
  dplyr::select(
    cohort_quinquenal_i, 
    high_edu_i,
    starts_with("type_area_"), 
    starts_with("ethnic_rec_")
  )

#estadísticos importantes
y_s <- log(strata_obs$y_h)                                                      #variable respuesta
X <- as.matrix(X_preparada)                                                     #matriz de covariables
l <- length(X[,1])                                                              #número de estratos
S <- cov(X)
maha_dist <- matrix(0, l, l)

# preparación de información para la fórmula, unir y_s con X_preparada
X_preparada <- X_preparada %>%
  dplyr::select(-ethnic_rec_otros, -type_area_rural)

# X para Mahalanobis conserva todas las categorías
# Posteriormente se eliminan categorías base únicamente para el ajuste SAR


### Construcción de la matriz de vecindad en el espacio de atributos ###
maha_dist <- pairwise_maha(X, S)
d0 <- find_best_d0(maha_dist, y_s, l)

I_m <- d0$best_I
d0 <-  d0$best_d0

W <- matriz_pesos_fun(maha_dist, d0) 
sum(W)
n_i <- strata_obs$n_h

deg <- rowSums(W)
which(deg == 0)

lonely_strats <- strata_obs[which(deg == 0),]
lonely_strats                                                                   #solo 1 estrato no tiene vecinos, no se modifica la estructura espacial

# estandarización de W
W_listw <- mat2listw(W, style = "W", zero.policy = TRUE)                        #"W" para estandarización por filas
as.matrix(W_listw$weights)

W_est <- listw2mat(W_listw)



### Ajuste del modelo Fay-Herriot espacial o SAR frecuentista ###
# definición de datos para el modelo
datos_modelo <- data.frame(
  y_log = df_small_area$y_dir,
  vardir = df_small_area$psi_i,
  X_preparada
)


#Implementación del modelo
mod_SARf <- sae::mseSFH(
  formula = y_log ~ cohort_quinquenal_i + high_edu_i + type_area_urbana +
    ethnic_rec_indigena + ethnic_rec_blanco.mestizo + ethnic_rec_afro,
  vardir = vardir,
  proxmat = W_est,        
  method = "REML",                                                              #Método de Máxima Verosimilitud Restringida 
  data = datos_modelo
)



### Diagnóstico y verificación de supuestos ###
# Diagnóstico gráfico de residuos y efectos aleatorios
est_sar <- mod_SARf$est$eblup                                                   #estimaciones espaciales
residuos <- datos_modelo$y_log - est_sar
rho_est <- mod_SARf$est$fit$spatialcorr                                         #Coeficiente de autocorrelación espacial


# Prueba de normalidad de Kolmogorov-Smirnov
residuos_std <- (residuos - mean(residuos)) / sd(residuos)
test_norm <- ks.test(residuos_std, "pnorm")
print(test_norm)

# Histogramas
hist(residuos_std, main = "Distribución de Residuos", xlab = "Residuos de Pearson", col = "lightblue")

# Gráficos Q-Q norm
qqnorm(residuos_std)
qqline(residuos_std, col = "red")

# Residuos vs. Valores Ajustados
plot(est_sar, residuos, main = "Residuos vs. Ajustados", 
     xlab = "Valores Predichos", ylab = "Residuos")
abline(h = 0, col = "red")


# Prueba de autocorrelación espacial
test_moran <- moran.test(residuos, W_listw, alternative = "two.sided")
print(test_moran)                                                               #p-valor > alpha, la estructura espacial ha sido capturada por el modelo



### Cálculo del factor de encogimiento ###

A_est <- mod_SARf$est$fit$refvar                                                #Varianza de efectos aleatorios
I_mat <- diag(l)

# construcción de la matriz de covarianza de los efectos aleatorios SAR
# G = A * [(I - rho*W_est)' * (I - rho*W_est)]^-1

C <- I_mat - rho_est * W_est
G <- A_est * solve(t(C) %*% C)                                                  #Matriz de covarianza espacial 


# Matriz de varianza de muestreo (Psi) y varianza total (Sigma)
psi <- diag(datos_modelo$vardir)                                                #Varianzas directas conocidas 
Sigma <- G + psi                                                                #Varianza total del modelo


# construcción de la Matriz de Pesos EBLUP 
# K = G * Sigma^-1
K <- G %*% solve(Sigma)


# Crear dataframe con los 3 componentes de peso
SAR_pesos <- data.frame(
  n_i = n_i,                                                                    # Tamaño de muestra por área
  Peso_Propio = diag(K),
  Peso_Vecindario = rowSums(K) - diag(K),
  Peso_Sintetico = 1 - rowSums(K)
) %>%
  mutate(stratum_id = row_number()) 

SAR_pesos$num_vecinos <- rowSums(W)
table(SAR_pesos$num_vecinos)


# creación de grupo facet según la cantidad de vecinos
SAR_pesos <- SAR_pesos %>%
  mutate(grupo_facet = case_when(
    num_vecinos %in% c(0, 1) ~ "0 o 1 Vecino (Baja conectividad)",
    num_vecinos == 2        ~ "2 Vecinos (Alta conectividad)"
  ),
  # Asegurar el orden de los niveles
  vecinos_cat = factor(num_vecinos, levels = c(0, 1, 2),
                       labels = c("0 vecinos", "1 vecino", "2 vecinos")))


### Visualización de pesos ###
#vecindario vs. global
p1 <- ggplot(SAR_pesos, aes(x = Peso_Vecindario, 
                            y = Peso_Sintetico, 
                            color = n_i,
                            alpha = n_i
)) +
  # Aumentamos un poco el alpha para que las estrellas sean legibles
  geom_point(size = 4) +
  facet_wrap(~ grupo_facet)+
  scale_color_gradientn(colors = weights_cols, limits = c(min(SAR_pesos$n_i), max(SAR_pesos$n_i))) +
  scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") + 
  
  labs(
    #title = "Análisis de Pesos Sintético vs. vecindario",
    #subtitle = "Color: Tamaño de muestra (n_i)",
    x = "Peso vecindario",
    y = "Peso de la Media Global",
    color = expression(n[i])
  ) +
  theme_minimal() +
  theme(legend.position = "right")

#propio vs. global
p2 <- ggplot(SAR_pesos, aes(x = Peso_Propio, 
                            y = Peso_Sintetico, 
                            color = n_i,
                            alpha = n_i
)) +
  geom_point(size = 4) +
  facet_wrap(~ grupo_facet)+
  scale_color_gradientn(colors = weights_cols, limits = c(min(SAR_pesos$n_i), max(SAR_pesos$n_i))) +
  scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") + 
  
  labs(
    x = "Peso propio",
    y = "Peso de la Media Global",
    color = expression(n[i]) 
  ) +
  theme_minimal() +
  theme(legend.position = "right")


(p1 / p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Diagnóstico de Factores de Contracción del Modelo Espacial",
    subtitle = "Comparativa según conectividad y tamaño de muestra"
  )



### Coeficientes de Variación ###
# Se extraen los errores cuadráticos medios (MSE) y medidas relativas de precisión 
#para cada estrato - CV calculado sobre la escala logarítmica del modelo
mse_sar <- mod_SARf$mse

# Cálculo del CV
RSE_log_sar <- sqrt(mse_sar) / abs(est_sar)


resultados_estratos_SARf <- strata_obs %>%
  mutate(
    edad_estimada = est_sar,
    error_estandar = sqrt(mse_sar),
    RSE_log = RSE_log_sar * 100
  )

resultados_estratos_SARf <- resultados_estratos_SARf %>%
  left_join(strats, by = "strata") %>%
  arrange(n_h) %>%
  mutate(order_n = row_number())


ggplot(resultados_estratos_SARf, aes(x = order_n, y = RSE_log)) +
  geom_point() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Error Relativo Estándar por Estrato", 
       x = "Estratos (ordenados por tamaño)", y = "RSE (%)") +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )



### Criterios de Información ###
k <- length(mod_SARf$est$fit$estcoef$beta)+2                                    #Número de parámetros
log_L <- as.numeric(mod_SARf$est$fit$goodness[1]) 
log_L

KIC <- -2*log_L + 3*k
KIC

mod_SARf$est$fit$goodness


##### Visualización de medias #####
# df para predicciones
df_predicciones_SARf <- df_indiv %>%
  group_by(stratum_id, strata, cohort_quinquenal, high_edu, type_area, ethnic_rec) %>%
  summarise(n = n(), .groups = "drop")

# predicción en escala log
df_predicciones_SARf$pred_log <- est_sar

# volver a la escala real (años de edad)
df_predicciones_SARf <- df_predicciones_SARf %>%
  mutate(edad_real = exp(pred_log))


### Mapa de calor - combinaciones de variables ###
ggplot(df_predicciones_SARf, aes(x = cohort_quinquenal, y = high_edu, fill = edad_real)) +
  
  geom_tile(color = "white") +
  
  # divide matriz por área (filas) y Etnia (columnas)
  facet_grid(type_area ~ ethnic_rec) + 
  
  scale_fill_gradientn(
    colours = grunge_fun(100),
    name = "Edad \nPrimer \nHijo",
    limits = c(min(df_predicciones_SARf$edad_real),
               max(df_predicciones_SARf$edad_real)),
    oob = scales::squish
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Cohorte de Nacimiento",
    y = "Nivel Educativo") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )


### diagrama de dispersión de la edad promedio predicha por estrato - ordenado 
#por tamaño ###
df_plot_SARf <- df_predicciones_SARf %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n) %>%
  mutate(order_n = row_number())


ggplot(df_plot_SARf) +
  geom_point(aes(x = order_n, y = edad_real),
             color = "black",
             shape = "✦",
             alpha = 0.6,      
             size = 5) +     
  labs(x = "Estratos (ordenados por tamaño)",
       y = "Mediana de Edad Estimada") +
  theme(
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)  
  )


### Exportación de resultados ###
# Se almacenan los efectos aleatorios estimados y las medidas de precisión para su 
#uso posterior en comparación de modelos.
write.csv(resultados_estratos_SARf, "resultados_estratos_SARf.csv", row.names = FALSE)



# ------------------------------------------------------------------------------- #
##### Ajuste del modelo Fay-Herriot espacial-HB o SAR bayesiano jerárquico #####

set.seed(151724)  #Semilla para reproducibilidad
formula_sae <- y_log ~ cohort_quinquenal_i + high_edu_i + type_area_urbana +
  ethnic_rec_indigena + ethnic_rec_blanco.mestizo + ethnic_rec_afro

modelo_sar_hb <- sar.normal (
  formula     = formula_sae,
  vardir      = "vardir",                                                        
  proxmat     = W_est,                 
  iter.mcmc   = 100000,
  burn.in     = 20000,
  thin        = 20,
  iter.update = 10,
  data        = datos_modelo
)

autocorr(modelo_sar_hb$plot[[3]])
x11()
plot(modelo_sar_hb$plot[[3]])


### Diagnóstico y verificación de supuestos ###
# estadísticos puntuales de rho
rho_mean   <- modelo_sar_hb$coefficient["rho", "Mean"]
rho_sd     <- modelo_sar_hb$coefficient["rho", "SD"]
rho_median <- modelo_sar_hb$coefficient["rho", "50%"]

# límites del Intervalo de Credibilidad al 95%
rho_ic_inf <- modelo_sar_hb$coefficient["rho", "2.5%"]
rho_ic_sup <- modelo_sar_hb$coefficient["rho", "97.5%"]

cat("--- Extracción de Parámetro Espacial Rho ---\n",
    "Media Posterior:     ", rho_mean, "\n",
    "Mediana Posterior:   ", rho_median, "\n",
    "Desviación Estándar: ", rho_sd, "\n",
    "IC 95% Inferior:     ", rho_ic_inf, "\n",
    "IC 95% Superior:     ", rho_ic_sup, "\n")


varianza_u <- modelo_sar_hb$refVar                                              #Varianza de los efectos aleatorios 

betas <- modelo_sar_hb$coefficient[-8, "Mean"]                                  #se saca rho
betas <- as.numeric(betas)                                                      #Coeficientes estimados
rho_est <- modelo_sar_hb$coefficient["rho", "Mean"]                             #rho estimado

X_mat <- model.matrix(formula_sae, data = datos_modelo)                         #matriz de diseño
pred_fija <- X_mat %*% betas   
theta_hb <- modelo_sar_hb$Est[, "MEAN"]   
efectos_aleatorios <- theta_hb - pred_fija
residuos_directos <- y_s - theta_hb 


# Histograma
hist(efectos_aleatorios, main = "Histograma de Efectos Aleatorios", 
     xlab = "u_i", col = "lightblue", border = "darkgray")

# QQ-Plot
qqnorm(efectos_aleatorios, main = "QQ-Plot de Efectos Aleatorios")
qqline(efectos_aleatorios, col = "red", lwd = 2)


# Test estadístico formal de Shapiro-Wilk
shapiro.test(efectos_aleatorios)


# Gráfico de Residuos frente a Valores Predichos
plot(theta_hb, residuos_directos, 
     main = "Residuos Directos vs. Valores Predichos",
     xlab = "Predicción Fija", 
     ylab = "Residuos Brutos", 
     pch = 19, col = "darkblue")
abline(h = 0, col = "red", lwd = 2, lty = 2)



### Cálculo del factor de encogimiento ###
# cálculo de los pesos asegurando que sumen 1 exactamente
# G = A * [(I - rho*W_est)' * (I - rho*W_est)]^-1

C <- I_mat - rho_est * W_est
G <- A_est * solve(t(C) %*% C)                                                  #Matriz de covarianza espacial 


# Matriz de varianza de muestreo (Psi) y varianza total (Sigma)
psi <- diag(datos_modelo$vardir)                                                #Varianzas directas conocidas 
Sigma <- G + psi                                                                #Varianza total del modelo


# construcción de la Matriz de Pesos EBLUP 
# K = G * Sigma^-1
K <- G %*% solve(Sigma)



SAR_pesos_bay <- data.frame(
  n_i = strata_obs$n_h,                                                         #Tamaño de muestra por área
  psi = datos_modelo$vardir,                                                    #Varianza directa de muestreo
  sig2_v = varianza_u                                                           #Varianza espacial calculada localmente
) %>%
  mutate(
    
    Denominador = psi + sig2_v,                                                 #variabilidad local total del área
    Peso_Propio = sig2_v / Denominador,
    Peso_Vecindario = (psi / Denominador) * abs(rho_est),
    Peso_Sintetico = 1 - Peso_Propio - Peso_Vecindario
  ) %>%
  # Corrección marginal de precisión decimal por si algún valor da levemente menor a 0
  mutate(
    Peso_Sintetico = ifelse(Peso_Sintetico < 0, 0, Peso_Sintetico),
    stratum_id = row_number()
  )

# incorporar la conectividad geográfica con W estandarizada
SAR_pesos_bay$num_vecinos <- rowSums(W)

# Clasificación de conectividad, los facets según la cantidad de vecinos
SAR_pesos_bay <- SAR_pesos_bay %>%
  mutate(
    grupo_facet = case_when(
      num_vecinos %in% c(0, 1) ~ "0 o 1 Vecino (Baja conectividad)",
      num_vecinos == 2        ~ "2 Vecinos (Alta conectividad)"
    ),
    vecinos_cat = factor(num_vecinos, levels = c(0, 1, 2),
                         labels = c("0 vecinos", "1 vecino", "2 vecinos"))
  )


# vecindario vs global
p1 <- ggplot(SAR_pesos_bay, aes(x = Peso_Vecindario, y = Peso_Sintetico, color = n_i, alpha = n_i)) +
  geom_point(size = 4) +
  facet_wrap(~ grupo_facet) +
  scale_color_gradientn(colors = weights_cols_hb, limits = c(min(SAR_pesos_bay$n_i), max(SAR_pesos_bay$n_i))) +
  scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") + 
  labs(x = "Peso vecindario", y = "Peso de la Media Global", color = expression(n[i])) +
  theme_minimal() +
  theme(legend.position = "right")

# propio vs global
p2 <- ggplot(SAR_pesos_bay, aes(x = Peso_Propio, y = Peso_Sintetico, color = n_i, alpha = n_i)) +
  geom_point(size = 4) +
  facet_wrap(~ grupo_facet) +
  scale_color_gradientn(colors = weights_cols_hb, limits = c(min(SAR_pesos_bay$n_i), max(SAR_pesos_bay$n_i))) +
  scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") + 
  labs(x = "Peso propio", y = "Peso de la Media Global", color = expression(n[i])) +
  theme_minimal() +
  theme(legend.position = "right")


(p1 / p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Diagnóstico de Factores de Contracción del Modelo Espacial (Enfoque Bayesiano)",
    subtitle = "Comparativa según conectividad y tamaño de muestra"
  )



### Coeficientes de Variación ###
# Se extraen los errores cuadráticos medios (MSE) y medidas relativas de precisión 
#para cada estrato - CV calculado sobre la escala logarítmica del modelo

estimaciones_hb <- modelo_sar_hb$Est$MEAN
sd_posterior    <- modelo_sar_hb$Est$SD

# cálculo del Coeficiente de Variación (%)
cv_hb <- (sd_posterior / estimaciones_hb) * 100

# Unir resultados para los estratos
resultados_SARhb <- data.frame(
  strata_obs,
  mu_EBP = estimaciones_hb,
  mu_PSE = sd_posterior
)

resultados_SARhb <- resultados_SARhb %>%
  dplyr::select(stratum_id, strata, n_h, y_h, var_h, mu_EBP, mu_PSE) 


resultados_SARhb$mu_RSE_log <- (resultados_SARhb$mu_PSE / abs(resultados_SARhb$mu_EBP)) * 100

resultados_SARhb <- resultados_SARhb %>% 
  arrange(n_h) %>%              
  mutate(order_n = row_number())


ggplot(resultados_SARhb, aes(x = order_n, y = mu_RSE_log)) +
  geom_point() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Error Relativo Estándar por Estrato", 
       x = "Estratos (ordenados por tamaño)", y = "RSE (%)") +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )


##### Visualización de medias #####
# df para predicciones
df_predicciones_SARhb <- strata_obs

df_predicciones_SARhb$edad_real <- exp(estimaciones_hb)

df_predicciones_SARhb$ethnic_rec <- fct_relevel(
  df_predicciones_SARhb$ethnic_rec,
  "otros", "indigena", "blanco/mestizo", "afro"
)

df_predicciones_SARhb$type_area <- fct_relevel(
  df_predicciones_SARhb$type_area, "urbana", "rural"
)

df_predicciones_SARhb$high_edu <- fct_relevel(
  df_predicciones_SARhb$high_edu, 
  "superior", "secundaria",  "primaria", "sin educación"
)


### Mapa de calor - combinaciones de variables ###
ggplot(df_predicciones_SARhb, aes(x = cohort_quinquenal, y = high_edu, fill = edad_real)) +
  
  geom_tile(color = "white") +
  
  # divide matriz por área (filas) y Etnia (columnas)
  facet_grid(type_area ~ ethnic_rec) + 
  
  scale_fill_gradientn(
    colours = grunge_fun(100),
    name = "Edad \nPrimer \nHijo",
    limits = c(min(df_predicciones_SARhb$edad_real),
               max(df_predicciones_SARhb$edad_real)),
    oob = scales::squish
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Cohorte de Nacimiento",
    y = "Nivel Educativo") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )



### diagrama de dispersión de la edad promedio predicha por estrato - ordenado 
#por tamaño ###

df_plot_SARhb <- df_predicciones_SARhb %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n_h) %>%
  mutate(order_n = row_number())


ggplot(df_plot_SARhb) +
  geom_point(aes(x = order_n, y = edad_real),
             color = "black",
             shape = "✦",
             alpha = 0.6,      
             size = 5) +     
  labs(x = "Estratos (ordenados por tamaño)",
       y = "Mediana de Edad Estimada") +
  theme(
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)  
  )


### Exportación de resultados ###
# Se almacenan los efectos aleatorios estimados y las medidas de precisión para su 
#uso posterior en comparación de modelos.
write.csv(resultados_SARhb, "resultados_estratos_SARhb.csv", row.names = FALSE)
