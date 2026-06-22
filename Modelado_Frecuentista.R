################################################################################
# Modelado_Frecuentista.R
#
# Objetivo:
# Ajustar y comparar modelos de estimación para áreas pequeñas
# bajo un enfoque frecuentista utilizando los estratos construidos
# a partir de variables sociodemográficas.
#
# El script realiza:
#   1. Construcción de estimadores directos y varianzas de muestreo.
#   2. Ajuste del modelo Fay-Herriot.
#   3. Ajuste del modelo Battese-Harter-Fuller (Nested Error Regression).
#   4. Diagnóstico y evaluación de supuestos.
#   5. Obtención de predicciones SAE por estrato.
#   6. Visualización y almacenamiento de resultados.
#   7. Exportación de resultados para análisis posteriores.
#
# Autor: Selene Alvarado Rosario
# Fecha: Diciembre 2025
################################################################################


# Ajustar ruta de trabajo según la ubicación del proyecto


library(dplyr)
library(ggplot2)
library(sae)
library(hbsae)


df_areas <- read.csv("complete_birth_rec.csv", stringsAsFactors = FALSE)        #Importo el dataset que se había preparado, ya con ethnic_rec

df_areas$type_area <- factor(
  df_areas$type_area,
  levels = c("urbana", "rural")
)
df_areas$high_edu <- factor(
  df_areas$high_edu,
  levels = c("superior", "secundaria", "primaria", "sin educación")
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

stratum_id <- df_indiv$stratum_id


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

##### Ajuste y calibración de pesos muestrales #####
y_i <- log(df_indiv$age_1st_b)                                                  #variable respuesta en escala log

v <- df_indiv$sam_weight                                                        #extracción de pesos

par(mfrow = c(1,1))
boxplot(v)

# cuartiles
Q1 <- quantile(v, 0.25)
Q3 <- quantile(v, 0.75)

# Definir el umbral v0
v0 <- Q3 + 1.5 * (Q3 - Q1)

# Identificar outliers
outlier <- v >= v0

N_hat <- sum(v)                                                                 #Suma de pesos muestrales (para estimador HT)
r <- (N_hat - sum(outlier) * v0) / sum(v[!outlier])  

v_estrellita <- ifelse(outlier, v0, r * v)

boxplot(v_estrellita)

# Verificar suma (debe ser igual a N_hat)
sum(v_estrellita) == N_hat                                                      #Se mantiene

df_indiv$sam_weight_adj <- v_estrellita

N_hat_by_stratum <- df_indiv %>%
  group_by(stratum_id) %>%
  summarise(
    N_i_hat = sum(sam_weight_adj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stratum_id)

N_vec <- N_hat_by_stratum$N_i_hat


##### Modelo Fay-Herriot #####

### Varianza global de respaldo ###
# Se utiliza para imputar la varianza de muestreo en estratos con n_i = 1, donde la varianza 
#muestral no puede calcularse.
S_global <- as.numeric(var(log(df_indiv$age_1st_b), na.rm = TRUE))


### Construcción de las áreas pequeñas (estratos observados) ###
# Se calcula el estimador directo de la media logarítmica, la varianza de muestreo aproximada 
#y el tamaño muestral.
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


### Ajuste del modelo Fay-Herriot clásico ###
# utilizando las estimaciones directas y sus varianzas de muestreo.
fh_age_log <- mseFH(y_dir ~ cohort_quinquenal + high_edu + type_area + ethnic_rec, 
                    psi_i, data = df_small_area)

str(fh_age_log)

plot(fh_age_log$est$eblup)


### Diagnóstico y verificación de supuestos ###
# Diagnóstico gráfico de residuos y efectos aleatorios.

# residuos de los modelos
ajustados <- fh_age_log$est$eblup
residuos <- df_small_area$y_dir - ajustados 


# Histogramas
hist(residuos, main = "Distribución de Residuos", xlab = "Residuos de Pearson", col = "lightblue")


# Gráficos Q-Q norm
qqnorm(residuos)
qqline(residuos, col = "red")


# Residuos vs. Valores Ajustados
plot(ajustados, residuos, main = "Residuos vs. Ajustados", 
     xlab = "Valores Predichos", ylab = "Residuos")
abline(h = 0, col = "red")

plot(df_small_area$y_dir, ajustados,
     pch = 19,
     xlab = "Estimador Directo",
     ylab = "EBLUP Fay-Herriot")

abline(0,1,col="red",lwd=2)

fh_age_log$est$fit$goodness


# Efectos aleatorios - Matriz de diseño
# Se construye la matriz de covariables utilizada para separar la parte fija y
# aleatoria del modelo.
X <- model.matrix(~cohort_quinquenal + high_edu + type_area + ethnic_rec, df_small_area)


efectos_u <- fh_age_log$est$eblup - X %*% fh_age_log$est$fit$estcoef$beta
hist(efectos_u, main="Normalidad de Efectos Aleatorios", xlab="u_d")
qqnorm(efectos_u); qqline(efectos_u, col="blue")


### Cálculo del factor de encogimiento ###
# valores cercanos a 1 indican mayor dependencia del modelo,
# mientras que valores cercanos a 0 indican mayor peso del estimador directo.

residuos_lm <- df_small_area$y_dir - X %*% fh_age_log$est$fit$estcoef$beta
refvar <- fh_age_log$est$fit$refvar

df_small_area$gamma <- refvar / (refvar + df_small_area$psi_i)

plot(residuos_lm)
plot(df_small_area$gamma)

summary(df_small_area$gamma)

ggplot(df_small_area, aes(x = reorder(n_i, gamma), y = gamma)) +
  geom_point() +
  labs(title = "Coeficiente de Shrinkage por Estrato - FH", 
       x = "Estratos (ordenados por tamaño)", y = expression(gamma[i])) +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )


### Coeficientes de Variación ###
# Se extraen los errores cuadráticos medios (MSE) y medidas relativas de precisión 
#para cada estrato - CV calculado sobre la escala logarítmica del modelo
resultados_estratos_fh <- data.frame(
  strata = df_small_area$strata,
  n_i = df_small_area$n_i,
  u_h = efectos_u,
  psi = df_small_area$psi_i
) %>%
  mutate(
    se_u_h = sqrt((refvar * psi) / (refvar + psi))
  ) %>%
  dplyr::select(strata, n_i, u_h, se_u_h)


resultados_estratos_fh <- resultados_estratos_fh %>%
  left_join(strats, by = "strata") %>%
  mutate(
    u_h_real = exp(u_h),
    mse = fh_age_log$mse,
    
    edad_hat_log = fh_age_log$est$eblup,
    edad_hat = exp(edad_hat_log),
    
    RRMSE = 100 * sqrt(mse) / edad_hat,
    
    RSE_log = 100 * sqrt(mse) / abs(edad_hat_log)
  )

resultados_estratos_fh <- resultados_estratos_fh %>%
  arrange(n_i) %>%  
  mutate(order_n = row_number())


ggplot(resultados_estratos_fh, aes(x = order_n, y = RSE_log)) +
  geom_point() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Error Relativo Estándar por Estrato", 
       x = "Estratos (ordenados por tamaño)", y = "RSE (%)") +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )


##### Visualización de medias #####
# df para predicciones
fh_predicciones <- df_indiv %>%
  group_by(stratum_id, strata, cohort_quinquenal, high_edu, type_area, ethnic_rec) %>%
  summarise(n = n(), .groups = "drop")

# predicción en escala log
fh_predicciones$pred_log <- as.numeric(fh_age_log$est$eblup)

# volver a la escala real (años de edad)
fh_predicciones <- fh_predicciones %>%
  mutate(edad_real = exp(pred_log))


### Mapa de calor - combinaciones de variables ###
ggplot(fh_predicciones, aes(x = cohort_quinquenal, y = high_edu, fill = edad_real)) +
  
  geom_tile(color = "white") +
  
  # divide matriz por área (filas) y Etnia (columnas)
  facet_grid(type_area ~ ethnic_rec) + 
  
  scale_fill_gradientn(
    colours = grunge_fun(100),
    name = "Edad \nprimer \nHijo",
    limits = c(min(fh_predicciones$edad_real),
               max(fh_predicciones$edad_real)),
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
fh_plot <- fh_predicciones %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n) %>%
  mutate(order_n = row_number())


ggplot(fh_plot) +
  geom_point(aes(x = order_n, y = edad_real),
             color = "black",
             shape = "✦",   #si falla, cambiar por 18
             alpha = 0.6,      
             size = 5) +    
  labs(x = "Estratos (ordenados por tamaño)",
       y = "Media de Edad Estimada") +
  theme(
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)  
  )


### Exportación de resultados ###
# Se almacenan los efectos aleatorios estimados y las medidas de precisión para su 
#uso posterior en comparación de modelos.
write.csv(resultados_estratos_fh, "resultados_estratos_fh.csv", row.names = FALSE)


# ------------------------------------------------------------------------------- #
#####  Modelo de Regresión de Errores Anidados #####

### Ajuste del modelo NER ###
X_unit_df <- as.data.frame(model.matrix(~ type_area + high_edu + ethnic_rec + cohort_quinquenal, data = df_indiv))


# ID del estrato y los pesos 
X_unit_df$stratum_id <- df_indiv$stratum_id
X_unit_df$sam_weight_adj <- df_indiv$sam_weight_adj


# medias poblacionales (proporciones) por estrato
Xpop_mat <- X_unit_df %>%
  group_by(stratum_id) %>%
  summarise(across(
    .cols = -sam_weight_adj, 
    .fns = ~ weighted.mean(., w = sam_weight_adj)                               #Media ponderada, consistencia con N_vec
  )) %>%
  arrange(stratum_id) %>%
  as.matrix()
  

all(
  as.character(N_hat_by_stratum$stratum_id) ==
    rownames(Xpop_mat)
)

row.names(Xpop_mat) <- Xpop_mat[,1]                                             #columna de stratum id
Xpop_mat <- Xpop_mat[,-1]


X_ner <- X_unit_df %>%
  dplyr::select(-stratum_id, -sam_weight_adj) %>%
  as.matrix()


modelo_ner <- fSAE.Unit(
  y = y_i,
  X = X_ner,  
  area = stratum_id,                                                            #Variable que define los dominios
  Narea = N_vec,
  Xpop = Xpop_mat,
  method = "REML"                                                               #por Máxima verosimilitud 
)

print(modelo_ner, correlation = TRUE)


### Diagnóstico y verificación de supuestos ###
# Diagnóstico gráfico de residuos y efectos aleatorios.

# residuos del modelo
res_log <- residuals(modelo_ner, type = "pearson")
adj_log <- fitted(modelo_ner)


# Histogramas
hist(res_log, main = "Distribución de Residuos (log)", xlab = "Residuos de Pearson", col = "lightblue")


# Gráficos Q-Q norm
qqnorm(res_log)
qqline(res_log, col = "red")


# Residuos vs. Valores Ajustados
plot(adj_log, res_log, main = "Residuos vs. Ajustados (log)", 
     xlab = "Valores Predichos", ylab = "Residuos")
abline(h = 0, col = "red")


#Criterios de información
logL <- modelo_ner$llh.c
k <- length(modelo_ner$beta) + 2                                                #sigma_e y sigma_u
logn <- log(length(modelo_ner$y))                                               #tamaño de muestra total

AIC <- 2*k - 2*logL
BIC <- -2 * logL + k*logn
KIC <- -2 * logL + 3 * k
AIC; BIC; KIC


# Extracción de efectos aleatorios
efectos_u <- raneff(modelo_ner)
hist(efectos_u, main="Normalidad de Efectos Aleatorios", xlab="u_d")
qqnorm(efectos_u); qqline(efectos_u, col="blue")

### Cálculo del factor de encogimiento ###
# valores cercanos a 1 indican mayor dependencia del modelo,
# mientras que valores cercanos a 0 indican mayor peso del estimador directo.

# Extraer las varianzas del modelo
sig_u2_ner <- sv2(modelo_ner)                                                   #Varianza entre áreas
sig_e2_ner <- se2(modelo_ner)                                                   #Varianza residual
ratio <- sig_u2_ner / sig_e2_ner

rho <- ratio / (1+ratio)
rho


# Obtener el factor gamma por estrato
shrinkage_df_ner <- df_indiv %>%
  group_by(stratum_id) %>%
  summarise(n_i = n()) %>%
  mutate(
    gamma_i = modelo_ner$gamma
  )

# Ver los primeros resultados
head(shrinkage_df_ner)

summary(shrinkage_df_ner$gamma_i)

shrinkage_df_ner <- shrinkage_df_ner %>%
  arrange(n_i) %>%
  mutate(order_n = row_number()) %>%
  mutate(stratum_id = as.character(stratum_id))


ggplot(shrinkage_df_ner, aes(x = order_n, y = gamma_i)) +
  geom_point() +
  labs(title = "Coeficiente de Shrinkage por Estrato - NER", 
       x = "Estratos (ordenados por tamaño)", y = expression(gamma[i])) +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )


### Precisión de las estimaciones SAE ###
# Se utilizan los errores relativos obtenidos mediante relSE(), equivalentes 
#al coeficiente de variación (CV) de las predicciones por estrato.

# Se extraen los coeficientes de variación (CV) y medidas relativas de 
#precisión para cada estrato - CV calculado sobre la escala logarítmica del modelo
cv_reales <- relSE(modelo_ner, type = "sae")  


#Crear el dataframe de resultados por estrato
resultados_estratos_ner <- data.frame(
  stratum_id = as.character(modelo_ner$sampledAreaNames),
  u_h = as.numeric(efectos_u),
  se_u_h = sqrt(modelo_ner$Vraneff),
  CV = as.numeric(cv_reales) *100
)

resultados_estratos_ner <- resultados_estratos_ner %>%
  left_join(shrinkage_df_ner, by = "stratum_id") %>%
  left_join(strats, by = "stratum_id") %>%
  arrange(n_i) %>%
  mutate(order_n = row_number())

ggplot(resultados_estratos_ner, aes(x = order_n, y = CV)) +
  geom_point() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Coeficientes de Variación por Estrato", 
       x = "Estratos (ordenados por tamaño)", y = "CV (%)") +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )



##### Visualización de medias #####
# df para predicciones
df_predicciones_ner <- df_indiv %>%
  mutate(high_edu = factor(high_edu, levels = c("superior", "secundaria", "primaria", "sin educación")),
         ethnic_rec = factor(ethnic_rec, levels = c("otros", "indigena", "blanco/mestizo", "afro")),
         type_area = factor(type_area, levels = c("urbana", "rural"))) %>%
  group_by(stratum_id, strata, cohort_quinquenal, high_edu, type_area, ethnic_rec) %>%
  summarise(n = n(), .groups = "drop")

# predicción en escala log
df_predicciones_ner$pred_log <- modelo_ner$est


# Corrección de sesgo por retransfomación lognormal.
# Dado que el modelo se ajusta sobre log(age_1st_b), las predicciones se devuelven 
# a la escala original mediante E(Y)=exp(eta_hat)
# donde eta_hat corresponde al predictor EBLUP obtenido en escala logarítmica.

df_predicciones_ner <- df_predicciones_ner %>%
  mutate(edad_real = exp(pred_log))


### Mapa de calor - combinaciones de variables ###
ggplot(df_predicciones_ner, aes(x = cohort_quinquenal, y = high_edu, fill = edad_real)) +
  
  geom_tile(color = "white") +
  
  #divide matriz por área (filas) y Etnia (columnas)
  facet_grid(type_area ~ ethnic_rec) + 
  
  scale_fill_gradientn(
    colours = grunge_fun(100),
    name = "Edad \nprimer \nHijo",
    limits = c(min(df_predicciones_ner$edad_real),
               max(df_predicciones_ner$edad_real)),
    oob = scales::squish
  ) +
  
  theme_minimal() +
  
  labs(x = "Cohorte de Nacimiento",
       y = "Nivel Educativo") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )



### diagrama de dispersión de la edad promedio predicha por estrato - ordenado por tamaño ###
ner_plot <- df_predicciones_ner %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n) %>%
  mutate(order_n = row_number())


ggplot(ner_plot) +
  geom_point(aes(x = order_n, y = edad_real),
             color = "black",
             shape = "✦",
             alpha = 0.6,      
             size = 5) +     
  labs(x = "Estratos (ordenados por tamaño)",
       y = "Media de Edad Estimada") +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)  
  )

### Exportación de resultados ###
# Se almacenan los efectos aleatorios estimados y las medidas de precisión para su 
#uso posterior en comparación de modelos.
write.csv(resultados_estratos_ner, "resultados_estratos_ner.csv", row.names = FALSE)
