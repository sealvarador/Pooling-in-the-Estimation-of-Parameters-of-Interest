################################################################################
# Modelado_Bayesiano.R
#
# Objetivo:
# Ajustar y comparar modelos de estimación para áreas pequeñas
# bajo un enfoque bayesiano utilizando los estratos construidos 
# a partir de variables sociodemográficas.
#
# El script realiza:
#   1. Construcción de estimadores directos y varianzas de muestreo.
#   2. Ajuste del modelo Scott-Smith.
#   3. Ajuste del modelo Battese-Harter-Fuller bayesiano jerárquico.
#   4. Diagnóstico y evaluación de supuestos.
#   5. Obtención de predicciones SAE por estrato.
#   6. Visualización y almacenamiento de resultados.
#   7. Exportación de resultados para análisis posteriores.
#
# Autor: Selene Alvarado Rosario
# Fecha: Febrero 2025
################################################################################


# Ajustar ruta de trabajo según la ubicación del proyecto

library(MVN)
library(dplyr)
library(coda)
library(hbsae)
library(ggplot2)
library(scales)


df_areas <- read.csv("complete_birth_rec.csv", stringsAsFactors = FALSE)        #Importo el dataset que se había preparado, ya con ethnic_rec


df_areas <- df_areas %>%
  mutate(y_s = log(age_1st_b))

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
    y_h = mean(y_s, na.rm = TRUE),
    var_h = var(y_s, na.rm = TRUE) / n_h,
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

#parámetros y estadísticos necesarios
y_s <- strata_obs$y_h                                                           #solo incluye la variable observada
n_i <- strata_obs$n_h
n <- sum(n_i)
strats_ids  <- strata_obs$stratum_id
l <- length(strats_ids)



ssw_est <- df_areas %>%
  group_by(type_area, high_edu, ethnic_rec, cohort_quinquenal) %>%
  summarise(sum_sq_diff = sum((y_s - mean(y_s, na.rm = TRUE))^2, na.rm = TRUE), .groups = "drop") %>%
  pull(sum_sq_diff) %>%
  sum()

#####Modelo Scott-Smith#####
set.seed(42)                                                                    #semilla para reproducibilidad

y_sim <- function(rho, y_s, n_i){
  denom <- ((n_i - 1) * rho) + 1
  full <- n_i / denom
  y_sim <- sum(full * y_s) / sum(full)
  return(y_sim)
}

lambda_i <- function(n_i, rho){
  num <- n_i * rho
  denom <- (n_i - 1) * rho + 1
  return(num/denom)
}


#grilla
grid_length <- 100
rho_grid <- seq(1e-3, 1-1e-3, length.out = grid_length)


logpost_rho <- function(rho, n_i, y_s, ssw_est = 1e-3, l, n){
  ysim <- y_sim(rho, y_s, n_i)
  denom <- ((n_i - 1) * rho) + 1
  full <- n_i / denom
  
  p1 <- ((l-2)/2) * log(1-rho)
  
  p2 <- sum(log(full)) - log(sum(full))
  p2 <- 0.5 * p2
  
  p3_num <- full * (y_s - ysim) * (y_s - ysim)
  p3_num <- sum(p3_num)
  
  p3 <- p3_num / ssw_est
  p3 <- 1 + (1-rho) * p3
  p3 <- -((n-1)/2) * log(p3)
  
  return(p1 + p2 + p3)
}


### Muestreador ###
FS <- 1900                                                                      #Tamaño final 1900
I <- 20000                                                                      #Num. de iteraciones total 20k
C <- 1000                                                                       #Periodo de calentamiento 1k


scott_smith_gibbs <- function(y_s, n_i, l, n_total, ssw_est, rho_grid, FS, I, C, stratum_ids, thin = 10){
  # almacenamiento
  chain   <- matrix(NA, nrow = FS, ncol = 3)
  mu_samp <- matrix(NA, nrow = FS, ncol = length(n_i))
  LL      <- matrix(data = NA, nrow = FS, ncol = 1)
  Lambda   <- matrix(data = NA, nrow = FS, ncol = length(n_i))
  
  # inicialización
  rho_ss  <- rho_grid[2]
  
  fs <- 1
  set.seed(151724)
  for (t in 1:I) {
    # Muestreo de rho
    logdens_rho <- sapply(rho_grid, function(r)
      logpost_rho(r, n_i, y_s, ssw_est, l, n_total))
    dens_rho <- exp(logdens_rho - max(logdens_rho))
    rho_ss <- sample(rho_grid, size = 1, prob = dens_rho)
    
    # Parámetros intermedios
    alpha <- (1 - rho_ss) / rho_ss
    lambdai <- lambda_i(n_i, rho_ss)
    ysim <- y_sim(rho_ss, y_s, n_i)
    
    # Muestreo de sigma2
    suma_cuadrados_modelo <- alpha * sum(lambdai * (y_s - ysim)^2)
    shape_sig <- (n_total - 1) / 2 
    scale_sig <- (ssw_est + suma_cuadrados_modelo) / 2
    
    sigma2 <- 1 / rgamma(1, shape = shape_sig, rate = scale_sig)
    
    # Muestreo de theta
    var_theta <- (sigma2 * rho_ss) / ((1 - rho_ss) * sum(lambdai))
    theta <- rnorm(1, mean = ysim, sd = sqrt(var_theta))
    
    # Muestreo de mu_h
    res <- 1 - lambdai
    mean_mu <- (lambdai * y_s) + (res * theta)
    var_mu  <- (res * rho_ss * sigma2) / (1 - rho_ss)
    mu <- rnorm(l, mean = mean_mu, sd = sqrt(var_mu))
    
    # Almacenar muestras
    if (t > C && t %% thin == 0) {
      if (fs <= FS) {
        chain[fs, ]   <- c(rho_ss, theta, sigma2)
        mu_samp[fs, ] <- mu
        
        # Log-verosimilitud predictiva
        ll_h <- dnorm(x = y_s, mean = mu, sd = sqrt(sigma2/n_i), log = TRUE)
        LL[fs, ] <- sum(ll_h)
        Lambda[fs, ] <- lambdai
        fs <- fs + 1
      }
    }
    if (fs > FS) break 
  }
  
  # ASIGNACIÓN DE NOMBRES PARA EL RASTREO
  colnames(chain)   <- c("rho", "theta", "sigma2")
  colnames(mu_samp) <- stratum_ids
  colnames(Lambda) <- stratum_ids
  
  return(list(chain = chain, mu_samp = mu_samp, LL = LL, Lambdai = Lambda))
}

scott_smith_mcmc <- scott_smith_gibbs(
  y_s = y_s, 
  n_i = n_i, 
  l = l, 
  n_total = n, 
  ssw_est = ssw_est, 
  rho_grid = rho_grid, 
  FS = FS, I = I, C = C,
  stratum_ids = strats_ids
)



### Diagnóstico y verificación de supuestos ###
# Convierte la cadena a objeto mcmc
colnames(scott_smith_mcmc$mu_samp) <- paste0("mu", 1:l)
colnames(scott_smith_mcmc$chain) <- c("rho", "theta", "sigma2")
mcmc_ch_ss <- coda::mcmc(scott_smith_mcmc$chain)
mcmc_mu_ss <- coda::mcmc(scott_smith_mcmc$mu_samp)

ev_ss <- cbind(scott_smith_mcmc$chain, scott_smith_mcmc$mu_samp)

# Traceplots de los hiperparámetros y parámetros escalares
plot(mcmc_ch_ss[, "rho"], main="Traceplot rho")
plot(mcmc_ch_ss[, "theta"], main="Traceplot theta")
plot(mcmc_ch_ss[, "sigma2"], main="Traceplot sigma²")
plot(scott_smith_mcmc$LL, col = "#7CCD7C", main = "Log-verosimilitud de la predicción", ylab = "LL",
     pch = "⭐", cex =  1.2)


# Densidades posteriores
par(mfrow=c(2,2))
hist(scott_smith_mcmc$chain[, "rho"], main="Posterior rho", freq=FALSE, col = "#76EEC6", 
     xlab = "rho")
hist(scott_smith_mcmc$chain[, "theta"], main="Posterior theta", freq=FALSE, col = "#FF82AB",
     xlab = "theta")
hist(scott_smith_mcmc$chain[, "sigma2"], main="Posterior sigma²", freq=FALSE, col = "#D264A5",
     xlab = "sigma²")

# Autocorrelación
par(mfrow=c(1,1))
acf(scott_smith_mcmc$chain[, "rho"], main="ACF rho")
acf(scott_smith_mcmc$chain[, "theta"], main="ACF theta")
acf(scott_smith_mcmc$chain[, "sigma2"], main="ACF sigma²")
acf(scott_smith_mcmc$mu_samp[, c(10)])                                          #se puede verificar cualquiera entre 1-l

#rho
mean_rho_ss <- mean(scott_smith_mcmc$chain[,"rho"])
lower_ci_rho_ss <- quantile(scott_smith_mcmc$chain[,"rho"], 0.05)
median_rho_ss <- quantile(scott_smith_mcmc$chain[,"rho"], 0.5)
upper_ci_rho_ss <- quantile(scott_smith_mcmc$chain[,"rho"], 0.95)
print(c(mean_rho_ss, lower_ci_rho_ss, upper_ci_rho_ss))
sd(scott_smith_mcmc$chain[,"rho"])

#theta
theta_post_mean_ss <- mean(scott_smith_mcmc$chain[, "theta"])
theta_post_ci_ss   <- quantile(scott_smith_mcmc$chain[, "theta"], c(0.05, 0.95))
print(c(theta_post_mean_ss, theta_post_ci_ss))


#sigma
sigma2_muestras <- scott_smith_mcmc$chain[,"sigma2"]

min_sigma <- min(scott_smith_mcmc$chain[,"sigma2"])
mean_sigma <- mean(scott_smith_mcmc$chain[,"sigma2"])
median_sigma <- quantile(scott_smith_mcmc$chain[,"sigma2"], 0.5)
max_sigma <- max(scott_smith_mcmc$chain[,"sigma2"])
print(c(min_sigma, mean_sigma, median_sigma, max_sigma))


#tau
tau2_muestras <- (scott_smith_mcmc$chain[, "rho"] / (1 - scott_smith_mcmc$chain[, "rho"])) * 
  scott_smith_mcmc$chain[, "sigma2"]

mean_tau2   <- mean(tau2_muestras)
median_tau2 <- quantile(tau2_muestras, 0.5)
ic_tau2     <- quantile(tau2_muestras, c(0.05, 0.95))

print(c(mean_tau2, ic_tau2))

# Estratos seleccionados únicamente con fines ilustrativos
mardia(scott_smith_mcmc$mu_samp[, c(9,25)])
plot(scott_smith_mcmc$mu_samp[,9], scott_smith_mcmc$mu_samp[,25])

summary(scott_smith_mcmc$mu_samp)
acf(scott_smith_mcmc$mu_samp[,41], lag.max = 50)

sum(is.na(scott_smith_mcmc$chain))

neffSS <- coda::effectiveSize(ev_ss)
round(neffSS, 0)

EMC_SS <- apply(X = ev_ss, MARGIN = 2, FUN = sd)/sqrt(neffSS)
round(EMC_SS, 4)

coda::geweke.diag(mcmc_ch_ss)


### Población Finita Posterior ###
# Matriz para almacenar las muestras de la población finita
Y_pop_muestras <- matrix(NA, nrow = FS, ncol = l)
f_h <- n_i / N_vec


for (i in 1:l) {
  # Componentes de la distribución predictiva (Ecuación A.6) [2]
  # Media: f_i*y_bar_i + (1 - f_i)*mu_i
  media_pred <- f_h[i] * y_s[i] + (1 - f_h[i]) * mcmc_mu_ss[, i]
  
  # Varianza: (1 - f_i) * (sigma^2 / N_i)
  var_pred <- (1 - f_h[i]) * (sigma2_muestras / N_vec[i])
  
  # Generación de la muestra predictiva para cada iteración
  Y_pop_muestras[, i] <- rnorm(FS, mean = media_pred, sd = sqrt(var_pred))
}


# Estimaciones Puntuales (Media Posterior)
ebp_mu <- colMeans(scott_smith_mcmc$mu_samp[, 1:l])
ebp_theta <- mean(scott_smith_mcmc$chain[, "theta"])
ebp_sigma2 <- mean(scott_smith_mcmc$chain[, "sigma2"])
ebp_rho <- mean(scott_smith_mcmc$chain[, "rho"])
y_ebp_ss <- colMeans(Y_pop_muestras)



##### Análisis de resultados #####
par(mfrow = c(1,1))
hist(exp(y_ebp_ss), main="Medias posteriores Scott-Smith", freq=FALSE, col = "#76EEC6",
     xlab = "Medias posteriores")


# por covariable:
Yh_post_mean_ss   <- y_ebp_ss
mu_post_low_ss    <- apply(Y_pop_muestras, 2, quantile, 0.025)
mu_post_high_ss   <- apply(Y_pop_muestras, 2, quantile, 0.975)

df_mu_ss <- data.frame(
  stratum_id = strats_ids,
  mu_post = round(exp(Yh_post_mean_ss), 4)
)

strata_obs_t <- strata_obs %>% 
  left_join(df_mu_ss, by = "stratum_id") %>%
  left_join(N_hat_by_stratum, by = "stratum_id")

df_mu_cohort_ss <- strata_obs_t %>% 
  group_by(cohort_quinquenal) %>% 
  summarise(
    mu_cohort = round(weighted.mean(mu_post, w = N_i_hat),4),
    sd_cohort = round(sd(mu_post),4),
    cv_cohort = round(sd(mu_post)/weighted.mean(mu_post, w = N_i_hat),4),
    CrI_low_cohort = round(quantile(mu_post, probs = 0.05),4),
    CrI_up_cohort = round(quantile(mu_post, probs = 0.95),4),
    .groups = "drop"
  )

as.data.frame(df_mu_cohort_ss)

df_mu_etnia <- strata_obs_t %>% 
  group_by(ethnic_rec) %>% 
  summarise(
    mu_etnia = round(weighted.mean(mu_post, w = N_i_hat),4),
    sd_etnia = round(sd(mu_post),4),
    cv_etnia = round(sd(mu_post)/weighted.mean(mu_post, w = N_i_hat),4),
    CrI_low_eth = round(quantile(mu_post, probs = 0.05),4),
    CrI_up_eth = round(quantile(mu_post, probs = 0.95),4),
    .groups = "drop"
  )

as.data.frame(df_mu_etnia)

df_mu_area <- strata_obs_t %>% 
  group_by(type_area) %>% 
  summarise(
    mu_area = round(weighted.mean(mu_post, w = N_i_hat),4),
    sd_area = round(sd(mu_post),4),
    cv_area = round(sd(mu_post)/weighted.mean(mu_post, w = N_i_hat),4),
    CrI_low_area = round(quantile(mu_post, probs = 0.05),4),
    CrI_up_area = round(quantile(mu_post, probs = 0.95),4),
    .groups = "drop"
  )

as.data.frame(df_mu_area)

df_mu_edu <- strata_obs_t %>% 
  group_by(high_edu) %>% 
  summarise(
    mu_edu = round(weighted.mean(mu_post, w = N_i_hat),4),
    sd_edu = round(sd(mu_post),4),
    cv_edu = round(sd(mu_post)/weighted.mean(mu_post, w = N_i_hat),4),
    CrI_low_edu = round(quantile(mu_post, probs = 0.05),4),
    CrI_up_edu = round(quantile(mu_post, probs = 0.95),4),
    .groups = "drop"
  )

as.data.frame(df_mu_edu)


# Precisión (Error Estándar Posterior, PSE)
pse_mu <- apply(scott_smith_mcmc$mu_samp[, 1:l], 2, sd)
pse_rho <- sd(scott_smith_mcmc$chain[, "rho"])
pse_ss <- apply(Y_pop_muestras, 2, sd)

# Intervalos de Credibilidad del 95%
ic_mu <- apply(scott_smith_mcmc$mu_samp[, 1:l], 2, quantile, probs = c(0.025, 0.975))

# Unir resultados para los estratos
resultados_ss <- data.frame(
  strata_obs,
  mu_post_mean = ebp_mu,
  mu_PSE = pse_mu,
  IC_025 = ic_mu[1, ],
  IC_975 = ic_mu[2, ]
)

resultados_ss <- resultados_ss %>%
  dplyr::select(stratum_id, strata, n_h, y_h, var_h, mu_post_mean, mu_PSE, 
                IC_025, IC_975)

resultados_ss$mu_CV <- (resultados_ss$mu_PSE / resultados_ss$mu_post_mean) * 100
resultados_ss$mu_post_CV <- (pse_ss / y_ebp_ss) * 100                           #CV calculado sobre la escala logarítmica del modelo

resultados_ss <- resultados_ss %>%
  arrange(n_h) %>%              
  mutate(order_n = row_number())



ggplot(resultados_ss, aes(x = order_n, y = mu_post_CV)) +
  geom_point() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Coeficientes de Variación por Estrato", 
       x = "Estratos (ordenados por tamaño)", y = "CV (%)") +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )

# Resumen de coeficientes de variación
CVMC_SS_sum <- as.table(c(min(resultados_ss$mu_CV), mean(resultados_ss$mu_CV),
                          max(resultados_ss$mu_CV)))

CVMC_SS_sum <- round(CVMC_SS_sum, 4)
names(CVMC_SS_sum) <- c("Mínimo mu_j", "Media mu_j", "Máximo mu_j")
CVMC_SS_sum


### Criterios de información ###
LP_SS        <- as.numeric(scott_smith_mcmc$LL)
theta_hat_ss <- mean(scott_smith_mcmc$chain[,c("theta")])
sigma2_hat_ss <- mean(scott_smith_mcmc$chain[,c("sigma2")])
mu_hat_ss <- colMeans(scott_smith_mcmc$mu_samp)
lpyth_ss   <- sum(dnorm(x = strata_obs$y_h, mean = mu_hat_ss, 
                        sd = sqrt(sigma2_hat_ss/n_i), log = T))                 #verosimilitud evaluada en la media posterior
pDIC_ss    <- 2*(lpyth_ss - mean(LP_SS))
dic_ss     <- -2*lpyth_ss + 2*pDIC_ss
print(c(lpyth_ss,dic_ss))


# Resumen Global
cat("\n--- Resumen del Modelo Scott-Smith ---\n")
cat("Estimación EBP de la Media Global (theta):", round(ebp_theta, 3), "\n")
cat("Estimación EBP de la Varianza Individual (sigma^2):", round(ebp_sigma2, 3), "\n")
cat("Estimación EBP de la Correlación Intraclase (rho):", round(ebp_rho, 5), " (PSE:", round(pse_rho, 5), ")\n")
cat("\n--- Resultados por Estrato (EBP de mu_h) ---\n")
print(head(resultados_ss))



##### visualizar medias #####

# df para predicciones
df_predicciones_ss <- strata_obs

# predicción
df_predicciones_ss$edad_esc_real <- exp(y_ebp_ss)  #edad en escala real/original


# re-organización de niveles para visualización
df_predicciones_ss$ethnic_rec <- fct_relevel(
  df_predicciones_ss$ethnic_rec,
  "otros", "indigena", "blanco/mestizo", "afro"
)

df_predicciones_ss$type_area <- fct_relevel(
  df_predicciones_ss$type_area, "urbana", "rural"
)

df_predicciones_ss$high_edu <- fct_relevel(
  df_predicciones_ss$high_edu, 
  "superior", "secundaria",  "primaria", "sin educación"
)


ggplot(df_predicciones_ss, aes(x = cohort_quinquenal, y = high_edu, fill = edad_esc_real)) +
  
  geom_tile(color = "white") +
  
  # divide matriz por área (filas) y Etnia (columnas)
  facet_grid(type_area ~ ethnic_rec) + 
  
  scale_fill_gradientn(
    colours = grunge_fun(100),
    name = "Edad \nprimer \nHijo",
    limits = c(min(df_predicciones_ss$edad_esc_real),
               max(df_predicciones_ss$edad_esc_real)),
    oob = scales::squish
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Cohorte de Nacimiento",
    y = "Nivel Educativo") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )


df_plot_ss <- df_predicciones_ss %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n_h) %>%
  mutate(order_n = row_number())


ggplot(df_plot_ss) +
  geom_point(aes(x = order_n, y = edad_esc_real),
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


##### Factores Shrinkage #####

shrinkage_ss <- strata_obs %>%
  mutate(
    gammai = (colMeans(scott_smith_mcmc$Lambdai))
  )

ggplot(shrinkage_ss, aes(x = reorder(n_h, gammai), y = gammai)) +
  geom_point() +
  labs(title = "Coeficiente de Shrinkage por Estrato - Scott-Smith", 
       x = "Estratos (ordenados por tamaño)", y = expression(gamma[i])) +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )

### Exportación de resultados ###
# Se almacenan los efectos aleatorios estimados y las medidas de precisión para su 
#uso posterior en comparación de modelos.
write.csv(resultados_ss, "resultados_estratos_ss.csv", row.names = FALSE)


# ------------------------------------------------------------------------------- #
#####  Modelo Battese-Harter-Fuller #####
y_i <- log(df_indiv$age_1st_b)
stratum_id <- df_indiv$stratum_id

# hbsae permite especificar la estructura de área por separado
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



X_bhf <- X_unit_df %>%
  dplyr::select(-stratum_id, -sam_weight_adj) %>%
  as.matrix()


modelo_hb <- fSAE.Unit(
  y = y_i,
  X = X_bhf,  
  area = stratum_id,                                                            #Variable que define los dominios
  Narea = N_vec,
  Xpop = Xpop_mat,
  method = "HB"                                                                 #el método Bayesiano Jerárquico 
)

print(modelo_hb, correlation = TRUE)



### Diagnóstico y verificación de supuestos ###
# Diagnóstico gráfico de residuos y efectos aleatorios.

# residuos del modelo
residuos <- residuals(modelo_hb)
ajustados <- fitted(modelo_hb)

# Histogramas
hist(residuos, main = "Distribución de Residuos", xlab = "Residuos de Pearson", col = "lightblue")

# Gráficos Q-Q norm
qqnorm(residuos)
qqline(residuos, col = "red")


# Residuos vs. Valores Ajustados
plot(ajustados, residuos, main = "Residuos vs. Ajustados", 
     xlab = "Valores Predichos", ylab = "Residuos")
abline(h = 0, col = "red")


#Criterios de información
# DIC
devianza <- -2 * modelo_hb$llh.c
logL <- modelo_hb$llh.c

# DIC = D_hat + 2 * número efectivo de parámetros 
DIC_bhf <- devianza + 2 * modelo_hb$p.eff
DIC_bhf


# Extracción de efectos aleatorios
efectos_u <- raneff(modelo_hb)
hist(efectos_u, main="Normalidad de Efectos Aleatorios", xlab="u_d")
qqnorm(efectos_u); qqline(efectos_u, col="blue")


### Cálculo del factor de encogimiento ###
# valores cercanos a 1 indican mayor dependencia del modelo,
# mientras que valores cercanos a 0 indican mayor peso del estimador directo.

# Extraer las varianzas del modelo
sig_u2_hb <- sv2(modelo_hb)                                                     #Varianza entre áreas
sig_e2_hb <- se2(modelo_hb)                                                     #Varianza residual
ratio <- sig_u2_hb / sig_e2_hb

rho <- ratio / (1+ratio)
rho

# Obtener el factor gamma por estrato
shrinkage_df_bhf <- df_indiv %>%
  group_by(stratum_id) %>%
  summarise(n_i = n()) %>%
  mutate(
    lambda_i = modelo_hb$gamma
  )

# Ver los primeros resultados
head(shrinkage_df_bhf)

summary(shrinkage_df_bhf$lambda_i)

shrinkage_df_bhf <- shrinkage_df_bhf %>%
  arrange(n_i) %>%
  mutate(stratum_id = as.character(stratum_id),
    order_n = row_number())

ggplot(shrinkage_df_bhf, aes(x = order_n, y = lambda_i)) +
  geom_point() +
  labs(title = "Coeficiente de Shrinkage por Estrato - BHF", 
       x = "Estratos (ordenados por tamaño)", y = expression(gamma[i])) +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )


### Precisión de las estimaciones SAE ###
# Se utilizan los errores relativos obtenidos mediante relSE(), equivalentes 
#al coeficiente de variación (CV) de las predicciones por estrato.

# Se extraen los coeficientes de variación (CV) y medidas relativas de 
#precisión para cada estrato - CV calculado sobre la escala logarítmica del modelo
cv_reales <- relSE(modelo_hb, type = "sae")
se.bhf <- SE(modelo_hb)
postVar <- raneff.se(modelo_hb)
se_u_h <- sqrt(modelo_hb$Vraneff)

# Crear el dataframe de resultados por estrato
resultados_estratos_bhf <- data.frame(
  stratum_id = as.character(modelo_hb$sampledAreaNames),
  u_h = as.numeric(efectos_u),
  se_u_h = as.numeric(se_u_h),
  CV = as.numeric(cv_reales) *100
)

resultados_estratos_bhf <- resultados_estratos_bhf %>%
  left_join(shrinkage_df_bhf, by = "stratum_id") %>%
  left_join(strats, by = "stratum_id") %>%
  arrange(n_i) %>%
  mutate(order_n = row_number())


ggplot(resultados_estratos_bhf, aes(x = order_n, y = CV)) +
  geom_point() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  labs(title = "Coeficientes de Variación por Estrato", 
       x = "Estratos (ordenados por tamaño)", y = "CV (%)") +
  theme(axis.text.x = element_text(size = 0),
        axis.title.y = element_text(size = 16)
  )



##### Visualización de medias #####
# df para predicciones
df_predicciones_bhf <- df_indiv %>%
  group_by(stratum_id, strata, cohort_quinquenal, high_edu, type_area, ethnic_rec) %>%
  summarise(n = n(), .groups = "drop")


# estimaciones media por estrato
est.bhf <- EST(modelo_hb, type = "sae")

# predicción en escala log
df_predicciones_bhf$pred_log <- est.bhf

# Corrección de sesgo por retransfomación lognormal.
# Dado que el modelo se ajusta sobre log(age_1st_b), las predicciones se devuelven 
# a la escala original mediante E(Y)=exp(eta_hat).
# donde eta_hat corresponde al predictor EBLUP obtenido en escala logarítmica.
df_predicciones_bhf <- df_predicciones_bhf %>%
  mutate(edad_real = exp(pred_log))


# re-organización de niveles para visualización
df_predicciones_bhf$ethnic_rec <- fct_relevel(
  df_predicciones_bhf$ethnic_rec,
  "otros", "indigena", "blanco/mestizo", "afro"
)

df_predicciones_bhf$type_area <- fct_relevel(
  df_predicciones_bhf$type_area, "urbana", "rural"
)

df_predicciones_bhf$high_edu <- fct_relevel(
  df_predicciones_bhf$high_edu, 
  "superior", "secundaria",  "primaria", "sin educación"
)


### Mapa de calor - combinaciones de variables ###
ggplot(df_predicciones_bhf, aes(x = cohort_quinquenal, y = high_edu, fill = edad_real)) +
  
  geom_tile(color = "white") +
  
  # divide matriz por área (filas) y Etnia (columnas)
  facet_grid(type_area ~ ethnic_rec) + 
  
  scale_fill_gradientn(
    colours = grunge_fun(100),
    name = "Edad \nprimer \nHijo",
    limits = c(min(df_predicciones_bhf$edad_real),
               max(df_predicciones_bhf$edad_real)),
    oob = scales::squish
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Cohorte de Nacimiento",
    y = "Nivel Educativo") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )



### diagrama de dispersión de la edad promedio predicha por estrato - ordenado por tamaño ###
bhf_plot <- df_predicciones_bhf %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n) %>%
  mutate(order_n = row_number())


ggplot(bhf_plot) +
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
write.csv(resultados_estratos_bhf, "resultados_estratos_bhf.csv", row.names = FALSE)
