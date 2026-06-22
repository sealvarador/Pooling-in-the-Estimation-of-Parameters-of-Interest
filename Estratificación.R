################################################################################
# Estratificación.R
#
# Objetivo: 
# Construir los estratos de modelado a partir de las variables
# sociodemográficas seleccionadas, identificar ceros estructurales,
# calcular estadísticas directas por estrato y asignar el identificador
# de estrato a cada individuo.
#
# El script realiza:
#   1. Generación de todas las combinaciones posibles de estratos.
#   2. Identificación de estratos observados y ceros estructurales.
#   3. Cálculo de tamaños muestrales y estimaciones directas por estrato.
#   4. Asignación de identificadores de estrato a nivel individual.
#   5. Exploración descriptiva de los estratos construidos.
#
#
# Autor: Selene Alvarado Rosario
# Fecha: Septiembre 2025
################################################################################


# Ajustar ruta de trabajo según la ubicación del proyecto
df_areas <- read.csv("complete_birth_rec.csv", stringsAsFactors = FALSE)        #Importo el dataset que se había preparado, ya con ethnic_rec


library(dplyr)
library(tidyr)
library(ggplot2)


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


####Creación de estratos ----
strata_groups <- expand.grid(                                                   #Se expande con todas las posibles combinaciones de los niveles de las covariables
  type_area = levels(df_areas$type_area),
  high_edu = levels(df_areas$high_edu),
  ethnic_rec = levels(df_areas$ethnic_rec),
  cohort_quinquenal = levels(df_areas$cohort_quinquenal)
) #224 estratos


conteos_groups <- df_areas %>%                                                  #Se cuentan las combinaciones que tienen medición en la encuesta
  group_by(type_area, high_edu, ethnic_rec, cohort_quinquenal) %>%
  summarise(n_individuos = n(), .groups = "drop")
conteos_groups #202 estratos


strata_final <- strata_groups %>%
  left_join(conteos_groups, by = c(                                             #Dejo solo los estratos con 1 o más mediciones, los otros son 0 estructurales
    "type_area", "high_edu", "ethnic_rec", "cohort_quinquenal")) %>%
  mutate(n_individuos = ifelse(is.na(n_individuos), 0, n_individuos))


strata_final$strata <- interaction(
  strata_final$type_area,
  strata_final$high_edu,
  strata_final$ethnic_rec,
  strata_final$cohort_quinquenal,
  sep = "."
)

summary(strata_final$n_individuos)
table(strata_final$n_individuos == 0)                                           #22 ceros estructurales

which(strata_final$n_individuos == 0)
structural_zeros <- strata_final[which(strata_final$n_individuos == 0),]


write.csv(structural_zeros, "structural_zeros.csv", row.names = FALSE)          #Estratos teóricamente posibles sin observaciones en la muestra


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


###### Breve exploración de los estratos #####
df_plot_ede <- strata_obs %>%
  dplyr::select(-type_area, -high_edu, -cohort_quinquenal, -ethnic_rec) %>%
  arrange(n_h) %>%                                                          
  mutate(order_n = row_number())


### Edad promedio observada por estrato - organizado por tamaño ###
ggplot(df_plot_ede) +
  geom_point(aes(x = order_n, y = y_h),
             color = "black",
             shape = "✦",                                                       #Si falla, cambia por 18
             alpha = 0.6,                                              	  
             size = 5) +                                                
  labs(x = "Estratos (ordenados por tamaño)",
       y = "Media de Edad") +
  theme(
    axis.text.x = element_text(size = 14),                            
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)  
  )


write.csv(df_indiv, "df_indiv.csv", row.names = FALSE)                          #guardo la información a nivel individuo
write.csv(strata_obs, "estratos_nacimientos.csv", row.names = FALSE)                #guardo la partición de los estratos observados
