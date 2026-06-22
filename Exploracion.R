################################################################################
# Exploracion.R
#
# Objetivo: Realizar la exploración descriptiva de la variable respuesta y de las
# variables auxiliares candidatas para el modelado SAE.
#
#   El script realiza:
#   1. Corrección y homologación de categorías.
#   2. Construcción de cohortes de nacimiento.
#   3. Exploración gráfica de la variable respuesta.
#   4. Comparación entre grupos mediante pruebas no paramétricas.
#   5. Agrupación de categorías étnicas poco frecuentes.
#
# Autor: Selene Alvarado Rosario
# Fecha: Septiembre 2025
################################################################################


# Ajustar ruta de trabajo según la ubicación del proyecto


library(dplyr)
library(ggplot2)
library(FSA)
library(forcats)


df_areas <- read.csv("complete_birth.csv", stringsAsFactors = FALSE)            #Importo el dataset que se había preparado


df_areas <- df_areas %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(caseid = as.character(caseid))

### Corrección de caracteres especiales ###
# La codificación del archivo DHS genera caracteres especiales incompatibles con
# algunos procedimientos posteriores

df_areas$sdepto <- recode(df_areas$sdepto, 
                          "bogot?" = "bogota",
                          "atl?ntico" = "atlantico",
                          "bol?var" = "bolivar",
                          "boyac?" = "boyaca",
                          "c?rdoba" = "cordoba",
                          "caquet?" = "caqueta",
                          "choc?" = "choco",
                          "guain?a" = "guainia",
                          "nari?o" = "narino",
                          "quind?o" = "quindio",
                          "san andr?s y providencia" = "san andres",
                          "vaup?s" = "vaupes")

df_areas$ssubreg <- recode(df_areas$ssubreg,
                           "antioquia sin medellin" = "Antioquia s/Medellin",
                           "atlantico, san andres, bolivar norte" = "Caribe Norte",
                           "barranquilla a. m." = "Barranquilla AM",
                           "bogota" = "Bogota DC",
                           "bolivar sur, sucre, cordoba" = "Caribe Sur",
                           "boyaca, cmarca, meta" = "Centro-Oriente",
                           "caldas, risaralda, quindio" = "Eje Cafetero",
                           "cali a.m." = "Cali AM",
                           "cauca y nari?o sin litoral"  = "Cauca-Narino s/litoral",
                           "guajira, cesar, magdalena" = "Caribe nororiental",
                           "litoral pacifico" = "Pacifico litoral",
                           "medellin a.m." = "Medellin AM",
                           "orinoquia y amazonia" = "Orinoquia-Amazonia",
                           "santanderes" = "Santanderes",
                           "tolima, huila, caqueta" = "Sur Andino",
                           "valle sin cali ni litoral" = "Valle interior s/Cali")

df_areas$ethnic <- recode(df_areas$ethnic, 
                          "black/mulato/afro-colombian/afro-descendant" = "afro/negro",
                          "indigenous" = "indigena",
                          "none of the above" = "blanco/mestizo",
                          "gypsy (rom)" = "romani",
                          "palanquero from san basilio" = "palenquero",
                          "raizal from archipelago (san andres)" = "raizal")

df_areas$high_edu <- recode(df_areas$high_edu, 
                            "higher" = "superior",
                            "secondary" = "secundaria",
                            "no education" = "sin educación",
                            "primary" = "primaria")

df_areas$type_area <- recode(df_areas$type_area,
                             "urban" = "urbana",
                             "rural" = "rural")

# Verificación de la distribución
table(df_areas$birth_year)


# La última cohorte (1995–2002) agrupa todos los nacimientos posteriores para 
# garantizar tamaños muestrales adecuados e incluir toda la información
# disponible en la encuesta.
df_areas$cohort_quinquenal <- cut(
  df_areas$birth_year,
  breaks = c(1964, 1969, 1974, 1979, 1984, 1989, 1994, 2002),
  labels = c(
    "1965–1969",
    "1970–1974",
    "1975–1979",
    "1980–1984",
    "1985–1989",
    "1990–1994",
    "1995–2002"                                                             
  ),                                                                        
  right = TRUE,
  include.lowest = TRUE
)

# Verificación del tamaño de las cohortes construidas
table(df_areas$cohort_quinquenal)



##### Exploración y estadística descriptiva  ##### 
### definición como factores de las variables categóricas ###

df_areas$type_area <- factor(
  df_areas$type_area,
  levels = c("urbana", "rural")
)

df_areas$high_edu <- factor(
  df_areas$high_edu,
  levels = c("superior", "secundaria", "primaria", "sin educación")
)


##### Reordenamiento de categorías #####
# Las categorías se ordenan según la mediana observada de la edad al primer hijo
#para facilitar la interpretación visual de los gráficos

df_areas$ethnic <- reorder(df_areas$ethnic, df_areas$age_1st_b, median)
df_areas$ethnic <- forcats::fct_rev(df_areas$ethnic)

df_areas$sdepto <- reorder(df_areas$sdepto, df_areas$age_1st_b, median)
df_areas$sdepto <- forcats::fct_rev(df_areas$sdepto)

df_areas$ssubreg <- reorder(df_areas$ssubreg, df_areas$age_1st_b, median)
df_areas$ssubreg <- forcats::fct_rev(df_areas$ssubreg)

df_areas$region <- reorder(df_areas$region, df_areas$age_1st_b, median)
df_areas$region <- forcats::fct_rev(df_areas$region)


### número de grupos geográficos/territoriales ###
m_region <- nlevels(df_areas$region)
m_region

m_subreg <- nlevels(df_areas$ssubreg)
m_subreg

m_depto <- nlevels(df_areas$sdepto)
m_depto

#número de individuos
n <- nrow(df_areas)
n


##### Visualizaciones de los datos #####
### definición de la paleta ###
n_col <- length(unique(df_areas$sdepto))

lux_pal <- c( "#606C38", "#BC6C25", "#6C584C", "#3B596C",
              "#3F4812", "#DDA15E", "#A98467", "#32406C") 

grunge_fun <- colorRampPalette(
  c("#2F3E2E",  # verde bosque oscuro
    "#3A5A40",  # verde musgo
    "#588157",  # sage
    "#A3B18A",  # verde grisáceo
    "#B7B7A4",  # gris oliva
    "#7F5539",  # café tierra
    "#5E503F")  # marrón oscuro
)

grunge_cols <- grunge_fun(n_col)


### distribución de la variable respuesta ###
ggplot(df_areas, aes(x = age_1st_b)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 18, fill = "darkolivegreen4", color = "cornsilk1", alpha = 0.7) +
  geom_density(color = "darkolivegreen", lwd = 1.2) +
  
  labs(
    x = "Edad",
    y = "Densidad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18)
  )



#### edad al primer hijo según variables sociodemográficas ####
### edad al primer hijo según tipo de área de residencia ###
ggplot(df_areas, aes(x = type_area, y = age_1st_b, fill = type_area)) +
  geom_jitter(width = 0.35, alpha = 0.06, size = 0.6, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal)  +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


### edad al primer hijo según máximo grado educativo alcanzado por la madre ###
ggplot(df_areas, aes(x = high_edu, y = age_1st_b, fill = high_edu)) +
  geom_jitter(width = 0.35, alpha = 0.06, size = 0.6, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal)  +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


### edad al primer hijo según la etnia de la madre ###
ggplot(df_areas, aes(x = ethnic, y = age_1st_b, fill = ethnic)) +
  geom_jitter(width = 0.35, alpha = 0.06, size = 0.6, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal)  +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14,
                               angle = 30,
                               hjust = 1,
                               vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


### edad al primer hijo según la cohorte de nacimiento de la madre ###
ggplot(df_areas, aes(x = cohort_quinquenal, y = age_1st_b, fill = cohort_quinquenal)) +
  geom_jitter(width = 0.35, alpha = 0.05, size = 0.5, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal) +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12,
                               angle = 30,
                               hjust = 1,
                               vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )



##### Exploración de variables territoriales (no incluidas en el modelado final) #####
### edad al primer hijo según el departamento de residencia ###
ggplot(df_areas, aes(x = sdepto, y = age_1st_b, fill = sdepto)) +
  geom_jitter(width = 0.35, alpha = 0.05, size = 0.5, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = grunge_cols)  +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11,
                               angle = 45,
                               hjust = 1,
                               vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


### edad al primer hijo según la subregión de residencia ###
ggplot(df_areas, aes(x = ssubreg, y = age_1st_b, fill = ssubreg)) +
  geom_jitter(width = 0.35, alpha = 0.05, size = 0.5, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = grunge_cols)  +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11,
                               angle = 45,
                               hjust = 1,
                               vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


### edad al primer hijo según la región de residencia ###
ggplot(df_areas, aes(x = region, y = age_1st_b, fill = region)) +
  geom_jitter(width = 0.35, alpha = 0.05, size = 0.5, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal) +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


### edad al primer hijo según cohorte de nacimiento de la madre y su región de residencia ###
ggplot(df_areas, aes(x = region, y = age_1st_b, fill = cohort_quinquenal)) +
  geom_jitter(width = 0.35, alpha = 0.05, size = 0.5, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal) +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


##### Unión de categorias étnicas #####
### verificación de si hay diferencias significativas entre niveles de cada variable ###
kruskal.test(age_1st_b ~ ethnic, data = df_areas)                               #diferencias significativas, pero solo hay 13 romanies, y 43 palenquero
kruskal.test(age_1st_b ~ high_edu, data = df_areas)                             #diferencias significativas, sigue siendo necesario confirmar por pares
kruskal.test(age_1st_b ~ type_area, data = df_areas)                            #diferencias significativas, no se agrupa


## verificación de diferencias por pares ##
dunnTest(age_1st_b ~ ethnic, data = df_areas, method="bonferroni")              #se decide negro+palenquero+afro = negro/afro, romani+raizal = otras minorias
dunnTest(age_1st_b ~ high_edu, data = df_areas, method="bonferroni")            #no se agrupa, diferencias significativas 


## unión de categorías etnicas
df_areas <- df_areas %>%
  mutate(
    ethnic_rec = case_when(
      ethnic %in% c("afro/negro", "palenquero") ~ "afro",
      ethnic %in% c("raizal", "romani") ~ "otros",
      ethnic == "blanco/mestizo" ~ "blanco/mestizo",
      ethnic == "indigena" ~ "indigena",
      TRUE ~ NA_character_
    ),
    ethnic_rec = as.factor(ethnic_rec)
  )


## reorganización para ver los grupos por decrecimiento de la mediana ##
df_areas$ethnic_rec <- reorder(df_areas$ethnic_rec, df_areas$age_1st_b, median)
df_areas$ethnic_rec <- forcats::fct_rev(df_areas$ethnic_rec)


### edad al primer hijo según etnia de la madre luego de la unión ###
ggplot(df_areas, aes(x = ethnic_rec, y = age_1st_b, fill = ethnic_rec)) +
  geom_jitter(width = 0.35, alpha = 0.06, size = 0.6, color = "#83845B") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 1, color = "grey20") +
  
  scale_fill_manual(values = lux_pal)  +
  
  labs(
    x = "",
    y = "Edad"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14,
                               angle = 30,
                               hjust = 1,
                               vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )


kruskal.test(age_1st_b ~ region, data = df_areas)  #significativo, igual confirmo por pares
kruskal.test(age_1st_b ~ ssubreg, data = df_areas) #significativo, igual confirmo por pares
kruskal.test(age_1st_b ~ sdepto, data = df_areas) #significativo, igual confirmo por pares

dunnTest(age_1st_b ~ region, data = df_areas, method="bonferroni") 
dunnTest(age_1st_b ~ ssubreg, data = df_areas, method="bonferroni") #no se agrupa
##### finalmente las variables territoriales no se incluyen en el modelado final.


write.csv(df_areas, "complete_birth_rec.csv", row.names = FALSE)
