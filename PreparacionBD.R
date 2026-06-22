################################################################################
# PreparacionBD.R
#
# Objetivo:
# Preparar la base de datos DHS Colombia 2015 para el análisis de edad al
# primer hijo. El script realiza:
#   1. Importación de la base original.
#   2. Selección de variables de interés.
#   3. Limpieza y transformación de variables.
#   4. Exportación de la base final.
#
# Autor: Selene Alvarado Rosario
# Fecha: Septiembre 2025
################################################################################

library(dplyr)
library(foreign)
library(lubridate)


# Ajustar ruta de trabajo según la ubicación del proyecto
ruta_datos <- "Colombia/COIR72DT/COIR72FL.DTA"

if (!file.exists(ruta_datos)) {
  stop(
    paste(
      "No se encontró el archivo COIR72FL.DTA.\n",
      "Los datos DHS no pueden distribuirse en este repositorio.\n",
      "Solicite acceso en https://dhsprogram.com/data/\n",
      "y ubique el archivo en Colombia/COIR72DT/."
    )
  )
}

COIR <- read.dta(ruta_datos)


##### Limpieza de la base original  ##### 

#1. 2015####

vars_keep <- c(                                                                 ##Lista de las variables a mantener
  "caseid",                                                                     #id de individuo para la encuesta
  "sdepto",                                                                     #departamento de residencia
  "ssubreg",                                                                    #subregión de residencia
  "v005",                                                                       #peso de muestreo
  "v009", "v010",                                                               #Año y mes de nacimiento de la madre
  "v101",                                                                       #Región de residencia
  "v102",                                                                       #tipo de área de residencia
  "v106",                                                                       #Máximo nivel educativo cursado
  "v131",                                                                       #etnia a la que pertenece
  "v212")                                                                       #edad a la que tuvo su primer hijo

COIR <- COIR %>% 
  select(any_of(vars_keep))  %>%                                                #Conservar únicamente variables relevantes para el análisis
  select(where(~ !all(is.na(.)))) %>%                                           #Eliminación de las columnas que solo contengan NA
  transmute(                                                                    #Construyo y renombro las columnas según la variable
    caseid = gsub(" ", "", caseid),                                             #Eliminación de espacios
    sam_weight = v005,
    sdepto,
    ssubreg,
    birth_date = as.Date(paste(v010, v009, "01", sep = "-")),                   #Consolidación de fecha completa (une año y mes)
    age_1st_b = v212,                                                           #respuesta
    region = v101,
    type_area = v102,
    high_edu = v106,
    ethnic = v131
  ) %>%
  select(
    caseid, sam_weight, sdepto, ssubreg, birth_date, age_1st_b, 
    everything()
  )


#write.csv(COIR,"2015_def.csv")                                                  #guardado de la base limpia



##### Estructuración de la base para el estudio  #####

COIR <- COIR %>%                                                                #Elimino los registros que no tengan información de la variable de interés
  filter(!is.na(age_1st_b)) %>%
  mutate(across(where(is.character), as.factor)) %>%                            #cambio de todas las variables no numéricas a categóricas
  mutate(caseid = as.character(caseid))                                         #corrección de columnas que no son variables categóricas, vuelven a string


COIR$birth_year <- year(COIR$birth_date)


df_areas <- COIR %>%                                                            #Defino el data frame que se va a usar en el modelado
  select(caseid, sam_weight, age_1st_b, sdepto, region, ssubreg, type_area, high_edu, ethnic, birth_year)


# Base final para el modelado SAE.
# Cada fila representa una mujer entrevistada y contiene:
# - Variable respuesta: edad al primer hijo.
# - Variables auxiliares sociodemográficas.
# - Identificadores geográficos.


write.csv(df_areas, "complete_birth.csv", row.names = FALSE)                    #guardado del data set final
