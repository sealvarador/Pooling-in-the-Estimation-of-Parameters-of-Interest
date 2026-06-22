# Estimation en Áreas Pequeñas de la Edad al Primer Hijo en Colombia

Implementación de modelos clásicos, espaciales y bayesianos para la estimación de áreas pequeñas utilizando información de la Encuesta Nacional de Demografía y Salud (ENDS/DHS).


## Descripción

Este repositorio contiene los scripts utilizados para la construcción de estratos sociodemográficos y la estimación de la edad al primer hijo mediante modelos de áreas pequeñas (SAE).
Los enfoques implementados incluyen:

- Fay-Herriot clásico.
- Battese-Harter-Fuller (Nested Error Regression).
- Fay-Herriot espacial (SAR).
- Fay-Herriot espacial Bayesiano Jerárquico.

Las estimaciones se realizan sobre la escala logarítmica de la edad al primer hijo y posteriormente se retransforman a la escala original para su interpretación.


## Requisitos

El proyecto fue desarrollado en:
- R >= 4.4.0

Paquetes principales:
```r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "lubridate",
  "sae",
  "saeHB",
  "saeHB.spatial",
  "spdep",
  "spatialreg",
  "Matrix",
  "MASS",
  "fastDummies",
  "patchwork",
  "coda",
  "mvtnorm",
  "forcats"
))
```


## Datos

Los datos utilizados en este proyecto provienen del programa Demographic and Health Surveys (DHS).
Debido a las restricciones de licencia del DHS Program, los archivos originales no pueden ser distribuidos en este repositorio.
Los investigadores interesados pueden solicitar acceso directamente en:

https://dhsprogram.com/data/

Una vez aprobada la solicitud, descargue el archivo:

COIR72FL.DTA

y ubíquelo en:

Colombia/COIR72DT/

antes de ejecutar los scripts de preparación de datos.


## Orden de ejecución

1. PreparacionBD.R
2. Exploracion.R
3. Estratificacion.R
4. Modelado_Frecuentista.R
5. Modelado_Bayesiano.R
6. Modelado_Espacial.R

Los scripts deben ejecutarse en este orden para reproducir completamente los resultados.


## Metodologia

El flujo de trabajo consiste de:

1. Preparación y limpieza de los datos.
2. Construcción de estratos sociodemográficos.
3. Cálculos de estimadores directos y varianzas de muestreo.
4. Ajuste de modelos SAE frecuentistas, bayesianos y espaciales.
5. Comparación de modelos usando medidas de precisión y diagnóstico.
6. Obtención de estimaciones e indicadores de incertidumbre para las áreas pequeñas.


## Salidas

Los scripts generan:

- Estimaciones SAE por estrato.
- Errores estándar.
- Coeficientes de variación.
- Factores de contracción (shrinkage).
- Gráficos diagnósticos.


## Citación

Si utiliza este repositorio, por favor cite:

Alvarado Rosario, Selene (2026).
*Evaluación del Efecto del Pooling en la Estimación de Parámetros de Interés: Aplicación a la Edad Promedio de las Mujeres al Tener su Primer Hijo en Colombia*
Trabajo de grado, Universidad Nacional de Colombia.


## Licencia

Este repositorio es para uso académico e investigativo únicamente.




##################################-----------------------------------------------------------------
# Small Area Estimation of Age at First Birth in Colombia

Implementation of classical, spatial, and Bayesian small area estimation models using data from the Colombian Demographic and Health Survey (DHS).


## Description

This repository contains the scripts used for the construction of sociodemographic strata and the estimation of women's age at first birth using Small Area Estimation (SAE) models.
The implemented approaches include:

* Classical Fay–Herriot model.
* Battese–Harter–Fuller (Nested Error Regression) model.
* Spatial Fay–Herriot (SAR) model.
* Hierarchical Bayesian Spatial Fay–Herriot model.

Estimation is performed on the logarithmic scale of age at first birth and subsequently back-transformed to the original scale for interpretation.


## Requirements

The project was developed using:

* R >= 4.4.0

Main packages:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "lubridate",
  "sae",
  "saeHB",
  "saeHB.spatial",
  "spdep",
  "spatialreg",
  "Matrix",
  "MASS",
  "fastDummies",
  "patchwork",
  "coda",
  "mvtnorm",
  "forcats"
))
```


## Data

The data used in this project come from the Demographic and Health Surveys (DHS) Program.
Due to DHS licensing restrictions, the original survey files cannot be distributed through this repository.
Researchers interested in reproducing the analysis may request access directly from:

https://dhsprogram.com/data/

Once access has been granted, download the file:

COIR72FL.DTA

and place it in:

Colombia/COIR72DT/

before running the data preparation scripts.


## Execution Order

1. PreparacionBD.R
2. Exploracion.R
3. Estratificacion.R
4. Modelado_Frecuentista.R
5. Modelado_Bayesiano.R
6. Modelado_Espacial.R

The scripts should be executed in this order to fully reproduce the results.


## Methodology

The workflow consists of:

1. Data preparation and cleaning.
2. Construction of sociodemographic strata.
3. Calculation of direct estimators and sampling variances.
4. Fitting of classical, spatial, and Bayesian SAE models.
5. Model comparison using precision and diagnostic measures.
6. Production of small area estimates and uncertainty indicators.


## Outputs

The scripts generate:

* Small area estimates for each stratum.
* Standard errors.
* Coefficients of variation.
* Shrinkage factors.
* Diagnostic plots.


## Citation

If you use this repository, please cite:

Alvarado Rosario, Selene (2026).
*Evaluación del Efecto del Pooling en la Estimación de Parámetros de Interés: Aplicación a la Edad Promedio de las Mujeres al Tener su Primer Hijo en Colombia.*
Undergraduate Thesis, Universidad Nacional de Colombia.


## License

This repository is provided for academic and research purposes only.
