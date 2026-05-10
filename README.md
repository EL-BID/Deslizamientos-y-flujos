# Geoflow-SPH: Plataforma deslizamientos y flujos

![Mi Imagen](images/Logo-BW.png)
![version](https://img.shields.io/badge/version-0.1.0-blue)

> 🌍 **Read this in [English](README.en.md)**

**Geoflow-SPH** es una plataforma de código abierto, desarrollada en Fortran por la **Universidad Politécnica de Madrid** (UPM), para el cálculo de deslizamientos geotécnicos y su propagación espacial, basado en el método **S**moothed **P**article **H**ydrodynamics (SPH).

# ✨ Descripción General

Este repositorio contiene el código para simular deslizamientos de tierra y flujos de detritos combinando datos geoespaciales y técnicas de modelado numérico. La simulación tiene como objetivo proporcionar información sobre la dinámica de los deslizamientos, incluyendo la iniciación, el movimiento y la deposición de sedimentos. Esta herramienta es de gran utilidad para investigadores, geólogos, profesionales de la gestión de desastres y responsables de políticas públicas.

## Características

- **Modelado de Iniciación de Deslizamientos:** Estimación de la probabilidad de rotura de taludes bajo diversos escenarios y obtención de la Probabilidad de Falla (PoF).
- **Dinámica de Flujos y Detritos:** Modelado de la propagación y deposición de flujos de tierra y detritos.
- **Modelado Estocástico:** Un marco que utiliza el modelo dinámico GeoFlow junto con modelos probabilísticos (Monte Carlo o FOSM) para predecir extensiones e intensidades potenciales en regiones con parámetros reológicos y geotécnicos inciertos.
- **Visualización:** Visualizaciones en 2D y 3D de los eventos simulados.

# 💻 Requisitos e Instalación

Los requisitos de hardware para utilizar este programa son mínimos; cualquier ordenador relativamente moderno (Pentium III o superior) es suficiente. Tampoco existen grandes limitaciones respecto al sistema operativo, ya que tanto FORTRAN como GiD son compatibles con Windows y Linux.

**Geoflow-SPH** no requiere instalación, ya que funciona como un archivo ejecutable portátil. Los usuarios pueden ejecutar directamente el ejecutable precompilado disponible en la carpeta [Executive file](<EJEMPLOS/Picacho Landslide/Executive file/>). Los usuarios avanzados que deseen realizar cambios en los archivos fuente pueden seguir estas instrucciones para compilar sus propios ejecutables:

**1. Instalar Visual Studio 2022**
- Visite el sitio web de Visual Studio.
- Descargue la edición Visual Studio 2022 Community.
- Durante la instalación, asegúrese de seleccionar la carga de trabajo "Desarrollo para el escritorio con C++". Esto es necesario ya que el compilador Intel Fortran se integra con las herramientas de compilación de C++.

**2. Instalar el Compilador Fortran**
- Obtenga el [**Intel® oneAPI HPC Toolkit**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.mjodbt/).
- Descargue Intel® Fortran Essentials.
- Durante la instalación, asegúrese de seleccionar el componente "Intel® Fortran Compiler".

**3. Compilar el Código Fortran**
- Compile el código utilizando el compilador Intel Fortran para generar un archivo .exe.

**4. Visualización**

**Opción principal (recomendada): QGIS**

Los resultados pueden visualizarse utilizando QGIS, un sistema de información geográfica de código abierto y de libre acceso.

**Flujo de trabajo**
- GeoFlow actúa como **preprocesador**
QGIS se utiliza para **posprocesamiento y visualización**

La etapa de posprocesamiento representa la fase final del flujo de trabajo computacional, en la cual los resultados de la simulación son analizados e interpretados.

**Opción alternativa: GiD Simulation**

La visualización también puede realizarse mediante GiD (versión 7.5 o superior), una herramienta de pre y posprocesamiento desarrollada por CIMNE. Para visualizar los resultados, se recomienda el uso del software [GiD Simulation](https://www.gidsimulation.com/), versión 7.5 o superior.

- Soporta visualización **2D y 3D**
- Proporciona una interfaz gráfica para posprocesamiento avanzado

# 📄 Guía del Usuario

Este documento no es un manual exhaustivo del código GeoFlow, sino una guía para facilitar el aprendizaje progresivo de su uso.

## Preparación de Datos de Entrada
El programa requiere un conjunto de archivos de datos ASCII. Asumiendo que el nombre del problema es "miproyecto", los archivos de entrada incluirían:

- `miproyecto.MASTER.DAT`: Parámetros del cálculo numérico (algoritmo, pasos de tiempo), datos de control y escenarios probabilísticos.
- `miproyecto.TOP`: Detalles topográficos e identificación de áreas especiales con propiedades únicas.
- `miproyecto.DAT`: Información esencial para modelar el problema (partículas virtuales, parámetros de control, propiedades de materiales, etc.).
- `miproyecto.PTS`: Descripción de la masa movilizada (utilizado en la fase de propagación).

## Ejecución de la Simulación
- Coloque los archivos de datos y el ejecutable (*.exe) en el mismo directorio.
- Al ejecutar el programa, escriba el nombre del problema (miproyecto) cuando se le solicite.
- Luego, el programa pedirá el nombre del archivo de topografía (no es necesario incluir la extensión).

## Resultados (Output)
La simulación genera dos archivos principales:
- `miproyecto.NC`: Un archivo NetCDF que contiene una dimensión temporal y una estructura de rejilla espacial (es decir, datos tipo raster) puede importarse y visualizarse directamente en **QGIS**.
- `miproyecto.POST.MESH`: Contiene la topografía y los datos de las partículas.
- `miproyecto.POST.RES`: Archivo de resultados para ser leído por [**GiD**](https://www.gidsimulation.com/), con las variables calculadas en los intervalos especificados.

# 👓 Ejemplos

Para facilitar la adopción de la plataforma, se proporciona un ejemplo de un deslizamiento histórico importante. El objetivo es que los usuarios exploren cómo se incorporan las características detalladas en los archivos de entrada. El siguiente ejemplo se encuentra en la carpeta [EJEMPLOS](https://github.com/EL-BID/Deslizamientos-y-flujos/tree/v0.1_WIP/EJEMPLOS):

- Deslizamiento de El Picacho, El Salvador (1982)

# 🧑‍🍳 Autores

Geoflow-SPH ha sido desarrollado por el **Departamento de Matemática Aplicada, ETS Ingenieros de Caminos, Universidad Politécnica de Madrid** con el apoyo del **Banco Interamericano de Desarrollo**.

Equipo de desarrollo:
Manuel Pastor ([@ManuelPastor53](https://github.com/ManuelPastor53)) y Saeid Moussavi Tayyebi ([@SaeidMT](https://github.com/SaeidMT))

# 📚 Publicaciones

1. Pastor, M., Tayyebi, S. M., Stickle, M. M., Yagüe, Á., Molinos, M., Navas, P., & Manzanal, D. (2021). A depth integrated, coupled, two-phase model for debris flow propagation. Acta Geotechnica, 16(8), 2409–2433. doi: 10.1007/s11440-020-01114-4
2. Tayyebi, S. M., Pastor, M., & Stickle, M. (2021). Two-phase SPH numerical study of pore-water pressure effect on debris flows mobility: Yu Tung debris flow. Computers and Geotechnics, 132, 103973. doi: 10.1016/j.compgeo.2020.103973
3. Pastor, M., Tayyebi, S. M., Stickle, M. M., Molinos, M., Yague, A., Manzanal, D., & Navas, P. (2022). An Arbitrary Lagrangian Eulerian (ALE) finite difference (FD)‐SPH depth integrated model for pore pressure evolution on landslides over erodible terrains. International Journal for Numerical and Analytical Methods in Geomechanics, 46(6), 1127–1153. doi: 10.1002/nag.3339
4. Tayyebi, Saeid M., Pastor, M., Stickle, M. M., Yagüe, Á., Manzanal, D., Molinos, M., & Navas, P. (2022). Two-phase SPH modelling of a real debris avalanche and analysis of its impact on bottom drainage screens. Landslides, 19(2), 421–435. doi: 10.1007/s10346-021-01772-9
5. Tayyebi, S. M., Pastor, M., Stickle, M. M., Yagüe, Á., Manzanal, D., Molinos, M., & Navas, P. (2022). SPH numerical modelling of landslide movements as coupled two-phase flows with a new solution for the interaction term. European Journal of Mechanics - B/Fluids, 96, 1–14. doi: 10.1016/j.euromechflu.2022.06.002
6. Tayyebi, S. M., Pastor, M., Hernandez, A., Gao, L., Stickle, M. M., Yifru, A. L., & Thakur, V. (2022). Two-Phase Two-Layer Depth-Integrated SPH-FD Model: Application to Lahars and Debris Flows. Land, 11(10), 1629. doi: 10.3390/land11101629
7. Pastor, M., Tayyebi, S. M., Hernandez, A., Gao, L., Stickle, M. M., & Lin, C. (2023). A new two-layer two-phase depth-integrated SPH model implementing dewatering: Application to debris flows. Computers and Geotechnics, 153, 105099. doi: 10.1016/j.compgeo.2022.105099
8. Tayyebi, S. M., Pastor, M., Yifru, A. L., Thakur, V. K. S., & Stickle, M. M. (2023). Two-phase SPH–FD depth-integrated model for debris flows: application to basal grid brakes. Géotechnique, 73(3), 218–233. doi: 10.1680/jgeot.21.00080
9. Pastor, M., Hernández, A., Tayyebi, S. M., Trejos, G. A., Suárez, G., & Zheng, J. (2024). A depth‐integrated SPH framework for slow landslides. International Journal for Numerical and Analytical Methods in Geomechanics, 48(16), 3848–3875. doi: 10.1002/nag.3814
10. Pastor, M., Tayyebi, S. M., Hernández, A., Zheng, J., Suárez, G., & Reyes, M. E. (2024). Modeling fast flows with variable water content: A depth-integrated SPH approach. Computers and Geotechnics, 174, 106581. doi: 10.1016/j.compgeo.2024.106581


# 📑 Licencia

Copyright© 2025. Banco Interamericano de Desarrollo ("BID"). Uso autorizado [AM-331-A3](https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/LICENSE.md)

## Limitación de responsabilidad

El Banco Interamericano de Desarrollo no será responsable, bajo circunstancia alguna, de daño ni indemnización, moral o patrimonial; directo o indirecto; accesorio o especial; o por vía de consecuencia, previsto o imprevisto, que pudiese surgir:

i. Bajo cualquier teoría de responsabilidad, ya sea por contrato, infracción de derechos de propiedad intelectual, negligencia o bajo cualquier otra teoría; y/o

ii. A raíz del uso de la Herramienta Digital, incluyendo, pero sin limitación de potenciales defectos en la Herramienta Digital, o la pérdida o inexactitud de los datos de cualquier tipo. Lo anterior incluye los gastos o daños asociados a fallas de comunicación y/o fallas de funcionamiento de computadoras, vinculados con la utilización de la Herramienta Digital.
