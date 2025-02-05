# Geoflow-SPH: Plataforma deslizamientos y flujos


![My Image](images/Logo-BW.png)
![version](https://img.shields.io/badge/version-0.1.0-blue)

**Geoflow-SPH** es una plataforma de código abierto, desarrollada en Fortran por la **Universidad Politécnica de Madrid**, para el cálculo de deslizamientos geotécnicos y su propagación espacial, basado en el método **S**moothed **P**article **H**ydrodynamics (SPH).

# ✨ Descripción

En la plataforma que se ha almacenado, se incluyen dos modelos que permiten estudiar la iniciación o disparo del deslizamiento y su propagación

El módulo de iniciación permite a partir de los modelos geológicos geotécnicos bicapa propuestos para el área de estudio y de la caracterización de los factores disparadores, realizar el análisis de estabilidad de taludes considerando mecanismos de falla tipo traslacionales, de poca profundidad y con grandes velocidades de deformación. 
El modelo de propagación de deslizamientos   consta de un modelo integrado en profundidad, basado en un modelo tridimensional que permite describir un caso general en el que las fases (partículas sólidas y fluidos) que forman los materiales que deslizan pueden tener velocidades y desplazamientos relativos grandes, pudiendo llegarse a casos en los que el agua abandona la matriz sólida.


# 📄 Guía de Usuario

Incluir aquí la información


# 💻 Guía de Instalación

Paso a paso de cómo instalar la herramienta digital. En esta sección es recomendable explicar la arquitectura de carpetas y módulos que componen el sistema.

Según el tipo de herramienta digital, el nivel de complejidad puede variar. En algunas ocasiones puede ser necesario instalar componentes que tienen dependencia con la herramienta digital. Si este es el caso, añade también la siguiente sección.

La guía de instalación debe contener de manera específica:
- Los requisitos del sistema operativo para la compilación (versiones específicas de librerías, software de gestión de paquetes y dependencias, SDKs y compiladores, etc.).
- Las dependencias propias del proyecto, tanto externas como internas (orden de compilación de sub-módulos, configuración de ubicación de librerías dinámicas, etc.).
- Pasos específicos para la compilación del código fuente y ejecución de tests unitarios en caso de que el proyecto disponga de ellos.

# 🧑‍🍳 Autores

Geoflow-SPH es desarrollado por el **Department of Applied Mathematics, ETS Ingenieros de Caminos, Universidad Politécnica de Madrid** con apoyo del **Banco Interamericano de Desarrollo**

Equipo de desarrolladores:
Manuel Pastor, Saeid M. Tayyebi, Miguel M. Stickle, Ángel Yagüe, Miguel Molinos, Pedro Navas & Diego Manzanal

# 📚 Publicaciones

1. Pastor, M., Tayyebi, S.M., Stickle, M.M. et al. A depth integrated, coupled, two-phase model for debris flow propagation. Acta Geotech. 16, 2409–2433 (2021). https://doi.org/10.1007/s11440-020-01114-4
2. Saeid Moussavi Tayyebi, Manuel Pastor, Miguel Martin Stickle, Ángel Yagüe, Diego Manzanal, Miguel Molinos, Pedro Navas, (2022). SPH numerical modelling of landslide movements as coupled two-phase flows with a new solution for the interaction term. European Journal of Mechanics - B/Fluids. Volume 96. Pages 1-14, https://doi.org/10.1016/j.euromechflu.2022.06.002.

# 📑 Licencia

Copyright© 2025. Banco Interamericano de Desarrollo ("BID"). Uso autorizado [AM-331-A3](https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/LICENSE.md)

## Limitación de responsabilidades

El BID no será responsable, bajo circunstancia alguna, de daño ni indemnización, moral o patrimonial; directo o indirecto; accesorio o especial; o por vía de consecuencia, previsto o imprevisto, que pudiese surgir:

i. Bajo cualquier teoría de responsabilidad, ya sea por contrato, infracción de derechos de propiedad intelectual, negligencia o bajo cualquier otra teoría; y/o

ii. A raíz del uso de la Herramienta Digital, incluyendo, pero sin limitación de potenciales defectos en la Herramienta Digital, o la pérdida o inexactitud de los datos de cualquier tipo. Lo anterior incluye los gastos o daños asociados a fallas de comunicación y/o fallas de funcionamiento de computadoras, vinculados con la utilización de la Herramienta Digital.

