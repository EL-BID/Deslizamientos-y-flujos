# Geoflow-SPH: Landslide and Debris Flow Simulation 


![My Image](images/Logo-BW.png)
![version](https://img.shields.io/badge/version-0.1.0-blue)

**Geoflow-SPH** es una plataforma de código abierto, desarrollada en Fortran por la **Universidad Politécnica de Madrid**, para el cálculo de deslizamientos geotécnicos y su propagación espacial, basado en el método **S**moothed **P**article **H**ydrodynamics (SPH).

# ✨ Overview

This repository contains code for simulating landslides and debris flows using a combination of geospatial data and numerical modeling techniques. The simulation aims to provide insights into the dynamics of landslides, including the initiation, movement, and deposition of debris. This tool can be useful for researchers, geologists, disaster management professionals, practitioners, and policymakers.

## Features

- Landslide Initiation Modeling: Estimating the likelihood of slope failure under various scenarios and obtaining the Probability of Failure (PoF).
- Flow slide and Debris Flow Dynamics: Model the propagation and deposition of Flow slides and debris flows.
- Stochastic modelling: A framework utilizing a dynamic model, GeoFlow, and a probability model, Monte Carlo or FOSM, aims to predict potential run-out extents and intensities in regions with uncertain rheological and geotechnical parameters.
- Visualization: 2D and 3D visualizations of the simulated events.


# 💻 Requirements and Installation

The hardware requirements for using this program are minimal, which is one of its greatest advantages; any relatively modern computer (Pentium III or newer) is sufficient. There are also no major limitations regarding the operating system, as both FORTRAN and GID are compatible with Windows and Linux.

**Geoflow-SPH** itself requires no installation, as it runs as a portable executable file. Users can directly run the precompliled executable available in *INSERT LOCATION OF FILE IN REPO*. More advanced users who want to make changes to the source files, can follow these instructions to compile their own executables:

**1. Install Visual Studio 2022**
-	Go to the Visual Studio website.
-	Download the Visual Studio 2022 Community edition (or any other edition you prefer).
-	During installation, ensure you select the Desktop development with C++ workload. This is necessary because the Intel Fortran compiler integrates with the C++ build tools in Visual Studio.

**2. Install Fortran Compiler**
-	Go to the Get the [**Intel® oneAPI HPC Toolkit**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.mjodbt/).
-	Download the Intel® Fortran Essentials.
-	During installation, ensure you select the Intel® Fortran Compiler component. This is essential for compiling Fortran code.

**3. Compile Fortran Code**
-	Compile the Fortran code using the Intel Fortran compiler and generate an .exe file.

**4. Visualization**

Visualization of results is done externally from GeoFlow. To view the results, users are directed to [**GiD Simulation**](https://www.gidsimulation.com/) software, version 7.5 or higher.


# 📄 User Guide

Please be aware that this is not designed to be a comprehensive manual for the GeoFlow code, but rather a guide to facilitate learning its usage progressively.

## Preparing Input Data

The program requires a set of ASCII data files to be prepared before running a problem. The data files share a problem name, which we will assume to be “myproblem”. In this hypothetical scenario, the input files would include:

-	myproblem.MASTER.DAT: Containing the parameters of the numerical calculation (type of algorithm and time steps), the control data, and those related to the probabilistic scenarios.
-	myproblem.TOP: Containing topographic details and identification of potential special areas that possess unique properties distinct from the norm.
-	myproblem.DAT: Containing the essential information to model the problem, including the form and parameters of virtual particles, program control parameters, material properties, etc.
-	myproblem.PTS: Containing the description of the mobilized mass. The file is not used during the initiation phase but is utilized in the propagation phase.

## Running the Simulation
-	Place your input data files and the executable file (*.exe) in the directory.
-	Upon running the program (*.exe), the initial screen displays a note instructing that typing the problem name (myproblem) is required to begin.
-	Then, the program will request the name of the topography data file. There is no need to enter the file extension, as the program takes it by default.

## Output

The simulation outputs two files in the directory:

-	myproblem.POST.MESH: Containing the topography and particle data.
-	myproblem.POST.RES. This file will be read by [**GiD**](https://www.gidsimulation.com/), which contains the results of all the variables calculated at various times (the times are determined by the interval specified in the MASTER.DAT file for writing to the output file).

## Visualization
•	2D and 3D Visualization: Using user-friendly [**GiD**](https://www.gidsimulation.com/) software, developed by CIMNE in Barcelona.

Consequently, the pre-processor GeoFlow and the post-processor [**GiD**](https://www.gidsimulation.com/) will be utilized. The post-processing stage represents the concluding phase of a computational model, during which the results are assessed and visualized.


# 👓 Examples 

To assist in the adoption of the platform, examples of important historic landslides with their detailed information and input files are provided. 
The end goal of these examples is for the users to be able to explore how detailed characteristics are incorporated in the input files and how they are reflected in the obtained results. The following examples are included in the [EJEMPLOS](https://github.com/EL-BID/Deslizamientos-y-flujos/tree/v0.1_WIP/EJEMPLOS) folder:

- The Frank Slide, Alberta, Canada (1903)
- The Tsing Shan Debris Flow, Hong Kong (1990)
- El Picacho Landslide, El Salvador (1982)

# Case study: The Frank Slide

<img src="https://github.com/user-attachments/assets/374a571f-2250-42ac-b141-be9423561b15" alt="1-ezgif com-crop" width="500" height="400"><img src="https://github.com/user-attachments/assets/37b7c4df-11ed-4630-a077-d5d7c9648228" alt="1-ezgif com-crop" width="500" height="400">

The numerical analysis of the Frank Slide is performed through the one-phase SPH model.

# 🧑‍🍳 Authors

Geoflow-SPH is developed by the **Department of Applied Mathematics, ETS Ingenieros de Caminos, Universidad Politécnica de Madrid** with the support of the **Interamerican Development Bank**

Development team:
Manuel Pastor ([@ManuelPastor53](https://github.com/ManuelPastor53)) and Saeid Moussavi Tayyebi ([@SaeidMT](https://github.com/SaeidMT))

# 📚 Publications

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

# 📑 License

Copyright© 2025. Banco Interamericano de Desarrollo ("BID"). Uso autorizado [AM-331-A3](https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/LICENSE.md)

## Limitation of responsibility

El BID no será responsable, bajo circunstancia alguna, de daño ni indemnización, moral o patrimonial; directo o indirecto; accesorio o especial; o por vía de consecuencia, previsto o imprevisto, que pudiese surgir:

i. Bajo cualquier teoría de responsabilidad, ya sea por contrato, infracción de derechos de propiedad intelectual, negligencia o bajo cualquier otra teoría; y/o

ii. A raíz del uso de la Herramienta Digital, incluyendo, pero sin limitación de potenciales defectos en la Herramienta Digital, o la pérdida o inexactitud de los datos de cualquier tipo. Lo anterior incluye los gastos o daños asociados a fallas de comunicación y/o fallas de funcionamiento de computadoras, vinculados con la utilización de la Herramienta Digital.

