# The PICACHO LANDSLIDE, San Salvador,  El Salvador (1982)

# ✨ Description

The Picacho landslide in San Salvador in 1982 was caused by heavy rainfall. It is estimated that around 425,000 cubic meters of material were mobilized. Infiltrating a large amount of water into the soil altered its physical and mechanical properties, ultimately causing a surface slide. This slide then transformed into a debris flow, traveling at least 4 kilometers downhill. Tragically, this event claimed the lives of an estimated 300 to 400 people.

# 📄 Instruction Guide for Using the Disaster Risk Management IADB Toolbox in QGIS

## Step 1: Install the Plugin

1. If **netCDF4** is not installed, open the **OSGeo4W Shell** from the Start menu and run the following command:
   
- **python -m pip install netCDF4**

<img width="1057" height="572" alt="image" src="https://github.com/user-attachments/assets/6a99bb85-0065-459c-9159-f1c11e4c3b5a" />

2.	Open **QGIS**.
3.	Go to **Plugins** in the top menu and click on **Manage and Install Plugins**.
4.	 Check the **Show also experimental plugins** setting
   
<img width="1382" height="961" alt="Show also experimental plugins" src="https://github.com/user-attachments/assets/d5d81434-1892-4b13-8228-acdacf4ca504" />

5.	In the search bar, type **GeoFlow-SPH**.
6.	Select the plugin and click **Install**.

<img width="1382" height="961" alt="image" src="https://github.com/user-attachments/assets/7f5aa713-7a10-417d-be20-8046177af2c1" />

## Step 2: Configure the Plugin

1.	Go to **Setting**s in the top menu and click **Options**.
2.	Navigate to the **Processing** section.
3.	Upload the required **executable file** to enable processing.

<img width="1170" height="808" alt="image" src="https://github.com/user-attachments/assets/38e348d3-cbf2-411f-bfa6-f2d21d95540a" />

## Step 3: Upload Input Files

1.	Prepare the required input files:
- **dem.tif** (Raster file containing topography coordinates)
- **PTS.shp** (Shapefile with elevation data for slope failure)
2.	Store both files in the **INPUT** folder.
3.	Open QGIS and upload these files to the project.

![7](https://github.com/user-attachments/assets/35985b07-2184-467a-a52f-8b17b7d84751)

## Step 4: Data File Configuration

1.	Go to **Plugins** in the top menu.
2.	Navigate to **GeoFlow-SPH** → **Create Data-Master Files**

<img width="708" height="250" alt="image" src="https://github.com/user-attachments/assets/8a7a3ec9-fa89-434a-b7d6-c424fdf57b68" />

3. This tab contains all parameters necessary for defining the physical simulation.

<img width="1758" height="1007" alt="image" src="https://github.com/user-attachments/assets/5c8e425b-11e9-490a-adb7-2ebd27195f49" />

4. Clicking Generate **.DAT** File validates all required fields and writes the **.DAT** file to a user-chosen directory.
5. This message appears when the .dat file has been generated successfully.

<img width="622" height="156" alt="image" src="https://github.com/user-attachments/assets/deaec939-9ac2-4296-98e1-cab5c16ff9a0" />

6. For input parameter details, see the [Input Parameters Reference](https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/QGIS%20Plugin/README.DAT.md).

## Step 5: Master File Configuration

1.	Go to **Plugins** in the top menu.
2.	Navigate to **GeoFlow-SPH** → **Create Data-Master Files**

<img width="708" height="250" alt="image" src="https://github.com/user-attachments/assets/8a7a3ec9-fa89-434a-b7d6-c424fdf57b68" />

3. The second tab controls simulation runtime parameters and multi-run configurations.

<img width="1757" height="1007" alt="image" src="https://github.com/user-attachments/assets/506a8597-93ec-4259-9780-8a40d4146531" />

4. Clicking Generate **.MASTER.DAT** File validates all required fields and writes the **.MASTER.DAT** file to a user-chosen directory.
5. This message appears when the .dat file has been generated successfully.

<img width="621" height="180" alt="image" src="https://github.com/user-attachments/assets/3fd37877-76c7-4ad7-9d6a-86491ff57da7" />

6. For simulation runtime parameters and multi-run configurations, see the [Simulation Parameters Reference](https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/QGIS%20Plugin/README.MASTER.DAT.md). 

## Step 6: Convert DEM Raster to TOP

1.	Go to **Plugins** in the top menu.
2.	Navigate to **GeoFlow-SPH** → **Convert DEM Raster to TOP**

<img width="707" height="250" alt="image" src="https://github.com/user-attachments/assets/a0c42d39-ec1d-4847-9c8a-224bed97a9ac" />

3.	Select **dem.tif** as the input file.
4.	Choose an arbitrary location to save the output.
5.	Click **Run** to generate the **.TOP** file.

<img width="857" height="631" alt="image" src="https://github.com/user-attachments/assets/95534e90-d7d0-4eea-b53f-fa824ddf2b23" />

## Step 7: Convert Initiation Zone to PTS

1.	Go to **Plugins** in the top menu.
2.	Navigate to **GeoFlow-SPH** → **Convert Initiation Zone to PTS**

<img width="707" height="251" alt="image" src="https://github.com/user-attachments/assets/e4d0503d-3f68-43dc-85ba-502d3aa0a795" />

3.	Select **PTS.shp** as the input file.
4.	Set **z** as the **Height field**.
5.	Choose an arbitrary location to save the output.
6.	Click **Run** to generate the **.PTS** file.

<img width="856" height="1105" alt="image" src="https://github.com/user-attachments/assets/ec45db87-4644-4096-8e69-61787a9da845" />

## Step 8: Run the SPH Model

1. Go to **Plugins** in the top menu.
2.	Navigate to **GeoFlow-SPH** → **Run SPH Model**.

<img width="708" height="252" alt="image" src="https://github.com/user-attachments/assets/97b0b380-aa8d-4823-b022-6032211cf0a9" />

3.	Upload the following input files:
-	**Picacho.MASTER.DAT**
-	**Picacho.DAT**
-	**Picacho.PTS**
-	**Picacho.TOP**
4.	Choose an output folder to store results.
5.	Click **Run** to start the simulation.

<img width="857" height="672" alt="image" src="https://github.com/user-attachments/assets/8cc14452-eb1d-46df-b11d-959a96c11e55" />

## Step 9: Convert Results to NetCDF Format

1. In the chosen folder where the results are stored, you can see the output files such as **Picacho.QGIS_res**.
2. Go to **Plugins** in the top menu.
3. Navigate to **GeoFlow-SPH** → **Convert RES to netCDF**.

<img width="707" height="251" alt="image" src="https://github.com/user-attachments/assets/2634a203-bf4b-45b3-9584-decee2d6dfb4" />

5. Upload the **Picacho.QGIS_res** file.
6. Convert it to **.nc format**, suitable for visualizing flow propagation in QGIS.

<img width="853" height="671" alt="image" src="https://github.com/user-attachments/assets/0adba205-39c8-40e1-8dc9-e8d04015110e" />

## Step 10: Visualize Results in QGIS

1.	Navigate to the output folder and locate **Picacho.nc**.
2.	Open **Data Source Manager**.
3.	Select **Mesh** and upload **Picacho.nc**.
4.	Add it to the **QGIS layer**.
5.	Open **Temporal Controller**.
6.	Click **Play** to visualize the simulation.

If everything has been done correctly, this animation should be displayed:

![Untitled Project6](https://github.com/user-attachments/assets/bbe687e7-bcdf-44b9-8bbf-f2771ba81d3d)
![Untitled Project1](https://github.com/user-attachments/assets/fcbd0675-bd6f-4e2d-818e-e6372b63cb9f)

**Final Notes**:
- Ensure all required input files are correctly formatted before running the tools.
- Store intermediate files in an organized manner to avoid confusion.
- If errors occur, check the processing log for troubleshooting.

This guide provides a step-by-step approach to using the Disaster Risk Management IADB Toolbox efficiently in QGIS.

The numerical analysis of the PICACHO LANDSLIDE was conducted using a straightforward one-phase model. It is evident that substantial material is flowing out of the impacted area. To achieve more accurate results for this case study, we recommend adopting the advanced TWO-PHASE TWO-LAYER MODEL. If you are interested in implementing this advanced approach, we suggest reading the article by [**Pastor et al. (2024)**](https://doi.org/10.1016/j.compgeo.2024.106581).

# Visualization in GiD Simulation

The simulation outputs two files in the directory:

-	myproblem.POST.MESH: Containing the topography and particle data.
-	myproblem.POST.RES. This file will be read by [**GiD**](https://www.gidsimulation.com/), which contains the results of all the variables calculated at various times (the times are determined by the interval specified in the MASTER.DAT file for writing to the output file).

**2D and 3D Visualization**: Using user-friendly [**GiD**](https://www.gidsimulation.com/) software, developed by CIMNE in Barcelona.

Consequently, the preprocessor GeoFlow and the postprocessor [**GiD**](https://www.gidsimulation.com/) will be utilized. The postprocessing stage represents the concluding phase of a computational model, during which the results are assessed and visualized.

In the animation below, a numerical analysis of the PICACHO LANDSLIDE is conducted using a TWO-PHASE TWO-LAYER MODEL. The results show that the propagation material can spread throughout all observed impacted areas and cover almost the entire trimline.

![2](https://github.com/user-attachments/assets/793f6ab0-4d1a-4414-99bd-0e295ca6a213)













