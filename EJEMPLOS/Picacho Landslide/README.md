# The PICACHO LANDSLIDE, San Salvador,  El Salvador (1982)

# ✨ Description

The Picacho landslide in San Salvador in 1982 was caused by heavy rainfall. It is estimated that around 425,000 cubic meters of material were mobilized. Infiltrating a large amount of water into the soil altered its physical and mechanical properties, ultimately causing a surface slide. This slide then transformed into a debris flow, traveling at least 4 kilometers downhill. Tragically, this event claimed the lives of an estimated 300 to 400 people.

# 📄 Instruction Guide for Using the Disaster Risk Management IADB Toolbox in QGIS

## Step 1: Install the Plugin

1. If **netCDF4** is not installed, open the **OSGeo4W Shell** from the Start menu and run the following command:
   
- **python -m pip install netCDF4**

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

![3](https://github.com/user-attachments/assets/1bf7c4ef-0b10-4309-91a1-805ad7d79381)

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

4. Clicking Generate .DAT File validates all required fields and writes the .DAT file to a user-chosen directory.
5. This message appears when the .dat file has been generated successfully.

<img width="622" height="156" alt="image" src="https://github.com/user-attachments/assets/deaec939-9ac2-4296-98e1-cab5c16ff9a0" />

## Step 5: Master File Configuration

1.	Go to **Plugins** in the top menu.
2.	Navigate to **GeoFlow-SPH** → **Create Data-Master Files**

<img width="708" height="250" alt="image" src="https://github.com/user-attachments/assets/8a7a3ec9-fa89-434a-b7d6-c424fdf57b68" />

3. The second tab controls simulation runtime parameters and multi-run configurations.

<img width="1757" height="1007" alt="image" src="https://github.com/user-attachments/assets/506a8597-93ec-4259-9780-8a40d4146531" />

4. Clicking Generate .MASTER.DAT File validates all required fields and writes the .MASTER.DAT file to a user-chosen directory.
5. This message appears when the .dat file has been generated successfully.

<img width="621" height="180" alt="image" src="https://github.com/user-attachments/assets/3fd37877-76c7-4ad7-9d6a-86491ff57da7" />

## Step 4: Convert DEM to TOP Format

1.	Open the **Processing Toolbox**.
2.	Navigate to **IADB** → **Tools** → **DEM to TOP**.

<img width="707" height="255" alt="image" src="https://github.com/user-attachments/assets/cacfb362-f22c-4564-af20-474f864553ce" />

3.	Select **dem.tif** as the input file.
4.	Choose an arbitrary location to save the output.
5.	Click **Run** to generate the **.TOP** file.

![9](https://github.com/user-attachments/assets/49b928e1-7578-4bd7-a972-44ca26afc78d)

## Step 5: Convert Points to PTS Format

1.	Open the **Processing Toolbox**.
2.	Navigate to **IADB** → **Tools** → **Points to PTS**.
3.	Select PTS.shp as the input file.
4.	Set **z** as the **Height field**.
5.	Choose an arbitrary location to save the output.
6.	Click **Run** to generate the **.PTS** file.

![10](https://github.com/user-attachments/assets/603cc544-c0ee-4935-b946-482acb90337e)


## Step 6: Run the SPH Model

1.	Open the **Processing Toolbox**.
2.	Navigate to **IADB** → **Modeling** → **SPH model (simple mode)**.
3.	Upload the following input files:
-	**Picacho.MASTER.dat**
-	**Picacho.dat**
-	**Picacho.pts**
-	**Picacho.top**
4.	Choose an output folder to store results.
5.	Click **Run** to start the simulation.

![11](https://github.com/user-attachments/assets/a0e0b1e2-0036-4c65-9ab9-19a2acc9150e)

## Step 7: Convert Results to NetCDF Format

1. In the chosen folder where the results are stored, you can see the output files such as **Picacho.QGIS_res**.
2. Open the **Processing Toolbox**.
3. Navigate to **IADB** → **Tool**s → **RES to NetCDF**.
4. Upload the **Picacho.QGIS_res** file.
5. Convert it to .nc format, suitable for visualizing flow propagation in QGIS.

  ![5](https://github.com/user-attachments/assets/19500d47-fd08-4332-a708-8fade9a366fa)

## Step 8: Visualize Results in QGIS

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













