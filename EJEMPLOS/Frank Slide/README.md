# The PICACHO LANDSLIDE, San Salvador,  El Salvador (1982)

# ✨ Description

The Picacho landslide in San Salvador in 1982 was caused by heavy rainfall. It is estimated that around 425,000 cubic meters of material were mobilized. Infiltrating a large amount of water into the soil altered its physical and mechanical properties, ultimately causing a surface slide. This slide then transformed into a debris flow, traveling at least 4 kilometers downhill. Tragically, this event claimed the lives of an estimated 300 to 400 people.

# 📄 Instruction Guide for Using the Disaster Risk Management IADB Toolbox in QGIS

## Step 1: Install the Plugin

1.	Open **QGIS**.
2.	Go to **Plugins** in the top menu and click on **Manage and Install Plugins**.
3.	In the search bar, type **Disaster Risk Management IADB Toolbox**.
4.	Select the plugin and click **Install**.
5.	
![1](https://github.com/user-attachments/assets/a7a9d026-82d6-4172-b5b6-d072cd039831)

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


## Step 4: Convert DEM to TOP Format

1.	Open the **Processing Toolbox**.
2.	Navigate to **IADB** → **Tools** → **DEM to TOP**.
3.	Select **dem.tif** as the input file.
4.	Choose an arbitrary location to save the output.
5.	Click **Run** to generate the **.TOP** file.

![8](https://github.com/user-attachments/assets/e4aabf21-423d-49f1-8de0-6fad3cd34120)
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

## Step 7: Visualize Results in QGIS

1.	Navigate to the output folder and locate **Picacho.nc**.
2.	Open **Data Source Manager**.
3.	Select **Mes**h and upload **Picacho.nc**.
4.	Add it to the **QGIS layer**.
5.	Open **Temporal Controller**.
6.	Click **Play** to visualize the simulation.

![10](https://github.com/user-attachments/assets/7c7cb14e-0579-4a6c-a1f8-13f5dbc6933d)


![Untitled Project4](https://github.com/user-attachments/assets/b0665a18-f24f-4d7e-ac89-100d069de9c3)


**Final Notes**:
- Ensure all required input files are correctly formatted before running the tools.
- Store intermediate files in an organized manner to avoid confusion.
- If errors occur, check the processing log for troubleshooting.

This guide provides a step-by-step approach to using the Disaster Risk Management IADB Toolbox efficiently in QGIS.










