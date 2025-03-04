# The PICACHO LANDSLIDE, San Salvador,  El Salvador (1982)

# ✨ Description

The Picacho landslide in San Salvador in 1982 was caused by heavy rainfall. It is estimated that around 425,000 cubic meters of material were mobilized. Infiltrating a large amount of water into the soil altered its physical and mechanical properties, ultimately causing a surface slide. This slide then transformed into a debris flow, traveling at least 4 kilometers downhill. Tragically, this event claimed the lives of an estimated 300 to 400 people.

# 📄 Instruction Guide for Using the Disaster Risk Management IADB Toolbox in QGIS

## Step 1: Install the Plugin

1.	Open **QGIS**.
2.	Go to **Plugins** in the top menu and click on **Manage and Install Plugins**.
3.	In the search bar, type **Disaster Risk Management IADB Toolbox**.
4.	Select the plugin and click **Install**.

## Step 2: Configure the Plugin

1.	Go to **Setting**s in the top menu and click **Options**.
2.	Navigate to the **Processing** section.
3.	Upload the required **executable file** to enable processing.

## Step 3: Upload Input Files

1.	Prepare the required input files:
- **dem.tif** (Raster file containing topography coordinates)
- **PTS.shp** (Shapefile with elevation data for slope failure)
2.	Store both files in the **INPUT** folder.
3.	Open QGIS and upload these files to the project.

## Step 4: Convert DEM to TOP Format

1.	Open the **Processing Toolbox**.
2.	Navigate to **IADB** → **Tools** → **DEM to TO**P.
3.	Select **dem.tif** as the input file.
4.	Choose an arbitrary location to save the output.
5.	Click **Run** to generate the **.TOP** file.

## Step 5: Convert Points to PTS Format

1.	Open the **Processing Toolbox**.
2.	Navigate to **IADB** → **Tools** → **Points to PTS**.
3.	Select PTS.shp as the input file.
4.	Set **z** as the **Height field**.
5.	Choose an arbitrary location to save the output.
6.	Click **Run** to generate the **.PTS** file.

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

**Final Notes**:
- Ensure all required input files are correctly formatted before running the tools.
- Store intermediate files in an organized manner to avoid confusion.
- If errors occur, check the processing log for troubleshooting.
This guide provides a step-by-step approach to using the Disaster Risk Management IADB Toolbox efficiently in QGIS.










