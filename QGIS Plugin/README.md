# Disaster Risk Management IADB Toolbox

Processing provider that integrates various disaster risk management tools into QGIS.
about=The toolbox was developed for IADB ([Inter-American Development Bank](https://iadb.org)) as a part of the ES-T1343 project (Platforma innovadora para la reduccion del riesgo de deslizamentos y flujos de detritos en El Salvador).

For now the plugin allows to run a tool called SPH (Smooth Particle Hydrodynamics) developed at the Escuela de Ingenieros de Caminos, Canales y Puertos of the University of Madrid by prof. Manuel Tomas Pastor Perez and his team. The SPH tool allow to model surface quick landslides.

## Usage

The plugin is available QGIS Python plugins repository.

After installation of the plugin it is necessary to configure it, namely to set the setting "SPH executable" ("Settings → Options", then switch to the "Processing" tab and expand "Providers" section) so it points to the executable file of the SPH tool. The SPH tool has to be downloaded and installed separately.

It is also necessary to install required dependencies, namely netCDF4 Python package. Windows users can install it using the OSGeo4W installer.

### Tools

The plugin provides several tools:

 * **DEM to TOP** — helper tool used to convert raster layers representing DEM (Digital Elevation Model) from commong GIS formats to a `.TOP` file required by the SPH tool. This tool shold be used to prepare input DEM for **SPH model** tools.
 * **RES to netCDF**  — helper tool used to convert one of the outputs produced by the SPH tool into a netCDF format. The file produced by this tool then can be opened in QGIS.
 * **SPH model (simple mode)** — runs the SPH tool using pre-made configuration and input files. Please note that input DEM should be provided in a `.TOP` format, if your DEM is in a common GIS format, e.g. GeoTiff, use the **DEM to TOP** tool to prepare input file. This tool creates several files in the output directory, but those files can not be opened in QGIS directly, please use **RES to netCDF** tool to convert `.QGIS_res` file to a netCDF format.
 * **SPH model (advanced mode)** — allows power users to tweak various parameters affecting modelling process and then run the SPH tool. Please note that input DEM should be provided in a `.TOP` format, if your DEM is in a common GIS format, e.g. GeoTiff, use the **DEM to TOP** tool to prepare input file. This tool creates several files in the output directory, but those files can not be opened in QGIS directly, please use **RES to netCDF** tool to convert `.QGIS_res` file to a netCDF format.

If necessary **DEM to TOP**, **SPH model** and **RES to netCDF** tools can be combined together in a model, so that preparation of inputs and conversion of outputs is done automatically.
