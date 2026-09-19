# GeoFlow-SPH

**A QGIS plugin for setting up and running SPH-based numerical simulations of geophysical mass flows (landslides and debris flows).**

GeoFlow-SPH provides an end-to-end, GUI-driven workflow inside QGIS: preparing input data (terrain and initiation zone), configuring and generating the model's data files, running the SPH simulation engine, and converting the results into GIS-friendly formats for visualization and analysis.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Typical Workflow](#typical-workflow)
- [Tools](#tools)
  - [1. Create Data-Master Files](#1-create-data-master-files)
  - [2. Convert DEM Raster to TOP](#2-convert-dem-raster-to-top)
  - [3. Convert Initiation Zone to PTS](#3-convert-initiation-zone-to-pts)
  - [4. Run SPH Model](#4-run-sph-model)
  - [5. Convert RES to netCDF](#5-convert-res-to-netcdf)
- [File Formats Reference](#file-formats-reference)
- [Menu and Toolbar](#menu-and-toolbar)
- [License](#license)
- [Credits](#credits)
- [Support](#support)

---

## Overview

GeoFlow-SPH was developed in Madrid over almost a decade of research and has been applied to theoretical, experimental, and real landslide/debris-flow case histories. The underlying numerical framework implements a **two-phase, two-layer, depth-integrated model** — a general approach to two-phase mass-flow modeling capable of reproducing key physical processes such as pore-water pressure evolution, and the dynamic behavior of debris flows as a function of soil properties (permeability, volumetric stiffness, etc.).

This plugin wraps that numerical engine with an intuitive QGIS interface, so that geoscientists and engineers can go from raw GIS data (a DEM and an initiation zone) to a fully configured, running SPH simulation, and back to georeferenced results, without leaving QGIS.

The plugin was developed for the **Inter-American Development Bank (IADB)** as part of project **ES-T1343**.

| | |
|---|---|
| **Category** | Science |
| **Tags** | SPH, debris flows, landslides, geohazards, geotechnical engineering, numerical modelling, mass flows, risk assessment |
| **QGIS versions** | 3.16 – 3.44 |
| **Status** | Experimental |
| **Processing provider** | Yes |

---

## Key Features

- Graphical generation of the SPH model's `.DAT` and `.MASTER.DAT` configuration files — no manual text-file editing required.
- Conversion of raster DEMs into the model's native `.TOP` terrain format.
- Multiple ways to define the landslide/debris-flow **initiation zone** (single raster, DEM-differencing/change detection, or an existing vector points layer) and export it to `.PTS` format.
- One-click execution of the SPH simulation engine from within QGIS.
- Conversion of raw simulation results (`.QGIS_res`) into **netCDF**, ready for further analysis, time-series extraction, or visualization in QGIS, Python, or other GIS/scientific tools.
- All processing steps are also registered as standard **QGIS Processing** algorithms, so they can be run from the Processing Toolbox, chained in the **Graphical Modeler**, or scripted/batched via `qgis.processing` and PyQGIS.

---

## Requirements

- **QGIS**: 3.16 to 3.44
- **Python packages** (typically bundled with QGIS): `numpy`, `gdal`/`osgeo`
- An external **SPH simulation executable** accessible on the system (invoked via a generated Windows batch file — see [Run SPH Model](#4-run-sph-model))
- Input data:
  - A Digital Elevation Model (DEM) covering the area of interest
  - A definition of the initiation zone (raster, dual-date DEMs, or a point layer)

---

## Installation

1. Download or clone the plugin repository:
   `https://github.com/EL-BID/Deslizamientos-y-flujos`
2. Copy (or symlink) the plugin folder into your QGIS plugins directory, e.g.:
   - Windows: `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\`
   - Linux/macOS: `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
3. Restart QGIS, then enable **GeoFlow-SPH** via **Plugins → Manage and Install Plugins → Installed**.
4. A new **GeoFlow-SPH** menu appears under the **Plugins** menu (and its icon on the toolbar), giving access to all five tools described below. The same tools are also available in the **Processing Toolbox** under the *Modeling* and *Tools* groups.

---

## Typical Workflow

The tools are designed to be used in the following order:

```
DEM (raster)  ──────────────►  [2] DEM to TOP  ─────────────►  .TOP file
                                                                     │
Initiation zone data ────────► [3] Points to PTS ───────────►  .PTS file
                                                                     │
Simulation setup ────────────► [1] Create Data-Master Files ─► .DAT + .MASTER.DAT files
                                                                     │
                              .TOP + .PTS + .DAT + .MASTER.DAT
                                                                     │
                                                                     ▼
                                                      [4] Run SPH model ──► .QGIS_res
                                                                     │
                                                                     ▼
                                                [5] RES to netCDF ──► .nc (for analysis/visualization)
```

---

## Tools

### 1. Create Data-Master Files

**Menu action:** *Create Data-Master Files*
**Interface:** Custom dialog (not a Processing algorithm) with two tabs, **Data File** and **Master File**.

This is the configuration hub of the plugin. It walks the user through defining every parameter needed by the SPH engine and writes them out as the two plain-text control files the engine expects.

**Data File tab** — defines the physical/numerical model setup, written to a `.dat` file:
- **Modeling Type**: Deterministic Propagation Model, Probabilistic Propagation Model, or Probabilistic Initiation Model.
- **Problem Type**: Shallow-Water Equation or Fractional-Step Method.
- **Algorithm Type**: numerical time-integration scheme (e.g. 4th-order Runge-Kutta; additional schemes are present in the interface for future support).
- **Rheological / uncertain variables**: basal friction angle (tangent), mixture density, cohesion, Voellmy turbulent coefficient, and viscosity — each with mean and standard deviation for probabilistic runs.
- **Probabilistic seismic/height parameters** (for the Probabilistic modes): PGA (Peak Ground Acceleration) values and return periods for minor/average/major events, relative height activation, number of years, probability distribution (triangular, normal, log-normal), and number of random samples.
- **Zone and target numbers**, used to define regions of interest for probabilistic analyses.

**Master File tab** — defines the time-stepping and I/O control of the run, written to a `.MASTER.DAT` file:
- **Data file name** (must match the name used in the corresponding `.dat` file).
- **Time-step**, **total analysis time**, and **maximum number of time-steps**.
- **Print / save / plot** frequency (in steps).
- **Time-stepping technique**: Adaptive (recommended) or Constant.
- **Case description**, stored as a header comment in the file.

Clicking **Generate .MASTER.Dat File** validates all required fields and writes the `.MASTER.DAT` file to a user-chosen directory. The companion `.dat` file must be generated/named to match, as noted in the confirmation message.

> **Tip:** The `.MASTER.DAT` and `.dat` file names must correspond to each other, or the SPH engine will fail to locate the configuration.

A secondary utility dialog, reachable from within this workflow, lets the user define a **cylindrical soil-particle distribution** (center X/Y, radius, height, and number of nodes per diameter) — useful for constructing idealized/test initiation geometries.

---

### 2. Convert DEM Raster to TOP

**Menu action:** *Convert DEM Raster to TOP*
**Processing ID:** `geoflowsph:dem2top` (group: *Tools*)

Converts a standard DEM raster layer into a `.TOP` file — the terrain topography format read natively by the SPH engine.

| Parameter | Description |
|---|---|
| **DEM** | Input raster layer representing terrain elevation. |
| **Output** | Destination path for the generated `.TOP` file (`*.top` / `*.TOP`). |

---

### 3. Convert Initiation Zone to PTS

**Menu action:** *Convert Initiation Zone to PTS*
**Processing ID:** `geoflowsph:points2pts` (group: *Tools*)

Generates the `.PTS` file describing the initiation zone (the initial mass of material that will be simulated), the second key input alongside the `.TOP` terrain file. Three alternative, mutually exclusive input methods are supported:

| Method | Description |
|---|---|
| **① Single TIF** | A single GeoTIFF whose pixel values represent the elevation/z-coordinate of the initiation zone. Every valid (positive, non-NoData) pixel becomes a point. |
| **② DEM Subtraction** | A pre-event and a post-event DEM are differenced to detect elevation change above a configurable threshold. Can be filtered to erosion only, deposition only, or both — useful for deriving an initiation zone from observed pre/post-landslide topography. |
| **③ Vector Points** | An existing point vector layer is used directly, with height taken either from a numeric attribute field or from a single constant value applied to all points. |

**Common parameters** (apply to all methods):

| Parameter | Description |
|---|---|
| **Preferred X/Y spacing** | Target point spacing for the generated point set. |
| **Scaling factor for the smoothing length (k)** | SPH smoothing-length multiplier (allowed range 1.0–5.0). |
| **Output CRS** | Coordinate reference system tag associated with the output. Coordinates are written as-is (no reprojection is performed). |
| **Output** | Destination path for the generated `.PTS` file (`*.pts` / `*.PTS`). |

The algorithm validates that exactly one input method's required parameters (and, for the Vector Points method, exactly one height option) are provided before running.

---

### 4. Run SPH Model

**Menu action:** *Run SPH model*
**Processing ID:** `geoflowsph:sphsimplemode` (group: *Modeling*)

Executes the SPH simulation itself, given all the previously generated input files.

| Parameter | Description |
|---|---|
| **Problem name** | Identifier used for the simulation run and its output files. |
| **.MASTER.DAT file** | Master control file, from [Create Data-Master Files](#1-create-data-master-files). |
| **.DAT file** | Physical/numerical configuration file, from [Create Data-Master Files](#1-create-data-master-files). |
| **.PTS file** | Initiation zone points, from [Convert Initiation Zone to PTS](#3-convert-initiation-zone-to-pts). |
| **.TOP file** | Terrain topography, from [Convert DEM Raster to TOP](#2-convert-dem-raster-to-top). |
| **Output folder** | Destination folder for the simulation results. |

**What it does:**
1. Copies and arranges all inputs into a temporary working directory.
2. Generates a Windows batch file to invoke the SPH engine.
3. Executes the batch file (`cmd.exe`) and streams progress via the QGIS feedback log.
4. Copies the resulting output files to the chosen output folder.
5. Cleans up the temporary working directory.

**Output:** a `<problem_name>.QGIS_res` file in the output folder, ready for conversion to netCDF.

> **Note:** Simulation execution relies on `cmd.exe`, so this step currently requires a Windows environment with the SPH engine executable available on the system.

---

### 5. Convert RES to netCDF

**Menu action:** *Convert RES to NetCDF*
**Processing ID:** `geoflowsph:res2netcdf` (group: *Tools*)

Converts a raw `.QGIS_res` simulation result file into a georeferenced **netCDF** file, suitable for further analysis, time-series extraction, or visualization in QGIS, Python (e.g. `xarray`), or other scientific/GIS tools.

| Parameter | Description |
|---|---|
| **RES file** | The `.QGIS_res` file produced by [Run SPH Model](#4-run-sph-model). |
| **DEM** | Raster layer used as the reference terrain for georeferencing the output. |
| **Output** | Destination path for the generated `.nc` file. |

---

## File Formats Reference

| Extension | Produced by | Consumed by | Description |
|---|---|---|---|
| `.TOP` | DEM to TOP | Run SPH model | Terrain topography grid |
| `.PTS` | Points to PTS | Run SPH model | Initiation zone point set |
| `.dat` | Create Data-Master Files (Data File tab) | Run SPH model | Physical/numerical model configuration |
| `.MASTER.DAT` | Create Data-Master Files (Master File tab) | Run SPH model | Time-stepping and I/O control configuration |
| `.QGIS_res` | Run SPH model | RES to netCDF | Raw simulation results |
| `.nc` | RES to netCDF | GIS/analysis tools | Georeferenced simulation results |

---

## Menu and Toolbar

All tools are available from **Plugins → GeoFlow-SPH** in QGIS, and the plugin icon is also added to the QGIS toolbar (opening the *Create Data-Master Files* dialog by default). Every tool other than *Create Data-Master Files* is additionally registered with the **QGIS Processing framework**, so it can be found by searching the Processing Toolbox, used inside the Graphical Modeler, or called from the Python/PyQGIS console, e.g.:

```python
from qgis import processing

processing.run("geoflowsph:dem2top", {
    "INPUT": "path/to/dem.tif",
    "OUTPUT": "path/to/output.top",
})
```

---

## License

Licensed under **GPL-3.0-or-later**. See [LICENSE.md](https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/LICENSE.md) for full terms.

Copyright © 2025 Banco Interamericano de Desarrollo (BID). Authorized use under AM-331-A3.

---

## Credits

- **Original code**: (C) 2024 NaturalGIS
- **Modified and extended by**: Saeid Moussavi Tayyebi, Manuel Pastor
- **Developed for**: Inter-American Development Bank (IADB), project ES-T1343
- **Contact**: saeid.moussavita@alumnos.upm.es

## Support

- **Homepage**: https://github.com/EL-BID/Deslizamientos-y-flujos/tree/v0.1_WIP
- **Source repository**: https://github.com/EL-BID/Deslizamientos-y-flujos
- **Issue tracker**: https://github.com/EL-BID/Deslizamientos-y-flujos/issues
