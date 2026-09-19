# -*- coding: utf-8 -*-

"""
***************************************************************************
    provider.py
    ---------------------
    Date                 : July 2024
    Original Copyright   : (C) 2024 by NaturalGIS
    Modified by          : Saeid Moussavi Tayyebi
    Modification Date    : September 2026
    Email                : saeid.moussavita@alumnos.upm.es
***************************************************************************
*                                                                                     *
*   Copyright © 2025 Banco Interamericano de Desarrollo (BID).                        *
*   Authorized use under AM-331-A3.                                                   *
*   For license terms, see:                                                           *
*   https://github.com/EL-BID/Deslizamientos-y-flujos/blob/v0.1_WIP/LICENSE.md        *
*                                                                                     *
***************************************************************************
This script is based on the original code developed by NaturalGIS.

Modifications and additional functionality:
- Modified the original code for integration into the GeoFlow-SPH QGIS Plugin.
- Adapted the script for landslide and debris-flow simulation analysis.
- Added functionality for processing and visualizing simulation data.
"""

import os

from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtCore import QCoreApplication

from qgis.core import QgsProcessingProvider

from processing.core.ProcessingConfig import ProcessingConfig, Setting

from .algs.dem_raster_to_top import DemToTop
from .algs.initiation_zone_to_pts import PointsToPts
from .algs.res_to_netcdf import ResToNetcdf
from .algs.run_sph_model import SphSimpleMode
from .utils import PLUGIN_ROOT, SPH_EXECUTABLE, sph_executable


class GeoflowProvider(QgsProcessingProvider):
    def __init__(self):
        super().__init__()
        self.algs = list()

    def id(self):
        return "Geoflow-SPH"

    def name(self):
        return "Geoflow-SPH"

    def longName(self):
        return "Geoflow-SPH"

    def icon(self):
        return QIcon(os.path.join(PLUGIN_ROOT, "icons", "Logo1.png"))

    def load(self):
        ProcessingConfig.settingIcons[self.name()] = self.icon()
        ProcessingConfig.addSetting(
            Setting(
                self.name(),
                SPH_EXECUTABLE,
                self.tr("SPH executable"),
                sph_executable(),
                valuetype=Setting.FILE,
            )
        )
        ProcessingConfig.readSettings()
        self.refreshAlgorithms()
        return True

    def unload(self):
        ProcessingConfig.removeSetting(SPH_EXECUTABLE)

    def loadAlgorithms(self):
        self.algs = [
            DemToTop(),
            PointsToPts(),
            ResToNetcdf(),
            SphSimpleMode(),
        ]
        for a in self.algs:
            self.addAlgorithm(a)

    def supportsNonFileBasedOutput(self):
        return False

    def supportedOutputRasterLayerExtensions(self):
        return ["nc"]

    def tr(self, string):
        return QCoreApplication.translate(self.__class__.__name__, string)
