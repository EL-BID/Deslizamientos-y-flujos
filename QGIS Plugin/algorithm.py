# -*- coding: utf-8 -*-

"""
***************************************************************************
    algorithm.py
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

from qgis.core import QgsProcessingAlgorithm
from .utils import PLUGIN_ROOT


class GeoflowAlgorithm(QgsProcessingAlgorithm):
    def __init__(self):
        super().__init__()

    def createInstance(self):
        return type(self)()

    def icon(self):
        return QIcon(os.path.join(PLUGIN_ROOT, "icons", "Logo1.png"))

    def tr(self, text):
        return QCoreApplication.translate(self.__class__.__name__, text)
