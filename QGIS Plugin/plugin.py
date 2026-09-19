# -*- coding: utf-8 -*-

"""
***************************************************************************
    plugin.py
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

from qgis.PyQt.QtCore import QCoreApplication, QTranslator
from qgis.core import QgsApplication

from .provider import GeoflowProvider
from .utils import PLUGIN_ROOT


class GeoflowPlugin:
    def __init__(self, iface):
        locale = QgsApplication.locale()
        qm_path = os.path.join(PLUGIN_ROOT, "i18n", f"Geoflow_SPH1_{locale}.qm")

        if os.path.exists(qm_path):
            self.translator = QTranslator()
            self.translator.load(qm_path)
            QCoreApplication.installTranslator(self.translator)

        self.provider = GeoflowProvider()

    def initProcessing(self):
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()

    def unload(self):
        if hasattr(self, 'provider') and self.provider:
            QgsApplication.processingRegistry().removeProvider(self.provider)
