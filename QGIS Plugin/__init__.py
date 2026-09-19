# -*- coding: utf-8 -*-

"""
***************************************************************************
    __init__.py
    ---------------------
    Date                 : July 2024
    Copyright            : (C) 2024 by NaturalGIS
    Email                : info at naturalgis dot pt
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from .plugin import GeoflowPlugin
from .dat_master import Ui_GeoFlowPlugin


class CombinedPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.geoflow_plugin = GeoflowPlugin(iface)
        self.dat_master = Ui_GeoFlowPlugin(iface)

    def initGui(self):
        self.geoflow_plugin.initGui()
        self.dat_master.initGui()

    def unload(self):
        self.geoflow_plugin.unload()
        self.dat_master.unload()


def classFactory(iface):
    return CombinedPlugin(iface)