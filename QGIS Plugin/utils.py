# -*- coding: utf-8 -*-

"""
***************************************************************************
    utils.py
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
import math
import shutil
from itertools import groupby
from operator import itemgetter
from tempfile import mkdtemp
from typing import Any
from datetime import datetime, timedelta

from netCDF4 import Dataset
from cftime import date2num

import numpy

from qgis.PyQt.QtCore import QProcess
from qgis.core import (
    Qgis,
    QgsRectangle,
    QgsMessageLog,
    QgsRunProcess,
    QgsRasterLayer,
    QgsBlockingProcess,
    QgsProcessingFeedback,
    QgsProcessingException,
    QgsRasterFileWriter,
    QgsRasterBlock,
    QgsProcessingFeatureSource,
    QgsFeatureRequest,
    QgsFeature,
)
from processing.core.ProcessingConfig import ProcessingConfig


PLUGIN_ROOT = os.path.dirname(__file__)
SPH_EXECUTABLE = "SPH_EXECUTABLE"


def sph_executable() -> str:
    """
    Returns path to the SHP executable.
    """
    filePath = ProcessingConfig.getSetting(SPH_EXECUTABLE)
    return filePath if filePath is not None else "sph24"


def execute(commands: list[str], feedback: QgsProcessingFeedback = None):
    """
    Executes SPH tool
    """
    if feedback is None:
        feedback = QgsProcessingFeedback()

    fused_command = " ".join([str(c) for c in commands])
    QgsMessageLog.logMessage(fused_command, "Processing", Qgis.Info)
    feedback.pushInfo("SPH command:")
    feedback.pushCommandInfo(fused_command)
    feedback.pushInfo("SPH output:")

    def onStdOut(ba: bytes):
        val = ba.data().decode("utf-8")
        if "%" in val:
            onStdOut.progress = int(progressRegex.search(val).group(0))
            feedback.setProgress(onStdOut.progress)
        else:
            onStdOut.buffer += val

        if onStdOut.buffer.endswith(("\n", "\r")):
            feedback.pushConsoleInfo(onStdOut.buffer.rstrip())
            onStdOut.buffer = ""

    onStdOut.progress = 0
    onStdOut.buffer = ""

    def onStdErr(ba: bytes):
        val = ba.data().decode("utf-8")
        onStdErr.buffer += val

        if onStdErr.buffer.endswith(("\n", "\r")):
            feedback.reportError(onStdErr.buffer.rstrip())
            onStdErr.buffer = ""

    onStdErr.buffer = ""

    command, *arguments = QgsRunProcess.splitCommand(fused_command)
    proc = QgsBlockingProcess(command, arguments)
    proc.setStdOutHandler(onStdOut)
    proc.setStdErrHandler(onStdErr)

    res = proc.run(feedback)
    if feedback.isCanceled() and res != 0:
        feedback.pushInfo("Process was canceled and did not complete.")
    elif not feedback.isCanceled() and proc.exitStatus() == QProcess.CrashExit:
        raise QgsProcessingException("Process was unexpectedly terminated.")
    elif res == 0:
        feedback.pushInfo("Process completed successfully.")
    elif proc.processError() == QProcess.FailedToStart:
        raise QgsProcessingException(
            'Process "{}" failed to start. Either "{}" is missing, or you may have insufficient permissions to run the program.'.format(
                command, command
            )
        )
    else:
        feedback.reportError("Process returned error code {}".format(res))


def generate_batch_file(problem_name: str, work_dir: str) -> str:
    """
    Generates script to run SPH tool.

    Returns a full path to the generated script.
    """
    input_file_name = "files.txt"
    input_file = os.path.join(work_dir, input_file_name)
    with open(input_file, "w", encoding="utf-8") as f:
        for i in range(2):
            f.write(f"{problem_name}\n")

    binary_name = os.path.split(sph_executable())[1]
    batch_file = os.path.join(work_dir, f"{problem_name}.bat")
    with open(batch_file, "w", encoding="utf-8") as f:
        f.write("set CWDIR=%~dp0\n")
        f.write(f"cd {work_dir}\n")
        f.write(f"{binary_name} < {input_file_name}\n")
        f.write("cd %WDIR%\n")

    return batch_file


def copy_inputs(
    problem_name: str,
    dem: str,
    pts_file: str,
    master_file: str | None = None,
    config_file: str | None = None,
) -> str:
    """
    Copies SPH executable and input files into a separate directory to create an
    environment for performing analysis.

    Returns full path to the created directory.
    """
    work_dir = mkdtemp(prefix=f"sph-")

    shutil.copy(sph_executable(), work_dir)

    new_path = os.path.join(work_dir, f"{problem_name}.top")
    shutil.copyfile(dem, new_path)

    new_path = os.path.join(work_dir, f"{problem_name}.pts")
    shutil.copyfile(pts_file, new_path)

    if master_file is not None:
        new_path = os.path.join(work_dir, f"{problem_name}.master.dat")
        shutil.copyfile(master_file, new_path)

    if config_file is not None:
        new_path = os.path.join(work_dir, f"{problem_name}.dat")
        shutil.copyfile(config_file, new_path)

    return work_dir


def generate_master_file(file_name: str, params: dict[str, Any]):
    """
    Generates configuration file with the model inputs and generic parameters.
    """
    with open(file_name, "w", encoding="utf-8") as f:
        f.write("1\n")
        f.write(f"{params['problem_name']}\n")
        f.write("sph     gfl  Monte-Carlo\n")
        f.write("  1      0      0 \n")
        f.write("problem_type   Integ_Alg \n")
        f.write("       1               4\n")
        f.write("file.dat\n")
        f.write(f"{params['problem_name']}\n")
        f.write("dt      time_end   maxtimesteps\n")
        f.write(f"{params['dt']}      {params['time_end']}          1000000\n")
        f.write("print_step   save_step  plot_step\n")
        f.write(
            f"   {params['print_step']}          {params['print_step']}       {params['print_step']}\n"
        )
        f.write("dt_sph  ic_adapt    \n")
        f.write("0.1      1  \n")
        f.write("time_curves      max_pts\n")
        f.write("     0              6  \n")
        f.write("cases_win      ic_eros\n")
        f.write("   0              0 \n")


def generate_data_file(file_name: str, params: dict[str, Any]):
    """
    Generates configuration file with specific model parameters.
    """
    with open(file_name, "w", encoding="utf-8") as f:
        f.write("  1\n")
        f.write(f"{params['problem_name']}\n")
        f.write("SWalg\n")
        f.write("  0\n")
        f.write("nhist\n")
        f.write("  0\n")
        f.write("ndimn\n")
        f.write("   2\n")
        f.write("soil  water  vps  ic_abs\n")
        f.write(" 1      0     0     0\n")
        f.write("icunk     h_inf_SW  \n")
        f.write("  6         0.1\n")
        f.write("file.pts\n")
        f.write(f"  {params['problem_name']}\n")
        f.write("pa_sph,  nnps, sle, skf \n")
        f.write("  2       2     1    1         \n")
        f.write("sum_den,  av_vel, virt_part , nor_dens  \n")
        f.write("   T         T      F            F       \n")
        f.write(
            "cgra   dens  cmanning  eros_Coef    nfrict     Tauy0   constK  visco   tanfi8  hfrict0    c11   tanfi0   .Bfact     hrelpw   Comp  \n"
        )
        f.write(
            f" {params['cgra']}   {params['dens']}   {params['cmanning']}           {params['eros_coef']}          {params['nfrict']}         {params['tauy0']}     0.0     {params['visco']}    {params['tanfi8']}    1.e-3     0.0     0.      0.0       0.0      0.001\n"
        )
        f.write("K0\n")
        f.write("  0\n")
        f.write("icpwp \n")
        f.write("  0\n")
        f.write("coarse\n")
        f.write("  0\n")
        f.write("chk_pts\n")
        f.write("  0\n")
        f.write(
            "Gid_Mask_SW   1.hs  2.disp 3.v  4.Pwb  5 eros   6.Z  7.hrel  8.hw  9.eta  10.hs+hw  11.hsat   12.Pw\n"
        )
        f.write(
            "               1      1    0    0      0        0        0       0       0     0       0          0\n"
        )
        f.write("T_change_to_W\n")
        f.write("  1.e+12    \n")


def dem_to_top(layer: QgsRasterLayer, file_path: str):
    """
    Converts a single-band raster layer representing DEM to a text format (.top)
    required by SPH tool.

    The .top format is basically a raster in the XYZ format with the custom
    header and footer.
    """
    provider = layer.dataProvider()
    width = provider.xSize()
    height = provider.ySize()
    pixel_size = layer.rasterUnitsPerPixelX()
    pixel_count = width * height
    extent = layer.extent()

    valid_count = 0

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("ictop\n")
        f.write("  10\n")
        f.write("  np      deltx\n")
        f.write(f" {pixel_count}     {pixel_size}    \n")
        f.write("Topo_x Topo_y Topo_z\n")

        for row in range(height):
            x_min = extent.xMinimum()
            x_max = extent.xMaximum()
            y_min = extent.yMaximum() - row * pixel_size
            y_max = extent.yMaximum() - (row + 1) * pixel_size
            block_extent = QgsRectangle(x_min, y_min, x_max, y_max)
            block = provider.block(1, block_extent, width, height, None)

            for col in range(width):
                x = extent.xMinimum() + (col + 0.5) * pixel_size
                y = extent.yMaximum() - (row + 0.5) * pixel_size
                if not block.isNoData(row, col):
                    valid_count += 1
                    f.write(f"{x}\t{y}\t{block.value(row, col)}\n")

        f.write("topo_props\n")
        f.write("  0\n")

    bak_name = f"{file_path}.bak"
    os.rename(file_path, bak_name)
    with open(bak_name, "r", encoding="utf-8") as bak:
        with open(file_path, "w", encoding="utf-8") as f:
            for i, line in enumerate(bak):
                if i == 3:
                    f.write(f" {valid_count}     {pixel_size}    \n")
                    continue
                f.write(line)

    os.remove(bak_name)


def points_to_pts(
    source: QgsProcessingFeatureSource,
    field_name: str | None,
    use_z: bool,
    file_path: str,
    distx: float | None = None,
    disty: float | None = None,
    factshml: float | None = None,
):
    """
    Converts a point vector layer representing unstable material to a text format (.pts)
    required by SPH tool.

    Height of the unstable material is taken either from the field_name attribute or, if use_z
    is True from the Z coordinate. In the latter case values from the field_name attribute are
    ignored.

    The PTS format is basically a raster in the XYZ format with the custom header.
    
    Args:
        source: Source feature source with points
        field_name: Name of field containing height values (if use_z is False)
        use_z: If True, use Z coordinate from geometry instead of field
        file_path: Output file path
        distx: Preferred distance between points in X direction (optional)
        disty: Preferred distance between points in Y direction (optional)
        factshml: Coordinates and the area of influence of every SPH node (2 by default)
    """
    request = QgsFeatureRequest()
    if use_z:
        request.setNoAttributes()
    else:
        request.setSubsetOfAttributes([field_name], source.fields())

    f1 = QgsFeature()
    r = QgsFeatureRequest([1])
    ok = source.getFeatures(r).nextFeature(f1)
    f2 = QgsFeature()
    r = QgsFeatureRequest([2])
    ok = source.getFeatures(r).nextFeature(f2)
    dist = f1.geometry().asPoint().distance(f2.geometry().asPoint())
    
    # Handle default values for distx and disty
    if distx < 0.01:
        distx = dist
    if disty < 0.01:
        disty = dist

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("npoin source  grid spacing  facthsml\n")
        f.write(
            f"   {source.featureCount()}           {dist}      {distx}    {disty}       {factshml}\n"
        )
        f.write("---  X   ---------  Y   -------  h  -----\n")

        for ft in source.getFeatures(request):
            p = ft.geometry().constGet()
            f.write(f"{p.x()}\t{p.y()}\t{p.z() if use_z else ft[field_name]}\n")

def copy_outputs(work_dir: str, problem_name: str, output_dir: str):
    """
    Copies output files produced by the SPH tool to the output directory.
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for suffix in ("post.msh", "post.res", "QGIS_res"):
        output_name = os.path.join(work_dir, f"{problem_name}.{suffix}")
        if os.path.exists(output_name):
            shutil.copy(output_name, output_dir)


def res_to_netcdf(res_file: str, dem: QgsRasterLayer, output: str):
    """
    COnverts QGIS_res file produced by SPH to a netCDF4 format.
    """
    pixel_size = dem.rasterUnitsPerPixelX()
    extent = dem.extent()
    crs = dem.crs()

    provider = dem.dataProvider()
    raster_width = provider.xSize()
    raster_height = provider.ySize()

    date = datetime.now()

    ds = Dataset(output, "w", format="NETCDF4", clobber=True)
    ds.description = "Landslide model"
    ds.history = f"Created {date.ctime()}"
    ds.source = "IADB Toolbox QGIS plugin"

    lat = ds.createDimension("latitude", raster_height)
    lon = ds.createDimension("longitude", raster_width)
    time = ds.createDimension("time", None)

    latitude = ds.createVariable("latitude", "f8", ("latitude",), fill_value=-9999)
    latitude.units = "degrees north"
    latitude.long_name = "latitude"
    latitude.standard_name = "latitude"
    latitude.grid_mapping = "spatial_ref"

    longitude = ds.createVariable("longitude", "f8", ("longitude",), fill_value=-9999)
    longitude.units = "degrees east"
    longitude.long_name = "longitude"
    longitude.standard_name = "longitude"
    latitude.grid_mapping = "spatial_ref"

    grid_mapping = ds.createVariable("spatial_ref", "i8")
    grid_mapping.crs_wkt = crs.toWkt()
    grid_mapping.spatial_ref = crs.toWkt()
    grid_mapping.geographic_crs_name = crs.description()
    grid_mapping.grid_mapping_name = "latitude_longitude"

    time = ds.createVariable(
        "time",
        "i8",
        ("time",),
    )
    time.units = (
        f"seconds since {datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f %z').strip()}"
    )
    time.calendar = "gregorian"
    time.long_name = "time"

    height = ds.createVariable(
        "Height",
        "f8",
        (
            "time",
            "latitude",
            "longitude",
        ),
        fill_value=-9999,
    )
    height.units = "m"
    height.positive = "up"
    height.grid_mapping = "spatial_ref"
    height.grid_mapping_name = "latitude_longitude"

    vx = ds.createVariable(
        "Vx",
        "f8",
        (
            "time",
            "latitude",
            "longitude",
        ),
        fill_value=-9999,
    )
    vx.grid_mapping = "spatial_ref"
    vx.grid_mapping_name = "latitude_longitude"

    vy = ds.createVariable(
        "Vy",
        "f8",
        (
            "time",
            "latitude",
            "longitude",
        ),
        fill_value=-9999,
    )
    vy.grid_mapping = "spatial_ref"
    vy.grid_mapping_name = "latitude_longitude"

    vavg = ds.createVariable(
        "Vavg",
        "f8",
        (
            "time",
            "latitude",
            "longitude",
        ),
        fill_value=-9999,
    )
    vavg.grid_mapping = "spatial_ref"
    vavg.grid_mapping_name = "latitude_longitude"

    for row in range(raster_height):
        y = extent.yMaximum() - (row + 0.5) * pixel_size
        latitude[row] = y

    for col in range(raster_width):
        x = extent.xMinimum() + (col + 0.5) * pixel_size
        longitude[col] = x

    data = []
    dates = []
    h = numpy.full((raster_height, raster_width), -9999)
    v_x = numpy.full((raster_height, raster_width), -9999)
    v_y = numpy.full((raster_height, raster_width), -9999)
    v_avg = numpy.full((raster_height, raster_width), -9999)
    c = 0

    with open(res_file) as f:
        sph_time = None
        prev_time = None
        block_read = False

        for line in f:
            line = line.strip()
            if line.startswith("time"):
                if sph_time is not None:
                    prev_time = sph_time
                    block_read = True

                sph_time = float(line.split()[3])
                dates.append(date + timedelta(seconds=sph_time))
                continue

            if block_read:
                data = sorted(data, key=itemgetter(1))
                values = groupby(data, itemgetter(1))
                groups = {}
                for k, v in values:
                    groups[k] = list(v)

                for row in range(raster_height):
                    y = truncate(extent.yMaximum() - (row + 0.5) * pixel_size, 1)
                    if y in groups.keys():
                        for i in groups[y]:
                            col = math.trunc((i[0] - extent.xMinimum()) / pixel_size)
                            h[row, col] = i[2]
                            v_x[row, col] = i[3]
                            v_y[row, col] = i[4]
                            v_avg[row, col] = i[5]

                c += 1
                data[:] = []

                height[c, :, :] = h
                vx[c, :, :] = v_x
                vy[c, :, :] = v_y
                vavg[c, :, :] = v_avg

                h = numpy.full((raster_height, raster_width), -9999)
                v_x = numpy.full((raster_height, raster_width), -9999)
                v_y = numpy.full((raster_height, raster_width), -9999)
                v_avg = numpy.full((raster_height, raster_width), -9999)
                block_read = False

            values = [float(v) for v in line.split()]
            data.append(values)

    time[:] = date2num(dates, units=time.units, calendar=time.calendar)

    ds.close()


def truncate(number: float, digits: int) -> float:
    decimals = len(str(number).split(".")[1])
    if decimals <= digits:
        return number
    step = 10.0**digits
    return math.trunc(step * number) / step
