# Script name: CHaMP_fishCoverIndex_Rasters.py (v 1.0)
#
# Last updated: 1/9/2017
# Created by: Sara Bangen (sara.bangen@gmail.com)
#
# Assumptions:
#   Folder path(s) contains:
#       - WettedDepth.tif
#       - Channel_Units.shp
#       - FuzzyHSI_Inputs.csv

# Output:
#   Cover Index Raster
#   FIS Inputs *.csv

# -----------------------------------
# Set user-defined  input parameters
# -----------------------------------

data = r"C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\FISValidation\ChinookSpawner\UGR_ValidationSites.csv"

# -----------------------------------
# Start of script

#  import required modules and extensions
import numpy
import os
import pandas
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension('Spatial')


visits = pandas.read_csv(data)

#visits = visits.iloc[4:]

def hsiRaster(visitPath):

    print visitPath

    arcpy.Delete_management("in_memory")

    #  environment settings
    outPath = os.path.join(visitPath, 'Sims/FIS/Output')
    arcpy.env.workspace = outPath  # set workspace to pf
    arcpy.env.overwriteOutput = True  # set to overwrite output

    #  import required rasters + shapefiles
    d50 = Raster(visitPath + '/Sims/FIS/Inputs/D50.tif')

    #  set raster environment settings
    desc = arcpy.Describe(d50)
    arcpy.env.extent = desc.Extent
    arcpy.env.outputCoordinateSystem = desc.SpatialReference
    arcpy.env.cellSize = desc.meanCellWidth

    #  --- create raster --

    #  read in fuzzy model output csv
    inTbls = arcpy.ListFiles('FuzzyChinookSpawner_*.csv')
    print inTbls
    for inTbl in inTbls:
        print inTbl
        fields = arcpy.ListFields(inTbl)

        for field in fields:
            print field.name

        arr = arcpy.da.TableToNumPyArray(inTbl, '*')
        arcpy.da.NumPyArrayToFeatureClass(arr, 'in_memory/fuzzyPts', ('x', 'y'), desc.SpatialReference)
        arcpy.PointToRaster_conversion('in_memory/fuzzyPts', 'FuzzyHSI', 'in_memory/fuzzyHSI', 'MEAN', '', 0.1)
        fuzzyHSI = ExtractByMask('in_memory/fuzzyHSI', d50)
        outName = os.path.join(os.path.splitext(os.path.basename(inTbl))[0] + '.tif')
        fuzzyHSI.save(outName)
        arcpy.Delete_management("in_memory")

    # arr = arcpy.da.TableToNumPyArray(inTbls[0], '*')
    # arcpy.da.NumPyArrayToFeatureClass(arr, 'in_memory/fuzzyPts', ('x', 'y'), desc.SpatialReference)
    # arcpy.PointToRaster_conversion('in_memory/fuzzyPts', 'Vel', 'in_memory/Vel', 'MEAN', '', 0.1)
    # fuzzyHSI = ExtractByMask('in_memory/Vel', d50)
    # fuzzyHSI.save(r'C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\TestRuns\SmallSites\OJD03458-000016\2012\VISIT_733\Sims\FIS\Inputs\Vel.tif')
    #
    # arcpy.PointToRaster_conversion('in_memory/fuzzyPts', 'Depth', 'in_memory/Depth', 'MEAN', '', 0.1)
    # fuzzyHSI = ExtractByMask('in_memory/Depth', d50)
    # fuzzyHSI.save(r'C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\TestRuns\SmallSites\OJD03458-000016\2012\VISIT_733\Sims\FIS\Inputs\Depth.tif')

visits['visit.dir'].map(hsiRaster)
