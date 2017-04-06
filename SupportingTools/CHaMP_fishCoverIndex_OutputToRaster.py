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
import pandas
import os
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension('Spatial')

visits = pandas.read_csv(data)

visits = visits.iloc[27:28]

def coverRaster(visitPath):
    print visitPath

    arcpy.Delete_management("in_memory")

    #  environment settings
    outPath = os.path.join(visitPath, 'Sims/FIS/Inputs')
    arcpy.env.workspace = outPath  # set workspace to pf
    arcpy.env.overwriteOutput = True  # set to overwrite output

    #  import required rasters + shapefiles
    wd = Raster('Water_Depth.tif')
    units = 'Channel_Units.shp'

    #  set raster environment settings
    desc = arcpy.Describe(wd)
    arcpy.env.extent = desc.Extent
    arcpy.env.outputCoordinateSystem = desc.SpatialReference
    arcpy.env.cellSize = desc.meanCellWidth

    #  --- cover index raster --

    # read in cover index csv
    inTblCI = arcpy.ListFiles('FuzzyCoverIndex.csv')[0]
    arr = arcpy.da.TableToNumPyArray(inTblCI, '*')
    arcpy.da.NumPyArrayToFeatureClass(arr, 'in_memory/coverIndex_Pts', ('x', 'y'), desc.SpatialReference)
    arcpy.PointToRaster_conversion('in_memory/coverIndex_Pts', 'FuzzyCoverIndex', 'in_memory/coverIndex', 'MEAN', '', 0.1)
    coverIndex = ExtractByMask('in_memory/coverIndex', units)
    coverIndex.save('CoverIndex.tif')
    arcpy.Delete_management('in_memory/coverIndex_Pts')

    #  --- fuzzy hsi inputs csv --

    inTblFI = arcpy.ListFiles('FuzzyHSI_Inputs.csv')[0]
    print inTblFI
    fieldnames = [field.name for field in arcpy.ListFields(inTblFI)]
    print fieldnames
    # fieldtypes = [field.type for field in arcpy.ListFields(inTblFI)]
    # print fieldtypes
    fieldList = arcpy.ListFields(inTblFI)
    for field in fieldList:  # loop through each field
        if field.type != 'Double':  # look for the name elev
            print field.name

    if 'D50' in fieldnames:
        arr2 = arcpy.da.TableToNumPyArray(inTblFI, ['x', 'y', 'Vel', 'Depth', 'D50'])
        arcpy.da.NumPyArrayToFeatureClass(arr2, 'in_memory/fuzzyInputs_Pts', ('x', 'y'), desc.SpatialReference)
        # extract cover index values to points
        ExtractMultiValuesToPoints('in_memory/fuzzyInputs_Pts', [['CoverIndex.tif', 'CoverIndex']],
                                   'NONE')
    else:
        arr2 = arcpy.da.TableToNumPyArray(inTblFI, ['x', 'y', 'Vel', 'Depth'])
        arcpy.da.NumPyArrayToFeatureClass(arr2, 'in_memory/fuzzyInputs_Pts', ('x', 'y'), desc.SpatialReference)
        # extract cover index values to points
        ExtractMultiValuesToPoints('in_memory/fuzzyInputs_Pts', [['D50.tif', 'D50'], ['CoverIndex.tif', 'CoverIndex']], 'NONE')

    # remove na values (i.e., values < 0)
    with arcpy.da.UpdateCursor('in_memory/fuzzyInputs_Pts', ['D50', 'CoverIndex']) as cursor:
        for row in cursor:
            if row[0] < 0.0 or row[1] < 0.0:
                cursor.deleteRow()

    #print os.path.join(arcpy.env.workspace, 'tmp_fuzzyInputs_Pts.shp')
    #arcpy.CopyFeatures_management('in_memory/fuzzyInputs_Pts', os.path.join(arcpy.env.workspace, 'tmp_fuzzyInputs_Pts.shp'))

    #  convert to numpy array and save as txt file
    nparr = arcpy.da.FeatureClassToNumPyArray('in_memory/fuzzyInputs_Pts', ['SHAPE@X', 'SHAPE@Y', 'Vel', 'Depth', 'D50', 'CoverIndex'])
    numpy.savetxt(os.path.join(arcpy.env.workspace, 'FuzzyHSI_Inputs.csv'), nparr, fmt="%s", delimiter=",", header = str('x,y,Vel,Depth,D50,CoverIndex'), comments='')


visits['visit.dir'].map(coverRaster)

#visitPath = visits['visit.dir'].iloc[10]
#coverRaster(visitPath)