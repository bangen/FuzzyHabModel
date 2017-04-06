# Script name: CHaMP_gdbExport.py
#
# Last updated: 1/30/2017
# Created by: Sara Bangen (sara.bangen@gmail.com)
#
# Description:
#	- Outputs shapefiles and reaster from geodatabase
#
# Output:
#	- shapefile and rasters with same name as that in gdb
#
#
# Assumomptions:
#   - *gdb is named 'SurveyGeoDatabase'
#
# The user will need to specify the function arguments
#
# Args:
#   folderPath:  Parent directory folder.

# -----------------------------------
# Set user-defined  input parameters
# -----------------------------------

pf = r"C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data"
data = r"C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\FishData\ReddData\UGR\UGR_PotentialSitesWithRedds_Subsample.csv"

# -----------------------------------
# Start of script

# Import required modules
# Check out the ArcGIS Spatial Analyst extension license
import arcpy, os, fnmatch, pandas
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension('Spatial')

visits = pandas.read_csv(data)

def gdbExport(visitPath):
    print visitPath
    # Set workspace
    # Set environment settings to overwrite output

    arcpy.env.overwriteOutput = True
    env.qualifiedFieldNames = False

    # Search workspace folder for all polygon shapefiles that match searchName
    # Add to list
    gdbPath = os.path.join(visitPath, 'Topo', 'SurveyGeoDatabase.gdb')
    arcpy.env.workspace = gdbPath
    outPath = os.path.join(visitPath, 'Sims/FIS/Inputs')

    if len(arcpy.ListRasters('Water_Depth')) > 0:
        arcpy.CopyRaster_management("Water_Depth", os.path.join(outPath, "Water_Depth.tif"))
    else:
        print 'WARNING: No water depth raster for: ' + visitPath

    arcpy.FeatureClassToShapefile_conversion(["Channel_Units", "EdgeofWater_Points", "WaterExtent"], outPath)

    arcpy.env.workspace = os.path.join(visitPath, 'Topo')
    if len(arcpy.ListFiles('D50.tif')) > 0:
        arcpy.CopyRaster_management('D50.tif', os.path.join(outPath, 'D50.tif'))
    else:
        arcpy.env.workspace = outPath
        field = arcpy.ListFields('Channel_Units.shp', 'D50')
        if len(field) >= 1:
            d50_raw = arcpy.PolygonToRaster_conversion('Channel_Units.shp', 'D50', 'in_memory/d50_raw', 'CELL_CENTER', '', '0.1')
            d50 = ExtractByMask(d50_raw, 'Channel_Units.shp')
            d50.save('D50.tif')
        else:
            print 'WARNING: No D50 raster could be created for: ' + visitPath

    arcpy.Delete_management("in_memory")



visits['visit.dir'].map(gdbExport)

#visitPath = visits['visit.dir'].iloc[18]

#gdbExport(visitPath)