# Script name: CHaMP_fishCoverIndex_Rasters.py (v 1.0)
#
# Last updated: 12/12/2016
# Created by: Sara Bangen (sara.bangen@gmail.com)
#
# Assumptions:
#   Folder path(s) contains:
#       - Water_Depth.tif
#       - Channel_Units.shp (with field: 'Unit_Numbe')
#       - WaterExtent.shp
#       - EdgeofWater_Points.shp (with field: 'DESCRIPTIO' ['rw' + 'lw'])
#       - LargeWood*.xlsx (with fields: 'LargeWoodT' ['dry' + 'wet'], 'SumLWDCoun', 'ChannelU_1']
#       - UndercutBank*.xlsx (with fields: 'Bank' ['Right' + 'Left'], ChannelU_1']
#
# Output:
#   3 fish cover rasters:
#       - distance to deep pools
#       - distance to undercut bank
#       - large wood density (pieces per 100 m2)
#   CSV of values for each raster + [x,y] coords

# -----------------------------------
# Set user-defined  input parameters
# -----------------------------------

folderPath = r'C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\TestRuns\LargeSites\YFI00001-001971\2016\VISIT_4180\Sims\FIS\Inputs'
widthCategory = 16.0

# -----------------------------------
# Start of script

#  import required modules and extensions
import numpy
import os
import arcpy
from arcpy.sa import *
arcpy.CheckOutExtension('Spatial')


def main():

    arcpy.Delete_management("in_memory")

    #  environment settings
    arcpy.env.workspace = folderPath  # set workspace to pf
    arcpy.env.overwriteOutput = True  # set to overwrite output

    #  import required rasters + shapefiles
    wd = Raster('Water_Depth.tif')
    units = 'Channel_Units.shp'
    wePoly = 'WaterExtent.shp'
    wePoints = 'EdgeofWater_Points.shp'

    #  set raster environment settings
    desc = arcpy.Describe(wd)
    arcpy.env.extent = desc.Extent
    arcpy.env.outputCoordinateSystem = desc.SpatialReference
    arcpy.env.cellSize = desc.meanCellWidth

    #  --- lwd count --
    arcpy.MakeFeatureLayer_management(units, 'units_lyr')

    #  read in wood xlsx
    inTblLWD = arcpy.ListFiles('LargeWood*')[0]
    arcpy.ExcelToTable_conversion(inTblLWD, 'in_memory/tbl_lw')

    if 'Piece' in inTblLWD:
        #  remove 'dry' rows
        with arcpy.da.UpdateCursor('in_memory/tbl_lw', 'LargeWoo_1') as cursor:
            for row in cursor:
                if row[0] == 'Dry':
                    cursor.deleteRow()

        arcpy.AddField_management('in_memory/tbl_lw', 'Count', 'SHORT')
        with arcpy.da.UpdateCursor('in_memory/tbl_lw', 'Count') as cursor:
            for row in cursor:
                row[0] = 1
                cursor.updateRow(row)

        arcpy.Statistics_analysis('in_memory/tbl_lw', 'in_memory/tbl_lw2', [['Count', 'SUM']], 'ChannelU_1')

        # join wood count to units
        arcpy.JoinField_management('units_lyr', 'Unit_Numbe', 'in_memory/tbl_lw2', 'ChannelU_1', 'SUM_Count')

        #  add/calculate lwd density field
        arcpy.AddField_management('units_lyr', 'lwDen100m2', 'DOUBLE')
        fields = ['SHAPE@AREA', 'SUM_Count', 'lwDen100m2']
        with arcpy.da.UpdateCursor('units_lyr', fields) as cursor:
            for row in cursor:
                row[2] = float(row[1] / row[0]) * 100
                cursor.updateRow(row)
    else:
        #  remove 'dry' rows
        with arcpy.da.UpdateCursor('in_memory/tbl_lw', 'LargeWoodT') as cursor:
            for row in cursor:
                if row[0] == 'Dry':
                    cursor.deleteRow()

        # join wood count to units
        arcpy.JoinField_management('units_lyr', 'Unit_Numbe', 'in_memory/tbl_lw', 'ChannelU_1', 'SumLWDCoun')

        #  add/calculate lwd density field
        arcpy.AddField_management('units_lyr', 'lwDen100m2', 'DOUBLE')
        fields = ['SHAPE@AREA', 'SumLWDCoun', 'lwDen100m2']
        with arcpy.da.UpdateCursor('units_lyr', fields) as cursor:
            for row in cursor:
                row[2] = float(row[1] / row[0]) * 100
                cursor.updateRow(row)

    #  convert to raster
    arcpy.PolygonToRaster_conversion('units_lyr', 'lwDen100m2', 'in_memory/lwDensity', 'CELL_CENTER', '', 0.1)
    lwClip = ExtractByMask('in_memory/lwDensity', units)
    lwClip.save('lwDensity.tif')

    lw = Raster('lwDensity.tif')

    #  --- distance to deep pool raster --

    #  select all pool units
    arcpy.SelectLayerByAttribute_management('units_lyr', 'NEW_SELECTION', """ "Tier1" = 'Slow/Pool' """)

    if int(arcpy.GetCount_management('units_lyr').getOutput(0)) > 0:
        #  get max water depth value for each pool unit
        ZonalStatisticsAsTable('units_lyr', 'Unit_Numbe', wd, 'in_memory/tbl_wd', 'DATA', 'MAXIMUM')

        #  join max value back to pool units
        arcpy.JoinField_management('units_lyr', 'Unit_Numbe', 'in_memory/tbl_wd', 'Unit_Numbe', 'MAX')

        #  remove pools units with depth < 80 cm
        if widthCategory <= 10.0:
            arcpy.SelectLayerByAttribute_management('units_lyr', 'REMOVE_FROM_SELECTION', """ "MAX" < 0.43 """)
        else:
            arcpy.SelectLayerByAttribute_management('units_lyr', 'REMOVE_FROM_SELECTION', """ "MAX" < 0.86 """)

        if int(arcpy.GetCount_management('units_lyr').getOutput(0)) > 0:
            #  calculate euclidean distance
            poolDist = EucDistance('units_lyr')
        else:
            poolDist = Con(lw >= 0, 1000)
    else:
        poolDist = Con(lw >= 0, 1000)

    #  clip to wetted extent and save output
    poolDistClip = ExtractByMask(poolDist, units)
    poolDistClip.save('DeepPoolDist.tif')

    #  clear selection on units lyr
    arcpy.SelectLayerByAttribute_management('units_lyr', "CLEAR_SELECTION")

    #  --- undercuts --

    #  convert wetted extent poly to line
    arcpy.FeatureToLine_management(wePoly, 'in_memory/weLine', '', 'NO_ATTRIBUTES')

    #  split line using edge of water points
    arcpy.SplitLineAtPoint_management('in_memory/weLine', wePoints, 'in_memory/weLineSplit', '0.2 Meters')
    arcpy.MakeFeatureLayer_management('in_memory/weLineSplit', 'weLineSplit_lyr')

    #  make separate lyrs for rw + lw water points
    arcpy.MakeFeatureLayer_management(wePoints, 'rw_lyr', """ "DESCRIPTIO" = 'rw' """)
    arcpy.MakeFeatureLayer_management(wePoints, 'lw_lyr', """ "DESCRIPTIO" = 'lw' """)

    #  create separate rr/rw + rl/lw polylines using water point lyrs
    arcpy.SelectLayerByLocation_management('weLineSplit_lyr', 'WITHIN_A_DISTANCE', 'rw_lyr', '0.2 Meters', 'NEW_SELECTION')
    arcpy.SelectLayerByLocation_management('weLineSplit_lyr', 'WITHIN_A_DISTANCE', 'lw_lyr', '0.2 Meters', 'REMOVE_FROM_SELECTION')
    arcpy.UnsplitLine_management('weLineSplit_lyr', 'in_memory/weLine_rr')

    arcpy.SelectLayerByLocation_management('weLineSplit_lyr', 'WITHIN_A_DISTANCE', 'lw_lyr', '0.2 Meters', 'NEW_SELECTION')
    arcpy.SelectLayerByLocation_management('weLineSplit_lyr', 'WITHIN_A_DISTANCE', 'rw_lyr', '0.2 Meters','REMOVE_FROM_SELECTION')
    arcpy.UnsplitLine_management('weLineSplit_lyr', 'in_memory/weLine_rl')

    #  read in undercuts xlsx
    if len(arcpy.ListFiles('UndercutBank*')) > 0:

        inTblUC = arcpy.ListFiles('UndercutBank*')[0]
        arcpy.ExcelToTable_conversion(inTblUC, 'in_memory/tbl_uc')

        #  join undercuts to units
        arcpy.JoinField_management(units, 'Unit_Numbe', 'in_memory/tbl_uc', 'ChannelU_1', 'Bank')

        #  make separate lyrs for rr + lr undercuts
        arcpy.MakeFeatureLayer_management(units, 'rr_uc_lyr', """ "Bank" = 'Right' """)
        arcpy.MakeFeatureLayer_management(units, 'rl_uc_lyr', """ "Bank" = 'Left' """)

        #  clip undercuts to water extent
        arcpy.Intersect_analysis(['rr_uc_lyr', 'in_memory/weLine_rr'], 'in_memory/rr_uc', 'NO_FID', '', 'line')
        arcpy.Intersect_analysis(['rl_uc_lyr', 'in_memory/weLine_rl'], 'in_memory/rl_uc', 'NO_FID', '', 'line')

        #  merge rr + rl undercuts into single fc
        arcpy.Merge_management(['in_memory/rr_uc', 'in_memory/rl_uc'], 'in_memory/ucs')

        #  calculate euclidean distance
        ucDist = EucDistance('in_memory/ucs')
    else:
        ucDist = Con(lw >= 0, 1000)

    #  clip to wetted extent and save output
    ucDistClip = ExtractByMask(ucDist, units)
    ucDistClip.save('UndercutDist.tif')

    #  --- output table with values for each input raster --

    #  covert lwd count raster to points
    arcpy.RasterToPoint_conversion('lwDensity.tif', 'in_memory/ciInputs', 'VALUE')

    #  rename auto-assigned 'Value' field to 'lwdCount'
    arcpy.AlterField_management('in_memory/ciInputs', 'GRID_CODE', 'lwDensity')
    arcpy.CopyFeatures_management('in_memory/ciInputs', 'coverIndex_Inputs.shp')

    #  extract pool + undercut distance values to points
    ExtractMultiValuesToPoints('coverIndex_Inputs.shp', [['UndercutDist.tif', 'ucDist'], ['DeepPoolDist.tif', 'poolDist']], 'NONE')

    #  convert to numpy array and save as txt file
    nparr = arcpy.da.FeatureClassToNumPyArray('coverIndex_Inputs.shp', ['SHAPE@X', 'SHAPE@Y', 'lwDensity', 'ucDist', 'poolDist'])
    numpy.savetxt(os.path.join(arcpy.env.workspace, 'coverIndex_Inputs.csv'), nparr, fmt="%s", delimiter=",", header = str('x,y,lwdCount,ucDist,poolDist'))

if __name__ == '__main__':
    main()
