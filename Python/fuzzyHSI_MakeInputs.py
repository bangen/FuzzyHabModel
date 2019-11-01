# Script name: CHaMP_fishCoverIndex_Rasters.py (v 1.0)
#
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
#   Point shapefile with values for each of the 3 fish cover rasters
#   CSV of values for each raster + [x,y] coords -- turned off for now

# -----------------------------------
# Set user-defined  input parameters
# -----------------------------------
#data = r"C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\RunLists\visitList_NK15_AllInputs.csv"
visit_path = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029"

cover_points = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029\Hydro\dem_grid_results.csv"
depth_ras = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029\Sims\FIS\Inputs\Channel_Units.shp"
velocity_ras = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029\Sims\FIS\Inputs\ChannelUnit.csv"
d50_ras = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029\Sims\FIS\Inputs\WaterExtent.shp"

# -----------------------------------
# Start of script

#  import required modules and extensions
import numpy
import os
import arcpy
from arcpy.sa import *
from shutil import rmtree
from SupportingFunctions import find_file
import glob

arcpy.CheckOutExtension('Spatial')

def main(visit_path, bfw, water_depth, channel_units, channel_units_csv, water_extent, water_edge_points):

    # visit_path = row['visit.dir']
    # bfw = row['AveBFW']

    print visit_path
    print 'visit bankfull width: ' + str(bfw)
    
    # set output folder path
    inputs_path = visit_path + '/Sims/FIS/Inputs'

    # create temporary directory
    temp_dir = os.path.join(visit_path, 'Temp')
    if os.path.exists(temp_dir):
        rmtree(temp_dir)
    os.mkdir(temp_dir)

    # define workspace environment settings
    arcpy.env.workspace = temp_dir
    arcpy.env.scratchWorkspace = temp_dir
    arcpy.env.overwriteOutput = True

    #  import required rasters + shapefiles
    wd_ras = Raster(water_depth)
    units = arcpy.CopyFeatures_management(channel_units, os.path.join(temp_dir, 'units.shp'))

    #  set raster environment settings
    desc = arcpy.Describe(wd_ras)
    arcpy.env.extent = desc.Extent
    arcpy.env.outputCoordinateSystem = desc.SpatialReference
    arcpy.env.cellSize = desc.meanCellWidth

    # create layer from channel units shp
    arcpy.MakeFeatureLayer_management(units, 'units_lyr')

    # dissolve units shp -- used for clipping
    clip_shp = os.path.join(temp_dir, 'clip_shp.shp')
    arcpy.Dissolve_management(units, clip_shp, multi_part = "SINGLE_PART")

    #  --- channel unit tier 1/tier 2 attributes ---
    # if tier1/tier2 attributes exist delete them
    # [some cu fcs don't have fields, some do but are not populated, some do and are populated]
    tfield = arcpy.ListFields(units, 'Tier1')
    if len(tfield) >= 1:
        arcpy.DeleteField_management(units, ['Tier1', 'Tier2'])

    #  read in cu xlsx
    tbl_cus = os.path.join(temp_dir, 'tbl_cus.dbf')
    arcpy.TableToTable_conversion(channel_units_csv, temp_dir, 'tbl_cus.dbf')
    arcpy.JoinField_management(units, 'Unit_Numbe', tbl_cus, 'ChannelU_1', ['Tier1', 'Tier2'])

    #  --- lwd density --
    lw_csv = find_file(inputs_path, 'LargeWood*')

    if lw_csv is not None:

        #  read in wood table
        tbl_lw = os.path.join(temp_dir, 'tbl_lw.dbf')
        arcpy.TableToTable_conversion(lw_csv, temp_dir, 'tbl_lw.dbf')
        visit_year = arcpy.da.SearchCursor(tbl_lw, ('VisitYear',)).next()[0]
        print 'visit year: ' + str(visit_year)

        if visit_year >= 2014:
            #  remove 'dry' rows
            with arcpy.da.UpdateCursor(tbl_lw, 'LargeWoodType') as cursor:
                for row in cursor:
                    if row[0] == 'Dry':
                        cursor.deleteRow()

            arcpy.AddField_management(tbl_lw, 'Count', 'SHORT')
            with arcpy.da.UpdateCursor(tbl_lw, 'Count') as cursor:
                for row in cursor:
                    row[0] = 1
                    cursor.updateRow(row)

            tbl_lw_sum = os.path.join(temp_dir, 'tbl_lw_sum.dbf')
            arcpy.Statistics_analysis(tbl_lw, tbl_lw_sum, [['Count', 'SUM']], 'ChannelU_1')

            # join wood count to units
            arcpy.JoinField_management(units, 'Unit_Numbe', tbl_lw_sum, 'ChannelU_1', 'SUM_Count')

            #  add/calculate lwd density field
            arcpy.AddField_management(units, 'lwDen100m2', 'DOUBLE')
            fields = ['SHAPE@AREA', 'SUM_Count', 'lwDen100m2']
            with arcpy.da.UpdateCursor(units, fields) as cursor:
                for row in cursor:
                    if row[1] > 0:
                        row[2] = float(row[1] / row[0]) * 100
                    else:
                        row[2] = 0.0
                    cursor.updateRow(row)
        else:
            #  remove 'dry' rows
            with arcpy.da.UpdateCursor(tbl_lw, 'LargeWoodT') as cursor:
                for row in cursor:
                    if row[0] == 'Dry':
                        cursor.deleteRow()

            # join wood count to units
            arcpy.JoinField_management(units, 'Unit_Numbe', tbl_lw, 'ChannelU_1', ['SumLWDCoun'])

            #  add/calculate lwd density field
            arcpy.AddField_management(units, 'lwDen100m2', 'DOUBLE')
            fields = ['SHAPE@AREA', 'SumLWDCoun', 'lwDen100m2']
            with arcpy.da.UpdateCursor(units, fields) as cursor:
                for row in cursor:
                    if row[1] > 0:
                        row[2] = float(row[1] / row[0]) * 100
                    else:
                        row[2] = 0.0
                    cursor.updateRow(row)

        #  convert to raster
        lw_raw = os.path.join(temp_dir, 'lw_raw.tif')
        arcpy.PolygonToRaster_conversion(units, 'lwDen100m2', lw_raw, 'CELL_CENTER', '', 0.1)

    else:
        lw_raw = Con(wd_ras >= 0, 0)

    # clip to wetted extent and save output
    lw = os.path.join(inputs_path, 'lwDensity.tif')
    lw_clip = ExtractByMask(lw_raw, clip_shp)
    lw_clip.save(lw)

    lw_min = arcpy.GetRasterProperties_management(lw, "MINIMUM")
    lw_max = arcpy.GetRasterProperties_management(lw, "MAXIMUM")
    print 'Large wood min: ' + str(round(float(lw_min.getOutput(0)), 2))
    print 'Large wood max: ' + str(round(float(lw_max.getOutput(0)), 2))

    #  --- distance to deep pool raster --

    #  select all pool units
    arcpy.SelectLayerByAttribute_management('units_lyr', 'NEW_SELECTION', """ "Tier1" = 'Slow/Pool' """)

    if int(arcpy.GetCount_management('units_lyr').getOutput(0)) > 0:
        #  get max water depth value for each pool unit
        tbl_wd = os.path.join(temp_dir, 'tbl_wd.dbf')
        ZonalStatisticsAsTable('units_lyr', 'Unit_Numbe', wd_ras, tbl_wd, 'DATA', 'MAXIMUM')

        #  join max value back to pool units
        arcpy.JoinField_management('units_lyr', 'Unit_Numbe', tbl_wd, 'Unit_Numbe', ['MAX'])

        #  remove pools units with depth < 80 cm
        if bfw <= 10.0:
            arcpy.SelectLayerByAttribute_management('units_lyr', 'REMOVE_FROM_SELECTION', """ "MAX" < 0.43 """)
        else:
            arcpy.SelectLayerByAttribute_management('units_lyr', 'REMOVE_FROM_SELECTION', """ "MAX" < 0.86 """)

        if int(arcpy.GetCount_management('units_lyr').getOutput(0)) > 0:
            #  calculate euclidean distance
            pool_dist_raw = EucDistance('units_lyr')
            pool_dist_norm = pool_dist_raw / bfw
        else:
            pool_dist_norm = Con(lw >= 0, 100)
    else:
        pool_dist_norm = Con(lw >= 0, 100)

    #  clip to wetted extent and save output
    pool_dist = os.path.join(inputs_path, 'DeepPoolDist.tif')
    pool_dist_clip = ExtractByMask(pool_dist_norm, clip_shp)
    pool_dist_clip.save(pool_dist)

    pool_dist_min = arcpy.GetRasterProperties_management(pool_dist, "MINIMUM")
    pool_dist_max = arcpy.GetRasterProperties_management(pool_dist, "MAXIMUM")
    print 'Pool dist min: ' + str(round(float(pool_dist_min.getOutput(0)), 2))
    print 'Pool dist max: ' + str(round(float(pool_dist_max.getOutput(0)), 2))

    #  clear selection on units lyr
    arcpy.SelectLayerByAttribute_management('units_lyr', "CLEAR_SELECTION")

    #  --- undercuts --

    #  convert wetted extent poly to line
    we_line = os.path.join(temp_dir, 'we_line.shp')
    arcpy.FeatureToLine_management(water_extent, we_line, '', 'NO_ATTRIBUTES')

    #  split line using edge of water points
    we_line_split = os.path.join(temp_dir, 'we_line_split.shp')
    arcpy.SplitLineAtPoint_management(we_line, water_edge_points, we_line_split, '0.2 Meters')
    arcpy.MakeFeatureLayer_management(we_line_split, 'we_line_split_lyr')

    #  make separate lyrs for rw + lw water points
    arcpy.MakeFeatureLayer_management(water_edge_points, 'rw_lyr', """ "DESCRIPTIO" = 'rw' """)
    arcpy.MakeFeatureLayer_management(water_edge_points, 'lw_lyr', """ "DESCRIPTIO" = 'lw' """)

    #  create separate rr/rw + rl/lw polylines using water point lyrs
    rw_line = os.path.join(temp_dir, 'rw_line.shp')
    arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'rw_lyr', '0.2 Meters', 'NEW_SELECTION')
    arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'lw_lyr', '0.2 Meters', 'REMOVE_FROM_SELECTION')
    arcpy.UnsplitLine_management('we_line_split_lyr', rw_line)

    lw_line = os.path.join(temp_dir, 'lw_line.shp')
    arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'lw_lyr', '0.2 Meters', 'NEW_SELECTION')
    arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'rw_lyr', '0.2 Meters','REMOVE_FROM_SELECTION')
    arcpy.UnsplitLine_management('we_line_split_lyr', lw_line)

    #  read in undercuts xlsx
    uc_csv = find_file(inputs_path, 'UndercutBank*')
    if uc_csv is not None:

        tbl_uc = os.path.join(temp_dir, 'tbl_uc.dbf')
        arcpy.TableToTable_conversion(uc_csv, temp_dir, 'tbl_uc.dbf')

        #  join undercuts to units
        arcpy.JoinField_management(units, 'Unit_Numbe', tbl_uc, 'ChannelU_1', ['Bank'])

        #  make separate lyrs for rw + lw undercuts
        arcpy.MakeFeatureLayer_management(units, 'rw_uc_lyr', """ "Bank" = 'Right' """)
        arcpy.MakeFeatureLayer_management(units, 'lw_uc_lyr', """ "Bank" = 'Left' """)

        #  clip undercuts to water extent
        rw_uc = os.path.join(temp_dir, 'rw_uc.shp')
        lw_uc = os.path.join(temp_dir, 'lw_uc.shp')
        ucs = os.path.join(temp_dir, 'ucs.shp')

        arcpy.Intersect_analysis(['rw_uc_lyr', rw_line], rw_uc, 'NO_FID', '', 'line')
        arcpy.Intersect_analysis(['lw_uc_lyr', lw_line], lw_uc, 'NO_FID', '', 'line')

        #  merge rw + rl undercuts into single fc
        arcpy.Merge_management([rw_uc, lw_uc], ucs)

        #  calculate euclidean distance
        uc_dist_raw = EucDistance(ucs)
        uc_dist_norm = uc_dist_raw / bfw
    else:
        uc_dist_norm = Con(lw >= 0, 100)

    #  clip to wetted extent and save output
    uc_dist = os.path.join(inputs_path, 'UndercutDist.tif')
    uc_dist_clip = ExtractByMask(uc_dist_norm, clip_shp)
    uc_dist_clip.save(uc_dist)

    uc_dist_min = arcpy.GetRasterProperties_management(uc_dist, "MINIMUM")
    uc_dist_max = arcpy.GetRasterProperties_management(uc_dist, "MAXIMUM")
    print 'Undercut dist min: ' + str(round(float(uc_dist_min.getOutput(0)), 2))
    print 'Undercut dist max: ' + str(round(float(uc_dist_max.getOutput(0)), 2))

    #  --- output table with values for each input raster --

    #  covert lwd count raster to points
    arcpy.RasterToPoint_conversion(lw, 'in_memory/ciInputs', 'VALUE')

    #  rename auto-assigned 'Value' field to 'lwdCount'
    ci_inputs = os.path.join(inputs_path, 'coverIndex_Inputs.shp')
    arcpy.AlterField_management('in_memory/ciInputs', 'GRID_CODE', 'lwDensity')
    arcpy.CopyFeatures_management('in_memory/ciInputs', ci_inputs)

    #  extract pool + undercut distance values to points
    ExtractMultiValuesToPoints(ci_inputs, [[uc_dist, 'ucDist'], [pool_dist, 'poolDist']], 'NONE')

    # add point id field
    arcpy.AddField_management(ci_inputs, "PointID", "LONG")
    idx = 1
    with arcpy.da.UpdateCursor(ci_inputs, ["PointID"]) as cursor:
        for row in cursor:
            row[0] = idx
            idx += 1
            cursor.updateRow(row)

    # #  convert to numpy array and save as txt file
    # nparr = arcpy.da.FeatureClassToNumPyArray('coverIndex_Inputs.shp', ['SHAPE@X', 'SHAPE@Y', 'lwDensity', 'ucDist', 'poolDist'])
    # numpy.savetxt(os.path.join(arcpy.env.workspace, 'coverIndex_Inputs.csv'), nparr, fmt="%s", delimiter=",", header = str('x,y,lwdCount,ucDist,poolDist'), comments = '')

    # remove temp files, variables, environment settings
    del visit_path
    del bfw
    arcpy.Delete_management('units_lyr')

    temp_files = glob.glob(os.path.join(temp_dir, '*'))
    for temp_file in temp_files:
        try:
            arcpy.Delete_management(temp_file)
        except:
            pass

    try:
        rmtree(temp_dir)
    except:
        pass

    arcpy.ResetEnvironments()


if __name__ == '__main__':
    main(visit_path, bfw, water_depth, channel_units, channel_units_csv, water_extent, water_edge_points)