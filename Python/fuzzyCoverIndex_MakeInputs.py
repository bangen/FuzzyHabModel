# Script name: fuzzyCoverIndex_MakeInputs.py
#
# Created by: Sara Bangen (sara.bangen@gmail.com)
#
# Assumptions:
#   Folder path(s) contains:
#       - Water_Depth.tif
#       - Channel_Units.shp (with field: 'UnitNumber')
#       - WExtent.shp
#       - EdgeofWater_Points.shp (with field: 'Code' ['rw' + 'lw'])
#       - LargeWood*.xlsx (with fields: 'LargeWoodType' ['dry' + 'wet'], 'SumLWDCount', 'UnitNumber']
#       - UndercutBank*.xlsx (with fields: 'Bank' ['Right' + 'Left'], UnitNumber']
#
# Output:
#   3 fish cover rasters:
#       - distance to deep pools
#       - distance to undercut bank
#       - large wood density (pieces per 100 m2)

# -----------------------------------
# Set user-defined  input parameters
# -----------------------------------

visit_path = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\CBW05583-028079\2012\VISIT_1029"
bfw = 7.34

water_depth = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\CBW05583-028079\2012\VISIT_1029\Topo\Water_Depth.tif"
channel_units = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\CBW05583-028079\2012\VISIT_1029\Topo\Channel_Units.shp"
channel_units_csv = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\CBW05583-028079\2012\VISIT_1029\AuxData\ChannelUnit.csv"
water_extent = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\CBW05583-028079\2012\VISIT_1029\Topo\WExtent.shp"
water_edge_points = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\CBW05583-028079\2012\VISIT_1029\Topo\EdgeofWater_Points.shp"

# -----------------------------------
# Start of script

#  import required modules and extensions
import os
import arcpy
from arcpy.sa import *
from SupportingFunctions import find_file
from SupportingFunctions import make_folder
from SupportingFunctions import clean_temp_dir

arcpy.CheckOutExtension('Spatial')


def main(visit_path, bfw, water_depth, channel_units, channel_units_csv, water_extent, water_edge_points):

    bfw = float(bfw)

    print 'Creating fuzzy cover index input rasters for: ' + visit_path
    
    # set aux and inputs folder path
    aux_path = os.path.join(visit_path, 'AuxData')
    make_folder(visit_path, 'Habitat/CoverIndex/Inputs')
    proj_path = os.path.join(visit_path, 'Habitat/CoverIndex')
    inputs_path = os.path.join(visit_path, 'Habitat/CoverIndex/Inputs')

    # create temporary directory
    make_folder(proj_path, 'Temp')
    temp_dir = os.path.join(proj_path, 'Temp')

    # define workspace environment settings
    arcpy.env.workspace = temp_dir
    arcpy.env.scratchWorkspace = temp_dir
    arcpy.env.overwriteOutput = True

    #  import required rasters + shapefiles
    wd_ras = Raster(water_depth)
    units = arcpy.CopyFeatures_management(channel_units, 'in_memory/units')

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

    #  read in cu csv
    tbl_cus = os.path.join('in_memory', 'tbl_cus')
    arcpy.TableToTable_conversion(channel_units_csv, 'in_memory', 'tbl_cus')
    arcpy.JoinField_management(units, 'UnitNumber', tbl_cus, 'UnitNumber', ['Tier1', 'Tier2'])

    #  --- lwd density --
    lw_csv = find_file(aux_path, 'LargeWood*')

    if lw_csv is not None:

        #  read in wood table
        tbl_lw = os.path.join('in_memory', 'tbl_lw')
        arcpy.TableToTable_conversion(lw_csv, 'in_memory', 'tbl_lw')
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

            tbl_lw_sum = os.path.join('in_memory', 'tbl_lw_sum')
            arcpy.Statistics_analysis(tbl_lw, tbl_lw_sum, [['Count', 'SUM']], 'UnitNumber')

            # join wood count to units
            arcpy.JoinField_management(units, 'UnitNumber', tbl_lw_sum, 'UnitNumber', 'SUM_Count')

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
            with arcpy.da.UpdateCursor(tbl_lw, 'LargeWoodType') as cursor:
                for row in cursor:
                    if row[0] == 'Dry':
                        cursor.deleteRow()

            # join wood count to units
            arcpy.JoinField_management(units, 'UnitNumber', tbl_lw, 'UnitNumber', ['SumLWDCount'])

            #  add/calculate lwd density field
            arcpy.AddField_management(units, 'lwDen100m2', 'DOUBLE')
            fields = ['SHAPE@AREA', 'SumLWDCount', 'lwDen100m2']
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
    lw = os.path.join(inputs_path, 'LargeWoodDensity.tif')
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
        tbl_wd = os.path.join('in_memory', 'tbl_wd')
        ZonalStatisticsAsTable('units_lyr', 'UnitNumber', wd_ras, tbl_wd, 'DATA', 'MAXIMUM')

        #  join max value back to pool units
        arcpy.JoinField_management('units_lyr', 'UnitNumber', tbl_wd, 'UnitNumber', ['MAX'])

        #  remove pools units with depth < 80 cm
        if bfw <= 10.0:
            arcpy.SelectLayerByAttribute_management('units_lyr', 'REMOVE_FROM_SELECTION', """ "MAX" < 0.43 """)
        else:
            arcpy.SelectLayerByAttribute_management('units_lyr', 'REMOVE_FROM_SELECTION', """ "MAX" < 0.86 """)

        print 'Deep pool count: ' + str(int(arcpy.GetCount_management('units_lyr').getOutput(0)))

        if int(arcpy.GetCount_management('units_lyr').getOutput(0)) > 0:
            #  calculate euclidean distance
            pool_dist_raw = EucDistance('units_lyr')
        else:
            pool_dist_raw = Con(lw_clip >= 0, 1000)
            arcpy.CopyRaster_management(pool_dist_raw, os.path.join(temp_dir, 'lw_con.tif'))
    else:
        pool_dist_raw = Con(lw_clip >= 0, 1000)

    #  clip to wetted extent and save output
    pool_dist = os.path.join(inputs_path, 'DeepPoolDistance.tif')
    pool_dist_clip = ExtractByMask(pool_dist_raw, clip_shp)
    pool_dist_clip.save(pool_dist)

    pool_dist_min = arcpy.GetRasterProperties_management(pool_dist, "MINIMUM")
    pool_dist_max = arcpy.GetRasterProperties_management(pool_dist, "MAXIMUM")
    print 'Pool dist min: ' + str(round(float(pool_dist_min.getOutput(0)), 2))
    print 'Pool dist max: ' + str(round(float(pool_dist_max.getOutput(0)), 2))

    #  clear selection on units lyr
    arcpy.SelectLayerByAttribute_management('units_lyr', "CLEAR_SELECTION")

    #  --- undercuts --

    #  read in undercuts xlsx
    uc_csv = find_file(aux_path, 'UndercutBank*')

    if uc_csv is not None:

        #  convert wetted extent poly to line
        we_line = os.path.join(temp_dir, 'we_line.shp')
        arcpy.FeatureToLine_management(water_extent, we_line, '', 'NO_ATTRIBUTES')

        #  split line using edge of water points
        we_line_split = os.path.join(temp_dir, 'we_line_split.shp')
        arcpy.SplitLineAtPoint_management(we_line, water_edge_points, we_line_split, '0.2 Meters')
        arcpy.MakeFeatureLayer_management(we_line_split, 'we_line_split_lyr')

        #  make separate lyrs for rw + lw water points
        arcpy.MakeFeatureLayer_management(water_edge_points, 'rw_lyr', """ "Code" = 'rw' """)
        arcpy.MakeFeatureLayer_management(water_edge_points, 'lw_lyr', """ "Code" = 'lw' """)

        #  create separate rr/rw + rl/lw polylines using water point lyrs
        rw_line = os.path.join(temp_dir, 'rw_line.shp')
        arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'rw_lyr', '0.2 Meters',
                                               'NEW_SELECTION')
        arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'lw_lyr', '0.2 Meters',
                                               'REMOVE_FROM_SELECTION')
        arcpy.UnsplitLine_management('we_line_split_lyr', rw_line)

        lw_line = os.path.join(temp_dir, 'lw_line.shp')
        arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'lw_lyr', '0.2 Meters',
                                               'NEW_SELECTION')
        arcpy.SelectLayerByLocation_management('we_line_split_lyr', 'WITHIN_A_DISTANCE', 'rw_lyr', '0.2 Meters',
                                               'REMOVE_FROM_SELECTION')
        arcpy.UnsplitLine_management('we_line_split_lyr', lw_line)

        tbl_uc = os.path.join('in_memory', 'tbl_uc')
        arcpy.TableToTable_conversion(uc_csv, 'in_memory', 'tbl_uc')

        #  join undercuts to units
        arcpy.JoinField_management(units, 'UnitNumber', tbl_uc, 'UnitNumber', ['Bank'])

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
    else:
        uc_dist_raw = Con(lw_clip >= 0, 1000)

    # output some temporary data
    arcpy.CopyRaster_management(uc_dist_raw, os.path.join(temp_dir, 'uc_dist_raw.tif'))

    #  clip to wetted extent and save output
    uc_dist = os.path.join(inputs_path, 'UndercutDistance.tif')
    uc_dist_clip = ExtractByMask(uc_dist_raw, clip_shp)
    uc_dist_clip.save(uc_dist)

    uc_dist_min = arcpy.GetRasterProperties_management(uc_dist, "MINIMUM")
    uc_dist_max = arcpy.GetRasterProperties_management(uc_dist, "MAXIMUM")
    print 'Undercut dist min: ' + str(round(float(uc_dist_min.getOutput(0)), 2))
    print 'Undercut dist max: ' + str(round(float(uc_dist_max.getOutput(0)), 2))

    # remove temp files, variables, environment settings
    del visit_path
    del bfw

    clean_temp_dir(dir_path = temp_dir)
    arcpy.ResetEnvironments()
    arcpy.Delete_management('in_memory')


if __name__ == '__main__':
    main(visit_path, bfw, water_depth, channel_units, channel_units_csv, water_extent, water_edge_points)