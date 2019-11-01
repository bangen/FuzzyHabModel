# -------------------------------------------------------------------------------
# Name:        fuzzyCoverIndex
# Purpose:     Runs the fuzzy cover index FIS model
#
# Author:      Sara Bangen
#
# Created:     06/2019
# -------------------------------------------------------------------------------

proj_path = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\cover_inputs\Lemhi\CBW05583-026031\2012\VISIT_916\Habitat\CoverIndex"

import arcpy
import skfuzzy as fuzz
from skfuzzy import control as ctrl
from SupportingFunctions import make_folder
from SupportingFunctions import make_raster
from SupportingFunctions import make_cover_input_points
from SupportingFunctions import clean_temp_dir
import numpy as np
import os

arcpy.CheckOutExtension('Spatial')


def main(proj_path):

    # set folder paths
    inputs_folder = os.path.join(proj_path, "Inputs")
    output_folder = os.path.join(proj_path, "Output")
    make_folder(proj_path, "Output")
    out_ras = os.path.join(output_folder, "CoverIndex.tif")

    # create temporary directory
    make_folder(proj_path, 'Temp')
    temp_dir = os.path.join(proj_path, 'Temp')

    # define workspace environment settings
    arcpy.env.workspace = temp_dir
    arcpy.env.scratchWorkspace = temp_dir
    arcpy.env.overwriteOutput = True

    lw_ras = os.path.join(inputs_folder, "LargeWoodDensity.tif")
    uc_dist_ras = os.path.join(inputs_folder, "UndercutDistance.tif")
    pool_dist_ras = os.path.join(inputs_folder, "DeepPoolDistance.tif")

    inputs_list = [lw_ras, uc_dist_ras, pool_dist_ras]

    if all([os.path.isfile(file) for file in inputs_list]):

        in_points = os.path.join(inputs_folder, "CoverIndexInputs.shp")
        if not os.path.exists(in_points):
            tmp_points = make_cover_input_points(lw_ras, uc_dist_ras, pool_dist_ras, dir_path = temp_dir)
            arcpy.CopyFeatures_management(tmp_points, in_points)

        # run the cover index fis function for both potential and existing
        print 'Running fuzzy cover index model for: ' + proj_path
        ci_points = cover_index_fis(in_points)

        print 'Creating cover index raster...: '
        make_raster(ci_points, "CoverIndex", lw_ras, out_ras, dir_path = temp_dir)

        # remove temp files, variables, environment settings
        clean_temp_dir(dir_path = temp_dir)
        arcpy.ResetEnvironments()
        arcpy.Delete_management('in_memory')

    else:
        print 'Not all inputs exist.  Skipping: ' + proj_path


# cover index fis function
def cover_index_fis(in_points):

    import time
    start_time = time.time()

    # get list of all fields in the input shapefile
    fields = [f.name for f in arcpy.ListFields(in_points)]

    # set output field name
    out_field = "CoverIndex"

    # check for cover index field in the network attribute table and delete if exists
    if out_field in fields:
        arcpy.DeleteField_management(in_points, out_field)

    if "PointID" not in fields:
        arcpy.AddField_management(in_points, "PointID", "LONG")
        idx = 1
        with arcpy.da.UpdateCursor(in_points, ["PointID"]) as cursor:
            for row in cursor:
                row[0] = idx
                idx += 1
                cursor.updateRow(row)

    # get arrays for fields of interest
    pointid_np = arcpy.da.FeatureClassToNumPyArray(in_points, "PointID")
    lwd_np = arcpy.da.FeatureClassToNumPyArray(in_points, "lwDensity")
    ucdist_np = arcpy.da.FeatureClassToNumPyArray(in_points, "ucDist")
    pooldist_np = arcpy.da.FeatureClassToNumPyArray(in_points, "poolDist")

    pointid_array = np.asarray(pointid_np, np.int64)
    lwd_array = np.asarray(lwd_np, np.float64)
    ucdist_array = np.asarray(ucdist_np, np.float64)
    pooldist_array = np.asarray(pooldist_np, np.float64)

    # check that inputs are within range of fis
    # if not, re-assign the value to just within range
    lwd_array[lwd_array < 0] = 0
    lwd_array[lwd_array > 400] = 400
    ucdist_array[ucdist_array < 0] = 0
    ucdist_array[ucdist_array > 100] = 100
    pooldist_array[pooldist_array < 0] = 0
    pooldist_array[pooldist_array > 100] = 100

    # delete temp arrays
    items = [pointid_np, lwd_np, ucdist_np, pooldist_np]
    for item in items:
        del item

    # create antecedent (input) and consequent (output) objects to hold universe variables and membership functions
    lwd = ctrl.Antecedent(np.arange(0, 400, 0.01), 'input1')
    ucdist = ctrl.Antecedent(np.arange(0, 100, 0.01), 'input2')
    pooldist = ctrl.Antecedent(np.arange(0, 100, 0.01), 'input3')
    coverindex = ctrl.Consequent(np.arange(0, 1, 0.01), 'result')

    # build membership functions for each antecedent and consequent object
    lwd['low'] = fuzz.trimf(lwd.universe, [0, 0, 1])
    lwd['moderate'] = fuzz.trapmf(lwd.universe, [0, 1, 3, 4])
    lwd['high'] = fuzz.trapmf(lwd.universe, [3, 4, 400, 400])

    ucdist['near'] = fuzz.trapmf(ucdist.universe, [0, 0, 4, 7])
    ucdist['notfar'] = fuzz.trimf(ucdist.universe, [4, 7, 10])
    ucdist['far'] = fuzz.trapmf(ucdist.universe, [7, 10, 100, 100])

    pooldist['near'] = fuzz.trapmf(pooldist.universe, [0, 0, 4, 7])
    pooldist['notfar'] = fuzz.trimf(pooldist.universe, [4, 7, 10])
    pooldist['far'] = fuzz.trapmf(pooldist.universe, [7, 10, 100, 100])

    coverindex['none'] = fuzz.trimf(coverindex.universe, [0, 0, 0.1])
    coverindex['low'] = fuzz.trapmf(coverindex.universe, [0, 0.1, 0.3, 0.4])
    coverindex['moderate'] = fuzz.trapmf(coverindex.universe, [0.3, 0.4, 0.8, 0.9])
    coverindex['high'] = fuzz.trapmf(coverindex.universe, [0.8, 0.9, 1, 1])

    # build fis rule table
    rule1 = ctrl.Rule(lwd['high'], coverindex['high'])
    rule2 = ctrl.Rule(lwd['low'] & ucdist['near'] & pooldist['near'], coverindex['high'])
    rule3 = ctrl.Rule(lwd['low'] & ucdist['near'] & pooldist['notfar'], coverindex['moderate'])
    rule4 = ctrl.Rule(lwd['low'] & ucdist['near'] & pooldist['far'], coverindex['moderate'])
    rule5 = ctrl.Rule(lwd['low'] & ucdist['notfar'] & pooldist['near'], coverindex['moderate'])
    rule6 = ctrl.Rule(lwd['low'] & ucdist['notfar'] & pooldist['notfar'], coverindex['low'])
    rule7 = ctrl.Rule(lwd['low'] & ucdist['notfar'] & pooldist['far'], coverindex['low'])
    rule8 = ctrl.Rule(lwd['low'] & ucdist['far'] & pooldist['near'], coverindex['moderate'])
    rule9 = ctrl.Rule(lwd['low'] & ucdist['far'] & pooldist['notfar'], coverindex['low'])
    rule10 = ctrl.Rule(lwd['low'] & ucdist['far'] & pooldist['far'], coverindex['none'])
    rule11 = ctrl.Rule(lwd['moderate'] & ucdist['near'] & pooldist['near'], coverindex['high'])
    rule12 = ctrl.Rule(lwd['moderate'] & ucdist['near'] & pooldist['notfar'], coverindex['high'])
    rule13 = ctrl.Rule(lwd['moderate'] & ucdist['near'] & pooldist['far'], coverindex['moderate'])
    rule14 = ctrl.Rule(lwd['moderate'] & ucdist['notfar'] & pooldist['near'], coverindex['high'])
    rule15 = ctrl.Rule(lwd['moderate'] & ucdist['notfar'] & pooldist['notfar'], coverindex['moderate'])
    rule16 = ctrl.Rule(lwd['moderate'] & ucdist['notfar'] & pooldist['far'], coverindex['low'])
    rule17 = ctrl.Rule(lwd['moderate'] & ucdist['far'] & pooldist['near'], coverindex['moderate'])
    rule18 = ctrl.Rule(lwd['moderate'] & ucdist['far'] & pooldist['notfar'], coverindex['low'])
    rule19 = ctrl.Rule(lwd['moderate'] & ucdist['far'] & pooldist['far'], coverindex['low'])

    ctrl_sys = ctrl.ControlSystem([rule1, rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9, rule10,
                                          rule11, rule12, rule13, rule14, rule15, rule16, rule17, rule18, rule19])

    ctrl_fis = ctrl.ControlSystemSimulation(ctrl_sys)

    # run fuzzy inference system on inputs and defuzzify output
    out = np.zeros(len(lwd_array)) # todo: test this using nas instead of zeros
    for i in range(len(out)):
        ctrl_fis.input['input1'] = round(lwd_array[i], 3)
        ctrl_fis.input['input2'] = round(ucdist_array[i], 3)
        ctrl_fis.input['input3'] = round(pooldist_array[i], 3)
        ctrl_fis.compute()
        out[i] = round(ctrl_fis.output['result'], 3)

    # convert output value and point id to dictionary
    columns = zip(pointid_array, out)
    tblDict = {k: v for k, v in columns}

    # populate out field
    arcpy.AddField_management(in_points, out_field, 'DOUBLE')
    with arcpy.da.UpdateCursor(in_points, ["PointID", out_field]) as cursor:
        for row in cursor:
            aKey = row[0]
            row[1] = tblDict[aKey]
            cursor.updateRow(row)
    tblDict.clear()

    # # todo: turning off this functionality for now to be consistent with how model was run in the past
    # #       but would be a good idea to turn this on at some point
    # # calculate defuzzified centroid value for output 'none' MF group
    # # this will be used to re-classify output values that fall in this group
    # # important: will need to update the array (x) and MF values (mfx) if the
    # #            density 'none' values are changed in the model
    # x = np.arange(0, 1, 0.001)
    # mfx = fuzz.trimf(x, [0, 0, 0.1])
    # defuzz_centroid = round(fuzz.defuzz(x, mfx, 'centroid'), 6)
    #
    #
    # # set output cover index to 0 if output falls fully in 'none' category
    # with arcpy.da.UpdateCursor(in_points, [out_field]) as cursor:
    #     for row in cursor:
    #         if round(row[0], 6) == defuzz_centroid:
    #             row[0] = 0.0
    #         cursor.updateRow(row)

    # delete temporary tables and arrays
    items = [columns, out]
    # items = [columns, out, x, mfx, defuzz_centroid]
    for item in items:
        del item

    return in_points


if __name__ == '__main__':
    main(
    proj_path)
