# -------------------------------------------------------------------------------
# Name:        Comb_FIS
# Purpose:     Runs the combined FIS for the BRAT input table
#
# Author:      Jordan Gilbert
#
# Created:     09/2016
# Copyright:   (c) Jordan 2016
# Licence:     <your licence>
# -------------------------------------------------------------------------------

in_points = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029\Sims\FIS\Inputs\coverIndex_Inputs.shp"
out_name = "FuzzyChinookSteelhead_DVSC.tif"
bfw = 7.34
raster_template = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_CodeTest\CBW05583-028079\VISIT_1029\Sims\FIS\Inputs\Water_Depth.tif"

import arcpy
import skfuzzy as fuzz
from skfuzzy import control as ctrl
from SupportingFunctions import make_folder
from SupportingFunctions import make_raster
import numpy as np
import os

arcpy.CheckOutExtension("Spatial")


def main(
    in_points,
    bfw,
    out_name,
    raster_template):

    # create temporary directory
    make_folder(proj_path, 'Temp')
    temp_dir = os.path.join(proj_path, 'Temp')

    # define workspace environment settings
    arcpy.env.workspace = temp_dir
    arcpy.env.scratchWorkspace = temp_dir
    arcpy.env.overwriteOutput = True

    # parent_folder = os.path.dirname(os.path.dirname(in_points))
    # output_folder = make_folder(parent_folder, "Output")
    # if out_name.endswith('.shp'):
    #     out_points = os.path.join(output_folder, out_name)
    # else:
    #     out_points = os.path.join(output_folder, out_name + ".shp")
    #
    # if os.path.exists(out_points):
    #     arcpy.Delete_management(out_points)
    # arcpy.CopyFeatures_management(in_points, out_points)

    # run the combined fis function for both potential and existing
    fhm_points = steelhead_spawner_fis(in_points, bfw)

    out_ras = os.path.join(output_folder, os.path.splitext(out_name)[0])

    make_raster(fhm_points, "FuzzyHSI", raster_template, out_ras)

    # remove temp files, variables, environment settings
    clean_temp_dir(dir_path=temp_dir)
    arcpy.ResetEnvironments()
    arcpy.Delete_management('in_memory')


# steelhead spawner fis function
def steelhead_spawner_fis(in_points, bfw):

    arcpy.env.overwriteOutput = True

    # get list of all fields in the input shapefile
    fields = [f.name for f in arcpy.ListFields(in_points)]

    # set output field name
    out_field = "FuzzyHSI"

    # check for cover index field in the network attribute table and delete if exists
    if out_field in fields:
        arcpy.DeleteField_management(in_points, out_field)

    if "pointid" not in fields:
        arcpy.AddField_management(in_points, "pointid", "LONG")
        idx = 1
        with arcpy.da.UpdateCursor(in_points, ["pointid"]) as cursor:
            for row in cursor:
                row[0] = idx
                idx += 1
                cursor.updateRow(row)

    # get arrays for fields of interest
    pointid_np = arcpy.da.FeatureClassToNumPyArray(in_points, "pointid")
    depth_np = arcpy.da.FeatureClassToNumPyArray(in_points, "Depth")
    velocity_np = arcpy.da.FeatureClassToNumPyArray(in_points, "Vel")
    d50_np = arcpy.da.FeatureClassToNumPyArray(in_points, "D50")
    cover_np = arcpy.da.FeatureClassToNumPyArray(in_points, "CoverIndex")

    pointid_array = np.asarray(pointid_np, np.int64)
    depth_array = np.asarray(depth_np, np.float64)
    velocity_array = np.asarray(velocity_np, np.float64)
    d50_array = np.asarray(d50_np, np.float64)
    cover_array = np.asarray(cover_np, np.float64)

    # check that inputs are within range of fis
    # if not, re-assign the value to just within range
    if bfw > 10:
        depth_array[depth_array > 4] = 4
    else:
        depth_array[depth_array > 2] = 2
    velocity_array[velocity_array > 4] = 4
    d50_array[d50_array > 4000] = 4000
    cover_array[cover_array > 1] = 1

    # delete temp arrays
    items = [pointid_np, depth_np, velocity_np, d50_np, cover_np]
    for item in items:
        del item

    # create antecedent (input) and consequent (output) objects to hold universe variables and membership functions
    if bfw > 10:
        depth = ctrl.Antecedent(np.arange(0, 4, 0.01), 'depth')
    else:
        depth = ctrl.Antecedent(np.arange(0, 2, 0.01), 'depth')
    velocity = ctrl.Antecedent(np.arange(0, 4, 0.01), 'velocity')
    d50 = ctrl.Antecedent(np.arange(0, 4000, 0.01), 'd50')
    cover = ctrl.Antecedent(np.arange(0, 1, 0.01), 'cover')
    suitability = ctrl.Consequent(np.arange(0, 1, 0.01), 'suitability')

    # build membership functions for each antecedent and consequent object
    if bfw > 10:
        depth['veryshallow'] = fuzz.trapmf(depth.universe, [0, 0, 0.09, 0.19])
        depth['shallow'] = fuzz.trimf(depth.universe, [0.09, 0.19, 0.52])
        depth['moderate'] = fuzz.trimf(depth.universe, [0.19, 0.52, 0.77])
        depth['deep'] = fuzz.trapmf(depth.universe, [0.52, 0.77, 4, 4])
    else:
        depth['veryshallow'] = fuzz.trapmf(depth.universe, [0, 0, 0.06, 0.07])
        depth['shallow'] = fuzz.trimf(depth.universe, [0.06, 0.07, 0.23])
        depth['moderate'] = fuzz.trimf(depth.universe, [0.07, 0.23, 0.34])
        depth['deep'] = fuzz.trapmf(depth.universe, [0.23, 0.34, 2, 2])

    velocity['veryslow'] = fuzz.trapmf(velocity.universe, [0, 0, 0.03, 0.12])
    velocity['slow'] = fuzz.trimf(velocity.universe, [0.03, 0.12, 0.5])
    velocity['moderate'] = fuzz.trimf(velocity.universe, [0.12, 0.5, 0.69])
    velocity['fast'] = fuzz.trapmf(velocity.universe, [0.5, 0.69, 4, 4])

    d50['finessand'] = fuzz.trapmf(d50.universe, [0, 0, 2, 4])
    d50['finegravel'] = fuzz.trapmf(d50.universe, [2, 4, 16, 32])
    d50['coarsegravel'] = fuzz.trapmf(d50.universe, [16, 32, 64, 96])
    d50['cobbleboulder'] = fuzz.trapmf(d50.universe, [64, 96, 4000, 4000])

    cover['none'] = fuzz.trimf(cover.universe, [0, 0, 0.1])
    cover['low'] = fuzz.trapmf(cover.universe, [0, 0.1, 0.3, 0.4])
    cover['moderate'] = fuzz.trapmf(cover.universe, [0.3, 0.4, 0.8, 0.9])
    cover['high'] = fuzz.trapmf(cover.universe, [0.8, 0.9, 1, 1])

    suitability['poor'] = fuzz.trapmf(suitability.universe, [0, 0, 0.1, 0.2])
    suitability['low'] = fuzz.trapmf(suitability.universe, [0.1, 0.2, 0.4, 0.5])
    suitability['moderate'] = fuzz.trapmf(suitability.universe, [0.4, 0.5, 0.8, 0.9])
    suitability['high'] = fuzz.trapmf(suitability.universe, [0.8, 0.9, 1, 1])

    # build fis rule table
    rule1 = ctrl.Rule(depth['veryshallow'] | velocity['veryslow'] | d50['finessand'], suitability['poor'])
    rule2 = ctrl.Rule(~depth['veryshallow'] & velocity['fast'] & ~d50['finessand'], suitability['low'])
    rule3 = ctrl.Rule(~depth['veryshallow'] & velocity['fast'] & ~d50['cobbleboulder'], suitability['low'])
    rule4 = ctrl.Rule(d50['cobbleboulder'], suitability['poor'])
    rule5 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['finegravel'] & cover['none'], suitability['poor'])
    rule6 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['coarsegravel'] & cover['none'], suitability['low'])
    rule7 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['finegravel'] & cover['none'], suitability['moderate'])
    rule8 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['coarsegravel'] & cover['none'], suitability['moderate'])
    rule9 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['finegravel'] & cover['none'], suitability['poor'])
    rule10 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['coarsegravel'] & cover['none'], suitability['low'])
    rule11 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['finegravel'] & cover['none'], suitability['moderate'])
    rule12 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['coarsegravel'] & cover['none'], suitability['moderate'])
    rule13 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['finegravel'] & cover['none'], suitability['poor'])
    rule14 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['coarsegravel'] & cover['none'], suitability['poor'])
    rule15 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['finegravel'] & cover['none'], suitability['low'])
    rule16 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['coarsegravel'] & cover['none'], suitability['moderate'])
    rule17 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['finegravel'] & cover['low'], suitability['poor'])
    rule18 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['coarsegravel'] & cover['low'], suitability['low'])
    rule19 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['finegravel'] & cover['low'], suitability['moderate'])
    rule20 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['coarsegravel'] & cover['low'], suitability['moderate'])
    rule21 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['finegravel'] & cover['low'], suitability['poor'])
    rule22 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['coarsegravel'] & cover['low'], suitability['low'])
    rule23 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['finegravel'] & cover['low'], suitability['moderate'])
    rule24 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['coarsegravel'] & cover['low'], suitability['moderate'])
    rule25 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['finegravel'] & cover['low'], suitability['poor'])
    rule26 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['coarsegravel'] & cover['low'], suitability['poor'])
    rule27 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['finegravel'] & cover['low'], suitability['low'])
    rule28 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['coarsegravel'] & cover['low'], suitability['moderate'])
    rule29 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['finegravel'] & cover['moderate'], suitability['low'])
    rule30 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['coarsegravel'] & cover['moderate'], suitability['moderate'])
    rule31 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['finegravel'] & cover['moderate'], suitability['high'])
    rule32 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['coarsegravel'] & cover['moderate'], suitability['high'])
    rule33 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['finegravel'] & cover['moderate'], suitability['low'])
    rule34 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['coarsegravel'] & cover['moderate'], suitability['moderate'])
    rule35 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['finegravel'] & cover['moderate'], suitability['high'])
    rule36 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['coarsegravel'] & cover['moderate'], suitability['high'])
    rule37 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['finegravel'] & cover['moderate'], suitability['low'])
    rule38 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['coarsegravel'] & cover['moderate'], suitability['low'])
    rule39 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['finegravel'] & cover['moderate'], suitability['moderate'])
    rule40 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['coarsegravel'] & cover['moderate'], suitability['high'])
    rule41 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['finegravel'] & cover['high'], suitability['moderate'])
    rule42 = ctrl.Rule(depth['shallow'] & velocity['slow'] & d50['coarsegravel'] & cover['high'], suitability['high'])
    rule43 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['finegravel'] & cover['high'], suitability['high'])
    rule44 = ctrl.Rule(depth['shallow'] & velocity['moderate'] & d50['coarsegravel'] & cover['high'], suitability['high'])
    rule45 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['finegravel'] & cover['high'], suitability['moderate'])
    rule46 = ctrl.Rule(depth['moderate'] & velocity['slow'] & d50['coarsegravel'] & cover['high'], suitability['high'])
    rule47 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['finegravel'] & cover['high'], suitability['high'])
    rule48 = ctrl.Rule(depth['moderate'] & velocity['moderate'] & d50['coarsegravel'] & cover['high'], suitability['high'])
    rule49 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['finegravel'] & cover['high'], suitability['moderate'])
    rule50 = ctrl.Rule(depth['deep'] & velocity['slow'] & d50['coarsegravel'] & cover['high'], suitability['moderate'])
    rule51 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['finegravel'] & cover['high'], suitability['high'])
    rule52 = ctrl.Rule(depth['deep'] & velocity['moderate'] & d50['coarsegravel'] & cover['high'], suitability['high'])

    ctrl_sys = ctrl.ControlSystem([rule1, rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9, rule10, rule11, rule12,
                                    rule13, rule14, rule15, rule16, rule17, rule18, rule19, rule20, rule21, rule22, rule23,
                                    rule24, rule25, rule26, rule27, rule28, rule29, rule30, rule31, rule32, rule33, rule34,
                                    rule35, rule36, rule37, rule38, rule39, rule40, rule41, rule42, rule43, rule44, rule45,
                                    rule46, rule47, rule48, rule49, rule50, rule51, rule52])

    ctrl_fis = ctrl.ControlSystemSimulation(ctrl_sys)

    # run fuzzy inference system on inputs and defuzzify output
    out = np.zeros(len(depth_array)) # todo: test this using nas instead of zeros
    for i in range(len(out)):
        ctrl_fis.input['depth'] = depth_array[i]
        ctrl_fis.input['velocity'] = velocity_array[i]
        ctrl_fis.input['d50'] = d50_array[i]
        ctrl_fis.input['cover'] = cover_array[i]
        ctrl_fis.compute()
        out[i] = ctrl_fis.output['suitability']

    # convert output value and point id to dictionary
    columns = zip(pointid_array, out)
    tblDict = {k: v for k, v in columns}

    # populate out field
    arcpy.AddField_management(in_points, out_field, 'DOUBLE')
    with arcpy.da.UpdateCursor(in_points, ["pointid", out_field]) as cursor:
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
    in_points,
    bfw,
    out_name,
    raster_template)