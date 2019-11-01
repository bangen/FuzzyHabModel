# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Name: Supporting Functions
# Purpose: A series of useful functions, placed in one spot so they're easier to bug fix
#
# Author: Sara Bangen
# Created on: 10 June 2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import os
import arcpy
import glob
import shutil
import fnmatch


def find_folder(folder_location, folder_name):
    """
    If the folder exists, returns it. Otherwise, raises an error
    :param folder_location: Where to look
    :param folder_name: The folder to look for
    :return: Path to folder
    """
    folders = os.listdir(folder_location)
    for folder in folders:
        if folder.endswith(folder_name):
            return os.path.join(folder_location, folder)
    return None


def make_folder(path_to_location, new_folder_name):
    """
    Makes a folder and returns the path to it
    :param path_to_location: Where we want to put the folder
    :param new_folder_name: What the folder will be called
    :return: String
    """
    newFolder = os.path.join(path_to_location, new_folder_name)
    if not os.path.exists(newFolder):
        os.makedirs(newFolder)

    return None


def find_available_num_prefix(folder_root):
    """
    Tells us the next number for a folder in the directory given
    :param folder_root: Where we want to look for a number
    :return: A string, containing a number
    """
    taken_nums = [fileName[0:2] for fileName in os.listdir(folder_root)]
    POSSIBLENUMS = range(1, 100)
    for i in POSSIBLENUMS:
        string_version = str(i)
        if i < 10:
            string_version = '0' + string_version
        if string_version not in taken_nums:
            return string_version
    arcpy.AddWarning("There were too many files at " + folder_root + " to have another folder that fits our naming convention")
    return "100"


def find_available_num_suffix(folder_root):
    """
    Tells us the next number for a folder in the directory given
    :param folder_root: Where we want to look for a number
    :return: A string, containing a number
    """
    taken_nums = [fileName[-2:] for fileName in os.listdir(folder_root)]
    POSSIBLENUMS = range(1, 100)
    for i in POSSIBLENUMS:
        string_version = str(i)
        if i < 10:
            string_version = '0' + string_version
        if string_version not in taken_nums:
            return string_version
    arcpy.AddWarning("There were too many files at " + folder_root + " to have another folder that fits our naming convention")
    return "100"


def find_relative_path(path, project_root):
    """
    Looks for the relative path from the project root to the item in the path
    :param path:
    :param project_root:
    :return:
    """
    relative_path = ''
    while path != os.path.dirname(path): # While there are still
        if path == project_root:
            return relative_path
        path, basename = os.path.split(path)

        relative_path = os.path.join(basename, relative_path)
    raise Exception("Could not find relative path")


def find_file(proj_path, file_pattern):

    search_path = os.path.join(proj_path, file_pattern)
    if len(glob.glob(search_path)) > 0:
        file_path = glob.glob(search_path)[0]
    else:
        file_path = None

    return file_path


def recursive_find_file(root_dir, patterns, match_full = False):
    """Search recursively for files matching a specified pattern.
    example for full match:
    hydro_files = recursive_find_file(os.path.join(from_visit, 'VisitFolders/HydroModel'), ("dem_grid_results.csv", "Velocity.Magnitude.tif", "Depth.tif"), match_full = True)
    example to find all files in directory:
    hydro_files = recursive_find_file(os.path.join(from_visit, 'VisitFolders/HydroModel/HydroModel_' + str(hydro_dir)), (""), match_full = False)
    """
    matches = []

    for root, dirs, files in os.walk(root_dir):
        if match_full is True:
            for pattern in patterns:
                for file in fnmatch.filter(files, pattern):
                    matches.append(os.path.join(root, file))
        else:
            for file in files:
                if file.startswith(patterns):
                    matches.append(os.path.join(root, file))

    return matches


def copy_files(file_list, destination_path):

    if os.path.exists(destination_path):
        for file in file_list:
            shutil.copy2(file, destination_path)
    else:
        print 'Cannot copy files because destination folder does not exist: ' + destination_path

    return None


def make_raster(in_shp, field_name, raster_template, out_filepath, dir_path):

    desc = arcpy.Describe(raster_template)
    arcpy.env.overwriteOutput = True
    arcpy.env.extent = desc.Extent
    arcpy.env.outputCoordinateSystem = desc.SpatialReference
    arcpy.env.cellSize = desc.meanCellWidth

    desc_shp = arcpy.Describe(in_shp)

    tmp_ras = os.path.join(dir_path, 'tmp_ras.tif')
    if desc_shp.shapeType == "Point":
        arcpy.PointToRaster_conversion(in_shp, field_name, tmp_ras, 'MEAN')
    if desc_shp.shapeType == "Polygon":
        arcpy.PolygonToRaster_conversion(in_shp, field_name, tmp_ras, 'MAXIMUM_AREA')
    ras_clip = arcpy.sa.ExtractByMask(tmp_ras, raster_template)
    # arcpy.Clip_management(tmp_ras, out_raster = out_filepath, in_template_dataset = raster_template)
    ras_clip.save(out_filepath)

    arcpy.Delete_management(tmp_ras)


def make_cover_input_points(lw_ras, uc_dist_ras, pool_dist_ras, dir_path):

    ci_inputs = os.path.join(dir_path, 'ci_inputs.shp')

    #  covert lwd count raster to points
    arcpy.RasterToPoint_conversion(lw_ras, ci_inputs, 'VALUE')

    #  extract pool + undercut distance values to points
    arcpy.sa.ExtractMultiValuesToPoints(ci_inputs, [[lw_ras, 'lwDensity'], [uc_dist_ras, 'ucDist'], [pool_dist_ras, 'poolDist']], 'NONE')

    #  round the values
    with arcpy.da.UpdateCursor(ci_inputs, ['lwDensity', 'ucDist', 'poolDist']) as cursor:
        for row in cursor:
            row[0] = round(row[0], 3)
            row[1] = round(row[1], 3)
            row[2] = round(row[2], 3)
            cursor.updateRow(row)

    # add point id field
    arcpy.AddField_management(ci_inputs, "PointID", "LONG")
    idx = 1
    with arcpy.da.UpdateCursor(ci_inputs, ["PointID"]) as cursor:
        for row in cursor:
            row[0] = idx
            idx += 1
            cursor.updateRow(row)

    return(ci_inputs)


def make_fhm_input_points(depth_ras, vel_ras, d50_ras, cover_ras, dir_path):

    fhm_inputs = os.path.join(dir_path, 'fhm_inputs.shp')

    #  covert depth raster to points
    arcpy.RasterToPoint_conversion(depth_ras, fhm_inputs, 'VALUE')

    #  extract pool + undercut distance values to points
    arcpy.sa.ExtractMultiValuesToPoints(fhm_inputs, [[depth_ras, 'Depth'], [vel_ras, 'Vel'], [d50_ras, 'D50'], [cover_ras, 'CoverIndex']], 'NONE')

    #  round the values
    with arcpy.da.UpdateCursor(fhm_inputs, ['Depth', 'Vel', 'D50', 'CoverIndex']) as cursor:
        for row in cursor:
            row[0] = round(row[0], 3)
            row[1] = round(row[1], 3)
            row[2] = round(row[2], 3)
            row[3] = round(row[3], 3)
            cursor.updateRow(row)

    # add point id field
    arcpy.AddField_management(fhm_inputs, "PointID", "LONG")
    idx = 1
    with arcpy.da.UpdateCursor(fhm_inputs, ["PointID"]) as cursor:
        for row in cursor:
            row[0] = idx
            idx += 1
            cursor.updateRow(row)

    return(fhm_inputs)


def clean_temp_dir(dir_path):

    temp_files = glob.glob(os.path.join(dir_path, '*'))
    for temp_file in temp_files:
        try:
            arcpy.Delete_management(temp_file)
        except:
            pass
    try:
        shutil.rmtree(dir_path)
    except:
        pass


# function to check if 'D50' field in channel units shapefile
def check_d50_field(shp):
    field_names = [f.name for f in arcpy.ListFields(shp)]
    if 'D50' not in field_names:
        return False
    else:
        return True


def check_raster_overlap(ras_a, ras_b):

    shp_a = arcpy.RasterDomain_3d(ras_a, 'in_memory/shp_a', 'POLYGON')
    shp_b = arcpy.RasterDomain_3d(ras_b, 'in_memory/shp_b', 'POLYGON')

    area_a = sum([row[0] for row in arcpy.da.SearchCursor(shp_a, ['SHAPE@AREA'])])

    shp_overlap = arcpy.Intersect_analysis([shp_a, shp_b], 'in_memory/shp_overlap')
    area_overlap = sum([row[0] for row in arcpy.da.SearchCursor(shp_a, ['SHAPE@AREA'])])

    percent_overlap = area_overlap / area_a * 100

    print percent_overlap
    if percent_overlap < 90.0:

        return True
    else:
        return False


def remove_na_values(shp, field_name):

    with arcpy.da.UpdateCursor(shp, field_name) as cursor:
        for row in cursor:
            if row[0] < 0:
                cursor.deleteRow()
