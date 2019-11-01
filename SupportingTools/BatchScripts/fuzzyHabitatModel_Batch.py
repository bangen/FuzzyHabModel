#  define parent folder path and run folder name for directory search
root_path = r'C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data'
overwrite_output = False
overwrite_inputs = False
visit_csv = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_Projectwide\RunLists\visitList_15Oct2019_ReRun2.csv"


#  import required modules and extensionsp
import os
import sys
import imp
import arcpy
import csv
from collections import defaultdict
supporting_functions = imp.load_source('SupportingFunctions', "C:/etal/Shared/Projects/USA/CHaMP/HabitatSuitability/wrk_Code/06_FuzzyHabModel_v2/SupportingFunctions.py")

arcpy.CheckOutExtension('Spatial')
sys.path.append(r'C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Code\06_FuzzyHabModel_v2')
from fuzzyChinookSpawner import chinook_spawner_fis
from fuzzySteelheadSpawner import steelhead_spawner_fis
from fuzzyChinookJuvenile import chinook_juvenile_fis


def main(overwrite_inputs, overwrite_output):

    arcpy.env.overwriteOutput = True  # set to overwrite output

    # change directory to the parent folder path
    os.chdir(root_path)

    # read in csv parameters and convert to python dictionary
    visit_dict = defaultdict(dict)
    with open(visit_csv, "rb") as infile:
        reader = csv.reader(infile)
        headers = next(reader)[1:]
        for row in reader:
            visit_dict[row[0]] = {key: str(value) for key, value in zip(headers, row[1:])}

    visit_no = 1
    # run function for each huc8 folder
    for visit in visit_dict:

        print 'Working on visit ' + str(visit_no)

        basin_dir = visit_dict[visit]['WatershedName']
        site_dir = visit_dict[visit]['SiteName']
        year_dir = visit_dict[visit]['VisitYear']
        visit_dir = visit
        hydro_dir = visit_dict[visit]['HydroModel']

        visit_path = os.path.join(root_path, basin_dir, site_dir, year_dir, visit_dir)
        proj_path = os.path.join(visit_path, 'Habitat/FuzzyHM', hydro_dir)

        bfw = visit_dict[visit]['AverageBFWidth']

        depth_ras = os.path.join(proj_path, "Inputs/Hydro/Depth.tif")
        vel_ras = os.path.join(proj_path, "Inputs/Hydro/Velocity.Magnitude.tif")
        d50_ras = os.path.join(proj_path, "Inputs/GrainSize/D50.tif")
        cover_ras = os.path.join(proj_path, "Inputs/CoverIndex/CoverIndex.tif")

        inputs_list = [depth_ras, vel_ras, d50_ras, cover_ras]

        chinook_spawner_ras = os.path.join(proj_path, 'Output/Chinook/Spawner/FuzzyChinookSpawner_DVSC.tif')
        steelhead_spawner_ras = os.path.join(proj_path, 'Output/Steelhead/Spawner/FuzzySteelheadSpawner_DVSC.tif')
        chinook_juvenile_ras = os.path.join(proj_path, 'Output/Chinook/Juvenile/FuzzyChinookJuvenile_DVS.tif')

        outputs_list = [chinook_spawner_ras, steelhead_spawner_ras, chinook_juvenile_ras]

        if not all([os.path.isfile(file) for file in outputs_list]) or overwrite_output is True:

            if all([os.path.isfile(file) for file in inputs_list]):

                print visit_path

                # create temporary directory
                supporting_functions.make_folder(proj_path, 'Temp')
                temp_dir = os.path.join(proj_path, 'Temp')

                # define workspace environment settings
                arcpy.env.workspace = temp_dir
                arcpy.env.scratchWorkspace = temp_dir
                arcpy.env.overwriteOutput = True

                input_points = os.path.join(proj_path, 'Inputs/FHMInputs.shp')

                if not os.path.isfile(input_points) or overwrite_inputs is True:

                    tmp_points = supporting_functions.make_fhm_input_points(depth_ras, vel_ras, d50_ras, cover_ras, temp_dir)
                    arcpy.CopyFeatures_management(tmp_points, input_points)

                # remove na values
                supporting_functions.remove_na_values(input_points, 'D50')
                supporting_functions.remove_na_values(input_points, 'CoverIndex')

                if not os.path.isfile(chinook_spawner_ras) or overwrite_output is True:
                    try:
                        chinook_spawner_points = os.path.join(temp_dir, 'chinook_spawner.shp')
                        arcpy.CopyFeatures_management(input_points, chinook_spawner_points)
                        chinook_spawner_fis(chinook_spawner_points, bfw)
                        supporting_functions.make_raster(chinook_spawner_points, "FuzzyHSI", depth_ras, os.path.join(proj_path, 'Output/Chinook/Spawner/FuzzyChinookSpawner_DVSC.tif'), dir_path=temp_dir)
                    except:
                        print 'Ran into issue running chinook spawner model'
                        pass

                if not os.path.isfile(steelhead_spawner_ras) or overwrite_output is True:
                    try:
                        steelhead_spawner_points = os.path.join(temp_dir, 'steelhead_spawner.shp')
                        arcpy.CopyFeatures_management(input_points, steelhead_spawner_points)
                        steelhead_spawner_fis(steelhead_spawner_points, bfw)
                        supporting_functions.make_raster(steelhead_spawner_points, "FuzzyHSI", depth_ras,os.path.join(proj_path, 'Output/Steelhead/Spawner/FuzzySteelheadSpawner_DVSC.tif'), dir_path=temp_dir)
                    except:
                        print 'Ran into issue running steelhead spawner model'
                        pass

                if not os.path.isfile(chinook_juvenile_ras) or overwrite_output is True:
                    try:
                        chinook_juvenile_points = os.path.join(temp_dir, 'chinook_juvenile.shp')
                        arcpy.CopyFeatures_management(input_points, chinook_juvenile_points)
                        chinook_juvenile_fis(chinook_juvenile_points, bfw)
                        supporting_functions.make_raster(chinook_juvenile_points, "FuzzyHSI", depth_ras, os.path.join(proj_path, 'Output/Chinook/Juvenile/FuzzyChinookJuvenile_DVS.tif'), dir_path=temp_dir)
                    except:
                        print 'Ran into issue running chinook juvenile model'
                        pass

                    # remove temp files, variables, environment settings
                    supporting_functions.clean_temp_dir(dir_path=temp_dir)
                    arcpy.ResetEnvironments()
                    arcpy.Delete_management('in_memory')

            else:
                print 'Not all inputs exist.  Skipping: ' + visit_path

        else:
            print 'Model output already exists.  Skipping: ' + visit_path

        visit_no += 1

if __name__ == "__main__":
    main(overwrite_inputs, overwrite_output)
