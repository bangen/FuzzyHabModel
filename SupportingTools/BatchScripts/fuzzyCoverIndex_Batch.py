#  define parent folder path and run folder name for directory search
root_path = r'C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data'
overwrite_output = False
overwrite_inputs = False
visit_csv = r"C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Data\00_Projectwide\RunLists\visitList_15Oct2019_ReRun2.csv"


#  import required modules and extensionsp
import os
import sys
import arcpy
import csv
from collections import defaultdict
import imp
supporting_functions = imp.load_source('SupportingFunctions', "C:/etal/Shared/Projects/USA/CHaMP/HabitatSuitability/wrk_Code/06_FuzzyHabModel_v2/SupportingFunctions.py")


arcpy.CheckOutExtension('Spatial')
sys.path.append(r'C:\etal\Shared\Projects\USA\CHaMP\HabitatSuitability\wrk_Code\06_FuzzyHabModel_v2')
from fuzzyCoverIndex_MakeInputs import main as make_inputs
from fuzzyCoverIndex import main as run_fis


def main():

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

    # run function for each huc8 folder
    for visit in visit_dict:

        basin_dir = visit_dict[visit]['WatershedName']
        site_dir = visit_dict[visit]['SiteName']
        year_dir = visit_dict[visit]['VisitYear']
        visit_dir = visit

        visit_path = os.path.join(root_path, basin_dir, site_dir, year_dir, visit_dir)
        proj_path = os.path.join(visit_path, 'Habitat/CoverIndex')

        bfw = visit_dict[visit]['AverageBFWidth']
        water_depth = os.path.join(visit_path, 'Topo', 'Water_Depth.tif')
        channel_units = os.path.join(visit_path, 'Topo', 'Channel_Units.shp')
        channel_units_csv = os.path.join(visit_path, 'AuxData', 'ChannelUnit.csv')
        water_extent = os.path.join(visit_path, 'Topo', 'WExtent.shp')
        water_edge_points = os.path.join(visit_path, 'Topo', 'EdgeofWater_Points.shp')

        inputs_list = [water_depth, channel_units, channel_units_csv, water_extent, water_edge_points]

        lw_ras = os.path.join(proj_path, "Inputs", "LargeWoodDensity.tif")
        uc_dist_ras = os.path.join(proj_path, "Inputs", "UndercutDistance.tif")
        pool_dist_ras = os.path.join(proj_path, "Inputs", "DeepPoolDistance.tif")

        inputs_ras_list = [lw_ras, uc_dist_ras, pool_dist_ras]

        if not os.path.isfile(os.path.join(visit_path, 'Habitat/CoverIndex/Output/CoverIndex.tif')) or overwrite_output is True:

            if not all([os.path.isfile(file) for file in inputs_ras_list]) or overwrite_inputs is True:

                if all([os.path.isfile(file) for file in inputs_list]):
                    make_inputs(visit_path, bfw, water_depth, channel_units, channel_units_csv, water_extent, water_edge_points)
                else:
                    print 'Not all inputs exist.  Skipping: ' + visit_path

            run_fis(proj_path)

            supporting_functions.clean_temp_dir(os.path.join(proj_path, 'Temp'))

        else:
            print 'Model output already exists.  Skipping: ' + visit_path


if __name__ == "__main__":
    main()
