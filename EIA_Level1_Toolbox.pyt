"""
---------------------------------------------------------------------------------------------------------------------
Name: EIA Level 1 Toolbox
Purpose: This Python Toolbox allows users to calculate level 1 EIA metrics on selected polygons. This tool modifies the
original dataset.
Author: Molly Moore for Pennsylvania Natural Heritage Program
Created: 12/2/2024
Updates:
------------------------------------------------------------------------------------------------------------------------
"""

########################################################################################################################
## Import packages and define environment settings
########################################################################################################################

import arcpy
import os
import sys
import string
import pandas as pd
import datetime
from arcgis.gis import GIS
from arcgis.features import FeatureLayer
from arcgis.raster import ImageryLayer

arcpy.env.overwriteOutput = True
arcpy.env.transferDomains = True


########################################################################################################################
## Begin toolbox
########################################################################################################################

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "EIA Level 1 Toolbox"
        self.alias = "EIA Level 1 Toolbox"
        self.canRunInBackground = False
        self.tools = [CalculateLevel1Metrics]


########################################################################################################################
## Begin Calculate Level 1 EIA Metrics tool
########################################################################################################################

class CalculateLevel1Metrics(object):
    def __init__(self):
        self.label = "Calculate Level 1 EIA Metrics"
        self.description = "This tool modifies the input dataset"
        self.canRunInBackground = False

    def getParameterInfo(self):
        selected_input = arcpy.Parameter(
            displayName="Selected Communities (must be polygon layer)",
            name="selected_input",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        nlcd = arcpy.Parameter(
            displayName="NLCD Raster Layer (will perform faster if using NLCD raster already projected into custom Albers)",
            name="nlcd",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        nlcd.value = r"W:\\Heritage\\Heritage_Data\\Heritage_Data_Tools\\EIA_Level1_Tools\\EIA_Data.gdb\\NLCD_2019_pa_albers"

        scratch_workspace = arcpy.Parameter(
            displayName = "Scratch Geodatabase (geodatabase where intermediate geoprocessing files will be written)",
            name = "scratch_workspace",
            datatype = "DEWorkspace",
            parameterType = "Required",
            direction = "Input")

        params = [selected_input, nlcd, scratch_workspace]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, params):
        return

    def updateMessages(self, params):
        return

    def execute(self, params, messages):

        selected_input = params[0].valueAsText
        nlcd = params[1].valueAsText
        scratch_workspace = params[2].valueAsText

        arcpy.AddMessage("Preparing data for EIA Analysis.")
        arcpy.AddMessage("................................")

        # define rest endpoint path to nlcd lookup table that includes LAN2 coefficients and natural/non-natural land
        # cover designations
        lu_nlcd = r"https://gis.waterlandlife.org/server/rest/services/Hosted/lu_NLCD_EIA/FeatureServer/0"
        # define scratch workspace
        scratch_gdb = scratch_workspace

        # check to see if there is a selection on polygon input layer. If there isn't, error out. We want to include
        # this because FIND community layer will be used.
        desc = arcpy.Describe(selected_input)
        if not desc.FIDSet == '':
            pass
        else:
            arcpy.AddError("No community features are selected. Please make a selection and try again.")
            sys.exit()

        arcpy.AddMessage("Checking projections.")

        # define albers projection that we will use for all calculations
        albers_str = r'PROJCS["alber",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-78.0],PARAMETER["standard_parallel_1",40.0],PARAMETER["standard_parallel_2",42.0],PARAMETER["latitude_of_origin",39.0],UNIT["Meter",1.0]];-16085300 -8515400 279982320.962027;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision'
        albers_prj = arcpy.SpatialReference()
        albers_prj.loadFromString(albers_str)

        # check for transformations between input data and albers projection. if any exist, set them for output transformations
        input_sr = arcpy.Describe(selected_input).spatialReference
        transformations = arcpy.ListTransformations(input_sr, albers_prj)
        if len(transformations) == 0:
            transformation = ""
        else:
            transformation = transformations[0]

        # set output environment to albers projection and transformation if it exists
        arcpy.env.outputCoordinateSystem = albers_prj
        arcpy.env.geographicTransformations = transformation

        # get coordinate system for NLCD raster
        nlcd_sr = arcpy.Describe(nlcd).spatialReference
        nlcd_sr_string = nlcd_sr.exportToString()

        # if coordinate system for nlcd raster is not albers, project into albers
        if nlcd_sr_string != albers_str:
            arcpy.AddMessage("- Input NLCD raster is being projected into custom Albers projection.")
            nlcd = arcpy.ProjectRaster_management(nlcd, os.path.join(scratch_gdb,"nlcd_projected"), albers_prj)
        else:
            arcpy.AddMessage("- Projections look good. Moving on.")
            arcpy.AddMessage("................................")

        arcpy.AddMessage("Checking for EIA fields in input layer. If field is not found, it will be added.")

        # define list of fields to hold metrics we will be calculating
        eia_fields = ["lan1", "lan1_score", "lan2", "lan2_score", "buf1", "buf1_score", "buf2", "buf2_score"]
        eia_aliases = ["LAN1 Contiguous Natural Land Cover", "LAN1 Score", "LAN2 Land Use Index", "LAN2 Score",
                       "BUF1 Perimeter With Natural Edge", "BUF1 Score", "BUF2 Widge of Natural Edge", "BUF2 Score"]

        # get list of current fields in input feature class
        current_fields = [f.name for f in arcpy.ListFields(selected_input)]
        # loop through fields, if they exist, skip with message, if they don't exist, add the field with appropriate length and alias
        for field, alias in zip(eia_fields, eia_aliases):
            if field in current_fields:
                arcpy.AddMessage("- Field " + field + " already exists in input dataset.")
                pass
            else:
                arcpy.AddMessage("- Adding field: " + field)
                if field.endswith("_score"):
                    field_length = 5
                else:
                    field_length = 255
                try:
                    arcpy.AddField_management(selected_input, field, "TEXT", field_length=field_length, field_alias=alias)
                except arcpy.ExecuteError:
                    error_code = arcpy.GetMessages(2).split(":")[0]
                    if error_code == "ERROR 000499":
                        arcpy.AddError("Looks like the input layer does not have the appropriate EIA metric fields and that you do not have permissions to add fields to the dataset. Try exporting the features to a file geodatabase and try again with the exported features.")
                        sys.exit()
                    else:
                        arcpy.AddError(arcpy.GetMessages(2))

        arcpy.AddMessage("................................")

        # create list of object ids to loop through features and calculate EIA metrics
        oid_fieldname = arcpy.Describe(selected_input).OIDFieldName
        input_oids = sorted({row[0] for row in arcpy.da.SearchCursor(selected_input, oid_fieldname)})

        # begin loop to calculate metrics
        for oid in input_oids:
            arcpy.AddMessage("Calculating EIA Metrics for OID: "+str(oid))
            # make layer of input feature, filtered by OID
            input_lyr = arcpy.ExportFeatures_conversion(selected_input, os.path.join(scratch_gdb,"input_lyr_projected"), where_clause="{0} = {1}".format(oid_fieldname, oid))
            # create buffers at 100 and 500 meters around input feature for analysis buffer
            multi_buff = arcpy.MultipleRingBuffer_analysis(input_lyr, os.path.join(scratch_gdb, "multi_buff"),
                                                           [100, 500], "Meters", "distance", "ALL", "OUTSIDE_ONLY",
                                                           "PLANAR")
            # clip raster to buffer area to make it more manageable and quicker for analysis
            clip_raster = arcpy.Clip_management(in_raster=nlcd, out_raster=os.path.join(scratch_gdb, "clip_raster"),
                                                in_template_dataset=multi_buff)

            # convert clipped raster to polygon for line intersects
            nlcd_poly = arcpy.RasterToPolygon_conversion(clip_raster, os.path.join(scratch_gdb, "nlcd_poly"),
                                                         raster_field="Value", simplify="NO_SIMPLIFY")
            # join field for natural/non-natural designation
            arcpy.JoinField_management(nlcd_poly, "gridcode", lu_nlcd, "Value", "cover_type_value")

            # dissolve nlcd polygon layer to natural/non-natural categories
            natpoly = arcpy.analysis.PairwiseDissolve(nlcd_poly, os.path.join(scratch_gdb,"natpoly"), "cover_type_value")
            with arcpy.da.UpdateCursor(natpoly, "cover_type_value") as cursor:
                for row in cursor:
                    if row[0] == 0:
                        cursor.deleteRow()

########################################################################################################################
# begin section to calculate lan1 and lan2 metrics
########################################################################################################################

            # tabulate intersection of buffer layer with nlcd polygon for lan1 and lan2 calculations
            buff_intersect = arcpy.analysis.TabulateIntersection(multi_buff, "distance", nlcd_poly,
                                                                 os.path.join(scratch_gdb, "tab_intersect"), "gridcode")

            # convert tab intersect to pandas dataframe to do calculations
            final_fields = [f.name for f in arcpy.ListFields(buff_intersect)]
            tab_area_df = pd.DataFrame((row for row in arcpy.da.SearchCursor(buff_intersect, final_fields)),
                                       columns=final_fields)

            # load lu_nlcd table from rest endpoint
            lu_nlcd_flayer = FeatureLayer(lu_nlcd)
            lu_nlcd_query = lu_nlcd_flayer.query()
            lu_nlcd_df = lu_nlcd_query.sdf

            # join nlcd coefficient and natural/non-natural distinctions from lu nlcd table
            tab_area_df = pd.merge(tab_area_df, lu_nlcd_df, left_on="gridcode", right_on="value", how="left")

            # calculate percent area for each land cover class within polygon buffers
            tab_area_df['percent'] = tab_area_df['AREA'] / sum(tab_area_df['AREA']) * 100
            tab_area_df['coefficientXpercent'] = tab_area_df['coefficient'] * (tab_area_df['PERCENTAGE'] / 100)

            # calculate lan1 by adding percent of natural land cover classes
            lan1 = round(tab_area_df.loc[tab_area_df['cover_type_value'] == 1, 'percent'].sum(), 3)

            if lan1 >= 90:
                lan1_rank = "A"
            elif 60 <= lan1 <= 90:
                lan1_rank = "B"
            elif 20 <= lan1 <= 60:
                lan1_rank = "C"
            else:
                lan1_rank = "D"

            arcpy.AddMessage("- LAN1 value: " + str(lan1) + " LAN1 score: " + lan1_rank)

            # calculate lan2
            lan2 = round((tab_area_df.loc[tab_area_df['distance'] == 100, 'coefficientXpercent'].sum() * 0.6) + (
                        tab_area_df.loc[tab_area_df['distance'] == 500, 'coefficientXpercent'].sum() * 0.4), 3)

            if lan2 >= 9.5:
                lan2_rank = "A"
            elif 8.0 <= lan2 <= 9.5:
                lan2_rank = "B"
            elif 4.0 <= lan2 <= 8.0:
                lan2_rank = "C"
            else:
                lan2_rank = "D"

            arcpy.AddMessage("- LAN2 value: " + str(lan2) + " LAN2 score: " + lan2_rank)

########################################################################################################################
# begin section to calculate buf1
########################################################################################################################

            # convert feature to line
            input_line = arcpy.management.PolygonToLine(input_lyr, os.path.join(scratch_gdb, "input_line"),
                                                        "IGNORE_NEIGHBORS")

            # intersect input line with nlcd polygon to get length of line overlapping natural land cover
            line_intersect = arcpy.analysis.PairwiseIntersect("{0};{1}".format(input_line, nlcd_poly),
                                                              os.path.join(scratch_gdb, "line_intersect"), "ALL",
                                                              output_type="LINE")

            # convert line intersect to pandas dataframe to do calculations
            final_fields = ["gridcode", "Shape_Length"]
            line_intersect_df = pd.DataFrame((row for row in arcpy.da.SearchCursor(line_intersect, final_fields)),
                                             columns=final_fields)
            # join nlcd coefficient and natural/non-natural distinctions from lu nlcd table
            line_intersect_df = pd.merge(line_intersect_df, lu_nlcd_df, left_on="gridcode", right_on="value",
                                         how="left")

            # get percentage of line for each record
            line_intersect_df["percent"] = (line_intersect_df['Shape_Length'] / line_intersect_df[
                'Shape_Length'].sum()) * 100
            # calculate buf1 by adding percent of natural land cover classes overlapping perimeter of area of interest
            buf1 = round(line_intersect_df.loc[line_intersect_df['cover_type_value'] == 1, 'percent'].sum(), 3)

            if buf1 == 100:
                buf1_rank = "A"
            elif 75 <= buf1 < 100:
                buf1_rank = "B"
            elif 25 <= buf1 < 75:
                buf1_rank = "C"
            else:
                buf1_rank = "D"

            arcpy.AddMessage("- BUF1 value: " + str(buf1) + " BUF1 score: " + buf1_rank)

########################################################################################################################
# begin section to calculate buf2
########################################################################################################################

            # get length of area of interest line divided by 8 to segment into 8 sections between transects
            with arcpy.da.SearchCursor(input_line, "Shape_Length") as cursor:
                for row in cursor:
                    section_length = (row[0] / 8) - 1

            # generate 200m transects at equal intervals along area of interest line feature (so that there is 100m line extending to 100m buffer
            transects = arcpy.GenerateTransectsAlongLines_management(input_line, os.path.join(scratch_gdb, "transects"),
                                                                     section_length, "200 Meters", "NO_END_POINTS")

            # intersect input line with buffer polygon to get to keep 100m lines extending from area of interest
            transect_intersect = arcpy.analysis.PairwiseIntersect("{0};{1}".format(transects, multi_buff),
                                                                  os.path.join(scratch_gdb, "transect_intersect"),
                                                                  "ALL", output_type="LINE")

            # convert output to singlepart polygons
            transect_intersect = arcpy.MultipartToSinglepart_management(transect_intersect, os.path.join(scratch_gdb, "transects_singlepart"))

            # delete all transect lines that are under 99m to clean up overlapping transects
            with arcpy.da.UpdateCursor(transect_intersect,"Shape_Length") as cursor:
                for row in cursor:
                    if row[0] < 99:
                        cursor.deleteRow()

            # intersect transect lines with nlcd polygon to get length of transect lines overlapping natural land cover
            transect_natural = arcpy.analysis.PairwiseIntersect("{0};{1}".format(transect_intersect, natpoly),
                                                              os.path.join(scratch_gdb, "transect_natural"), "ALL",
                                                              output_type="LINE")

            transect_natural_lyr = arcpy.SelectLayerByLocation_management(transect_natural, "INTERSECT", input_line)

            # convert line transects fc to pandas dataframe to do calculations
            final_fields = ["Shape_Length"]
            transect_intersect_df = pd.DataFrame((row for row in arcpy.da.SearchCursor(transect_natural_lyr, final_fields)),
                                             columns=final_fields)

            # calculate buf1 by adding percent of natural land cover classes overlapping perimeter of area of interest
            if transect_intersect_df.empty:
                buf2 = 0
            else:
                buf2 = round((transect_intersect_df['Shape_Length'].sum()) / 8, 3)

            if buf2 >= 99:
                buf2_rank = "A"
            elif 75 <= buf2 < 99:
                buf2_rank = "B"
            elif 25 <= buf2 < 75:
                buf2_rank = "C"
            else:
                buf2_rank = "D"

            arcpy.AddMessage("- BUF2 value: " + str(buf2) + " BUF2 score: " + buf2_rank)
            arcpy.AddMessage("................................")

            # update values in the input dataset
            with arcpy.da.UpdateCursor(selected_input, eia_fields, where_clause="{0} = {1}".format(oid_fieldname, oid)) as cursor:
                for row in cursor:
                    row[0] = lan1
                    row[1] = lan1_rank
                    row[2] = lan2
                    row[3] = lan2_rank
                    row[4] = buf1
                    row[5] = buf1_rank
                    row[6] = buf2
                    row[7] = buf2_rank
                    cursor.updateRow(row)




