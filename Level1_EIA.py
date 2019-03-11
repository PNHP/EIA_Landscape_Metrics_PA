#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      ctracey
#
# Created:     2018/03/06
# Copyright:   (c) ctracey 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#import packages
import os, arcpy, datetime
from arcpy import env
from arcpy.sa import *

# set directories and env variables
arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = False

arcpy.env.workspace = "EIA_Level1_Tool.gdb"

# set input paramenters for Tool
site = arcpy.GetParameterAsText(0)
landcover = arcpy.GetParameterAsText(1)
#outGDB = arcpy.GetParameterAsText(1)

site_buffer = "out_MRbuff"

# create a multiple ring buffer on the site
arcpy.AddMessage("Creating the multiple ring buffer.")
arcpy.MultipleRingBuffer_analysis(site, site_buffer, [100,500], "meters", "", "ALL", "FULL")
#THIS ISN"T WORKING # site_buffer.save("SiteBuffer")

# Set Mask environment
#arcpy.env.mask = site_buffer

# Execute ExtractByMask
arcpy.AddMessage("Extracting the landcover data.")
outExtractByMask = ExtractByMask(landcover, site_buffer)

# load the various lookup tables
lu_reclass = "lu_NLCD2011_remapNatCov"
# Execute Reclassify
arcpy.AddMessage("Reclassify the landcover to natural cover.")
NatCov = ReclassByTable(outExtractByMask, lu_reclass, "Value", "Value", "CoverTypeval", "NODATA")
NatCov.save("Site_NatCover")

# region group
arcpy.AddMessage("Create Region Groups.")
outRgnGrp = "tmp_RegionGroups"
outRgnGrp = RegionGroup(NatCov, "FOUR")
# extract to points and extract from region group
tmp_points = "outSitePoints"
tmp_points1 = "outSitePoints_extract"
arcpy.FeatureToPoint_management(site, tmp_points, "INSIDE")
ExtractValuesToPoints(tmp_points, outRgnGrp, tmp_points1, "NONE", "ALL")
arcpy.AddMessage("Trying to get regiongroup id")

cursor = arcpy.da.SearchCursor(tmp_points1, ['RASTERVALU'])
for row in cursor:
    print(row[0])
    inSQLClause = '"' + "Value" + '"  = "' + str(row[0]) + '"'
    attRgnGrp = ExtractByAttributes(inRaster, inSQLClause)
    attExtract.save("tmp_ContNatCover")

#with arcpy.da.SearchCursor(tmp_points1, ['RASTERVALU']) as cursor:
#    for row in cursor:
#        x = {0}
#        #print('{0}'.format(row[0]))

#attRgnGrp = ExtractByAttributes(outRgnGrp, '"' + "Value" + '"  = "' + str(x) + '"' )





