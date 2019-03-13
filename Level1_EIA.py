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

################################################################################
## LAN1  Contiguous Natural Buffer #############################################

# region group
arcpy.AddMessage("Create Region Groups.")
outRgnGrp = RegionGroup(NatCov, "FOUR")
# extract to points and extract from region group
tmp_points = "outSitePoints"
tmp_points1 = "outSitePoints_extract"
arcpy.FeatureToPoint_management(site, tmp_points, "INSIDE")
ExtractValuesToPoints(tmp_points, outRgnGrp, tmp_points1, "NONE", "ALL")
arcpy.AddMessage("Trying to get regiongroup id")

with arcpy.da.SearchCursor(tmp_points1, ['RASTERVALU']) as cursor:
    for row in cursor:
        extractvalue = str(row[0])

inSQLClause = "Value ={0}".format(extractvalue)
attRgnGrp = ExtractByAttributes(outRgnGrp, inSQLClause)
attRgnGrp.save("tmp_ContNatCover")

# calculate the score
attRgnGrp = arcpy.BuildRasterAttributeTable_management(attRgnGrp, "Overwrite")
with arcpy.da.SearchCursor(attRgnGrp, ['Count']) as cursor:
    for row in cursor:
        NatCovCellCount = row[0]

outExtractByMask1 = arcpy.BuildRasterAttributeTable_management(outExtractByMask, "Overwrite")
AllCellCount = 0
with arcpy.da.SearchCursor(outExtractByMask1, ['Count']) as cursor:
    for row in cursor:
        AllCellCount += row[0]

NatCovProp = NatCovCellCount / AllCellCount
NatCovProp = round(NatCovProp, 3)

# calculate the rating
if NatCovProp >= 0.9:
    NatCovScore = "A"
elif NatCovProp >= 0.6 and NatCovProp < 0.9:
    NatCovScore = "B"
elif NatCovProp >= 0.2 and NatCovProp < 0.6:
    NatCovScore = "C"
else:
    NatCovScore = "D"

# report it out
arcpy.AddMessage("----------------------------------------")
arcpy.AddMessage("LAN1 - Contiguous Natural Buffer -------")
arcpy.AddMessage("Score: {0}".format(str(NatCovProp)))
arcpy.AddMessage("Rating: {0}".format(str(NatCovScore)))
arcpy.AddMessage("----------------------------------------")


################################################################################
## LAN2 - Land Use Index #######################################################




# report it out
arcpy.AddMessage("----------------------------------------")
arcpy.AddMessage("LAN2 - Land Use Index ------------------")
#arcpy.AddMessage("Score: {0}".format(str(NatCovProp)))
#arcpy.AddMessage("Rating: {0}".format(str(NatCovScore)))
arcpy.AddMessage("----------------------------------------")




################################################################################
## BUF1 - Perimeter with Natural Buffer ########################################




# report it out
arcpy.AddMessage("----------------------------------------")
arcpy.AddMessage("BUF1 -  Perimeter with Natural Buffer---")
#arcpy.AddMessage("Score: {0}".format(str(NatCovProp)))
#arcpy.AddMessage("Rating: {0}".format(str(NatCovScore)))
arcpy.AddMessage("----------------------------------------")







################################################################################
## BUF2 - Width of Natural Buffer ##############################################




# report it out
arcpy.AddMessage("----------------------------------------")
arcpy.AddMessage("BUF2 - Width of Natural Buffer ---------")
#arcpy.AddMessage("Score: {0}".format(str(NatCovProp)))
#arcpy.AddMessage("Rating: {0}".format(str(NatCovScore)))
arcpy.AddMessage("----------------------------------------")

