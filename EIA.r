#-------------------------------------------------------------------------------
# Name:        EIA.r
# Purpose:     This tool runs in ArcGIS Pro and creates a table of Ecological
#              Integrity Assessment (EIA) Scores and Ratings for wetland 
#              plant communities. The tool exports this data into a csv file for
#              import into EcoObs.
# Author:      Christopher Tracey
# Created:     2019-03-14
# Updated:     2019-03-16
#
# Updates:
# * 2019-03-16 comments in script, and misc cleanup
#
# To Do List/Future ideas:
# * vectorize to deal with more than one site at a time
# * ability to select the id field, but export it at whatever EcoObs uses
# * return buffer to ArcGIS???
#-------------------------------------------------------------------------------

tool_exec <- function(in_params, out_params)
{

# check and load required libraries  
 if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
require(dplyr)
if (!requireNamespace("arcgisbinding", quietly = TRUE)) install.packages("arcgisbinding")
  require(arcgisbinding)
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
  require(sf)
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  require(tidyr)
if (!requireNamespace("raster", quietly = TRUE)) install.packages("raster")
  require(raster)

#arc.check_product()
    
###########################################################
# define functions
###########################################################
# Function to create a multiple ring buffer
st_multibuffer <- function(x, from, to, by, nQuadSegs=30) { # get a multi buffered polygon. Require an sf object
  seq.buffer <-seq(from,to, by) # allocate a  vector of length equal to number of buffer required
  if(inherits(x,"sf")) {
    # create a list of buffers 
    k <- vector('list', length(seq.buffer))
    for (i in 1:length(seq.buffer)) {
      k[[i]] <- st_buffer(x =x, dist = seq.buffer[i], nQuadSegs = nQuadSegs)
      k[[i]]$idx <- i
      k[[i]]$distance <- seq.buffer[i]
      st_agr(k[[i]]) <- "constant"
    }
    # clip from the bigger to the lower 
    l <- vector('list', length(seq.buffer))
    for (i in length(seq.buffer):1) {
      if(i==1){
        l[[i]] <- st_difference(k[[i]], st_geometry(x))
      } 
      else {
        l[[i]] <- st_difference(k[[i]], st_geometry(k[[i - 1]]))
      }
    }
    # Join all multipolygon
    if (length(seq.buffer) == 2) {
      temp <- rbind(l[[1]], l[[2]])
    } else {
      temp <- rbind(l[[1]], l[[2]])
      for (m in 3:length(seq.buffer)) {
        temp <- rbind(temp, l[[m]])
      }
    }
    return(temp)
  } else {
    stop("x is not a sf object")
  }
}

# Function to tabulate land use by region and return a data.frame
tabFunc <- function(indx, extracted, region, regname) {
  dat <- as.data.frame(table(extracted[[indx]]))
  dat$name <- region[[regname]][[indx]]
  return(dat)
}

###########################################################
#
###########################################################
arc.progress_pos(0)
arc.progress_label("Getting data and variables...")

setwd("E:/Misc/EIA_Level1_Tool")
# Network Paths and such
eia_gdb <- "EIA_Level1_Tool.gdb"

# open the NHA feature class and select and NHA
print("Prepping the data...")
sitename <- in_params[[1]] # sitename <- here(eia_gdb,"_TestSite")
site <- arc.open(sitename) 
site <- arc.select(site)
siteID <- site$ID # use in_params[[2]] to select site id field
site_sf <- arc.data2sf(site) # convert to a simple features object

###########################################################
# Shared Data Prep - this preps the site buffers and clips out the landcover to the buffer
###########################################################

# multiple ring buffer
site_buffer <- st_multibuffer(site_sf, 100, 500, 400) # this produces two separate donuts with a hole in the center
site_buffer_all <- st_union(site_sf, site_buffer) # this fills in the hole in the center
site_buffer_all <-st_union(site_buffer_all, by_feature = FALSE) # dissolves the boundaries
fgdb_path <- file.path(eia_gdb)
arc.write(file.path(fgdb_path,paste("Site",siteID,"Buffer",sep="_")), data=site_buffer, overwrite=TRUE)


# extract landcover
landcover <- arc.open(paste0(eia_gdb, "/nlcd2011"))
landcover1 <- arc.raster(landcover)
landcover1 <- as.raster(landcover1)
landcover_mask <- crop(landcover1, as(site_buffer_all, 'Spatial'))
landcover_crop <- mask(landcover_mask, as(site_buffer_all, 'Spatial'))
print("works to here!")
# convert landcover to natural cover based on a lookup table
lu_nlcd2011_natcov <- arc.open(paste0(eia_gdb, "/lu_NLCD2011_remapNatCov"))
lu_nlcd2011_natcov <- arc.select(lu_nlcd2011_natcov)
lu_nlcd2011_natcov_matrix <- as.matrix(lu_nlcd2011_natcov[c("Value","CoverTypeval")])
natcov <- reclassify(landcover_crop, lu_nlcd2011_natcov_matrix) # reclassify to natural or not
rm(lu_nlcd2011_natcov, lu_nlcd2011_natcov_matrix)

###########################################################
# LAN1 - Contigous Natural Buffer
###########################################################
arc.progress_pos(20)
arc.progress_label("Calculating LAN1...")
print("Calculating LAN1...")

# convert the natural cover layer to a polygon
natcov_poly <- rasterToPolygons(natcov, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
natcov_poly <- st_as_sf(natcov_poly)
natcov_poly_sub <- st_cast(natcov_poly, "POLYGON")

# create a layer of continous natural cover by finding which features intersect with original site
natcov_cont <- natcov_poly_sub[site_sf,]

# caculate the area of each and then create the LAN1 score
area_total <- sum(as.numeric(st_area(natcov_poly)), na.rm=TRUE)
area_cont <- sum(as.numeric(st_area(natcov_cont)), na.rm=TRUE)
lan1_score <- round(area_cont/area_total, 3)
# calculate the rating
lan1_rating <- NA
if(lan1_score >= 0.9){
  lan1_rating <- "A"
} else if (lan1_score >= 0.6 && lan1_score < 0.9){
  lan1_rating <- "B"
} else if(lan1_score >= 0.2 && lan1_score < 0.6){
  lan1_rating <- "C"
} else {
  lan1_rating <- "D"
}
# return the results to the screen
print(paste("LAN1 Score:", lan1_score,sep=" "))
print(paste("LAN1 Rating:", lan1_rating,sep=" "))

###########################################################
# LAN2 - land Use Index
###########################################################
# a lot of this is based on http://zevross.com/blog/2015/03/30/map-and-analyze-raster-data-in-r/
arc.progress_pos(40)
arc.progress_label("Calculating LAN2...")
print("Calculating LAN2...")
# extract the landcover raster
ext <- raster::extract(landcover_crop, site_buffer, method='simple') # call the package first since dplyr is messing it up

# run through each region and compute a table of the count of raster cells by land use. Produces a list (see below)
tabs <- lapply(seq(ext), tabFunc, ext, site_buffer, "distance")
# assemble into one data frame
tabs <- do.call("rbind",tabs )

#colnames(tabs)[colnames(tabs)=="old_name"] <- "cellcount"
tabs1 <- tabs %>% group_by(name) %>% mutate(per = round(Freq/sum(Freq), 3)) 

# get lookup table for landuse index conversion
lu_lan2coeff <- arc.open(paste0(eia_gdb, "/lu_NLCD2011_LAN2coefficients"))
lu_lan2coeff <- arc.select(lu_lan2coeff)
lu_lan2coeff$OBJECTID <- NULL

tabs_final <- merge(tabs1, lu_lan2coeff[c("Value","coefficient")], by.x="Var1", by.y="Value")
tabs_final$score <- tabs_final$per * tabs_final$coefficient
tabs_final <- tabs_final[c("name","Var1", "score")] # subset to needed columns
# use the spread function from tidyr to make nicer
tabs_final <- tabs_final %>%
  group_by(Var1) %>% # group by region
  spread(key=name, value=score, fill=0) # make wide format

# calculate the score for the inner ring
lan2_100m_totscore <- sum(tabs_final$`100`)
# calculate the rating
lan2_100m_rating <- NA
if(lan2_100m_totscore >= 9.5){
  lan2_100m_rating <- "A"
} else if (lan2_100m_totscore >= 8 && lan2_100m_totscore < 9.5){
  lan2_100m_rating <- "B"
} else if(lan2_100m_totscore >= 4&& lan2_100m_totscore < 8){
  lan2_100m_rating <- "C"
} else {
  lan2_100m_rating <- "D"
}
lan2_500m_totscore <- sum(tabs_final$`500`)
# calculate the rating for the outer ring
lan2_500m_rating <- NA
if(lan2_500m_totscore >= 9.5){
  lan2_500m_rating <- "A"
} else if (lan2_500m_totscore >= 8 && lan2_500m_totscore < 9.5){
  lan2_500m_rating <- "B"
} else if(lan2_500m_totscore >= 4&& lan2_500m_totscore < 8){
  lan2_500m_rating <- "C"
} else {
  lan2_500m_rating <- "D"
}

# calculate the combined score for the inner and outer rings
lan2_combined_totscore <- round((lan2_100m_totscore*0.6) + (lan2_500m_totscore*0.4),3)
# calculate the rating
lan2_combined_rating <- NA
if(lan2_combined_totscore >= 9.5){
  lan2_combined_rating <- "A"
} else if (lan2_combined_totscore >= 8 && lan2_combined_totscore < 9.5){
  lan2_combined_rating <- "B"
} else if(lan2_combined_totscore >= 4&& lan2_combined_totscore < 8){
  lan2_combined_rating <- "C"
} else {
  lan2_combined_rating <- "D"
}
# return the results to the screen
print(paste("LAN2 Score:",lan2_combined_totscore, sep=" "))
print(paste("LAN2 Rating:",lan2_combined_rating, sep=" "))

###########################################################
# BUF1 - Perimeter with Natural Buffer
###########################################################
arc.progress_pos(60)
arc.progress_label("Calculating BUF1...")
print("Calculating BUF1...")
# conver to a line feature from a polygon
sf_ln <- st_cast(site_sf, 'LINESTRING')
# extract the raster cells that intersect the line
ext1 <- raster::extract(natcov, sf_ln, method='simple')

# i don't know that the next two lines do???????
a <- lapply(seq(ext1), tabFunc, ext1, sf_ln, "ID")
a2 <- do.call("rbind", a)

# calculate the rating
buf1_score <- a2$Freq[which(a2$Var1==1)] /sum(a2$Freq)
# return the results to the screen
buf1_rating <- NA
if(buf1_score >= 0.9){
  buf1_rating <- "A"
} else if (buf1_score >= 0.75 && buf1_score < 0.9){
  buf1_rating <- "B"
} else if(buf1_score >= 0.25 && buf1_score < 0.75){
  buf1_rating <- "C"
} else {
  buf1_rating <- "D"
}
# return the results to the screen
print(paste("BUF1 Score:",buf1_score,sep=" "))
print(paste("BUF1 Rating:",buf1_rating,sep=" "))

###########################################################
# BUF2 - 
###########################################################
arc.progress_pos(80)
arc.progress_label("Calculating BUF2...")
print("Calculating BUF2...")
# subset the 100m buffer from the multiple ring
site_100m <- site_buffer[which(site_buffer$distance==100),]
# clip the contigious natural cover buffer by the 100m buffer fill in the donut hole
natcov_100m <- st_intersection(site_100m, natcov_cont)
natcov_100m <- st_union(site_sf, natcov_100m)
natcov_100m <- st_union(natcov_100m, by_feature = FALSE)
arc.write(file.path(fgdb_path,paste("Site",siteID,"NatCover",sep="_")), data=natcov_100m, overwrite=TRUE)
# convert to polygon to points to use in the distance matrix
natcov_100m_pts <- st_cast(natcov_100m, to="MULTIPOINT")
natcov_100m_pts <- st_cast(natcov_100m_pts, to="POINT")  # may be good to thin these
rm(natcov_100m)
# calculate of vector of distances
natcov_dist <- st_distance(natcov_100m_pts, site_sf)
natcov_dist <- as.numeric(natcov_dist)
rm(natcov_100m_pts)
# calculate the BUF2 score
buf2_score <- round(mean(natcov_dist),3)
# calculate the rating
buf2_rating <- NA
if(buf2_score >= 100){
  buf2_rating <- "A"
} else if (buf2_score >= 75 && buf2_score < 100){
  buf2_rating <- "B"
} else if(buf2_score >= 25 && buf2_score < 75){
  buf2_rating <- "C"
} else {
  buf2_rating <- "D"
}
# return the results to the screen
print(paste("BUF2 Score:",buf2_score,sep=" "))
print(paste("BUF2 Rating:",buf2_rating,sep=" "))

###########################################################
# Assemble results into a csv table 
###########################################################
arc.progress_pos(90)
arc.progress_label("Assembing data into a table...")

# create a dataframe of the above results
EIA_results <- data.frame("Site ID"=siteID, "LAN1_Score"=lan1_score, "LAN1_Rating"=lan1_rating, "LAN2_Rating"=lan2_combined_totscore, "LAN2_Score"=lan2_combined_rating, "BUF1_Score"=buf1_score, "BUF1_Rating"=buf1_rating, "BUF2_Score"=buf2_score, "BUF2_Rating"=buf2_rating) 

# write the csv to the working directory
write.csv(EIA_results, paste("site",siteID,"_EIAresults.csv",sep=""), row.names=FALSE)

}