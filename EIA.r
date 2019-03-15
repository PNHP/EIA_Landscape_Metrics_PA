#-------------------------------------------------------------------------------
# Name:        EIA.r
# Purpose:     Create an empty, new COA databases
# Author:      Christopher Tracey
# Created:     2019-03-14
# Updated:     
#
# Updates:
# * 2019-
# To Do List/Future ideas:
#
#-------------------------------------------------------------------------------

tool_exec <- function(in_params, out_params)  #
{

# check and load required libraries  
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
require(dplyr)
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
require(here)
#if (!requireNamespace("arcgisbinding", quietly = TRUE)) install.packages("arcgisbinding")
#require(arcgisbinding)
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
require(sf)
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
require(tidyr)
if (!requireNamespace("raster", quietly = TRUE)) install.packages("raster")
require(raster)


#arc.check_product()

# define functions
st_multibuffer <- function(x, from, to, by, nQuadSegs=30) {
  # get a multi buffered polygon. Require an sf object
  library(sf)
  library(dplyr)
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

## Network Paths and such
eia_gdb <- "EIA_Level1_Tool.gdb"
#planning_units <- in_params[[2]]

# open the NHA feature class and select and NHA
print("Prepping the data...")
sitename <- in_params[[1]]
site <- arc.open(sitename) # test site
#site <- in_params[[1]]
site <- arc.select(site)
siteID <- site$ID
site_sf <- arc.data2sf(site)


###########################################################
# Shared Data Prep
###########################################################

# multiple ring buffer
site_buffer <- st_multibuffer(site_sf, 100, 500, 400) # this produces two separate donuts. Se below for a combined one
site_buffer_all <- st_union(site_sf, site_buffer)
site_buffer_all <-st_union(site_buffer_all, by_feature = FALSE)

# extract landcover
landcover <- arc.raster(arc.open(here(eia_gdb, "nlcd2011")))
landcover <- as.raster(landcover)
landcover_crop <- mask(crop(landcover, as(site_buffer_all, 'Spatial')), as(site_buffer_all, 'Spatial'))

# get lookup table for natcover conversion
lu_nlcd2011_natcov <- arc.open(here(eia_gdb, "lu_NLCD2011_remapNatCov"))
lu_nlcd2011_natcov <- arc.select(lu_nlcd2011_natcov)
lu_nlcd2011_natcov$OBJECTID <- NULL
lu_nlcd2011_natcov_matrix <- as.matrix(lu_nlcd2011_natcov[c("Value","CoverTypeval")])

natcov <- reclassify(landcover_crop, lu_nlcd2011_natcov_matrix)

###########################################################
# LAN1 - Contigous Natural Buffer
###########################################################
arc.progress_label("Calculating LAN1...")
print("Calculating LAN1...")

natcov_poly <- rasterToPolygons(natcov, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
natcov_poly <- st_as_sf(natcov_poly)
natcov_poly_sub <- st_cast(natcov_poly, "POLYGON")

natcov_cont <- natcov_poly_sub[site_sf,]

area_total <- sum(as.numeric(st_area(natcov_poly)), na.rm=TRUE)
area_cont <- sum(as.numeric(st_area(natcov_cont)), na.rm=TRUE)

lan1_score <- round(area_cont/area_total, 3)
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

print(paste("LAN1 Score:", lan1_score,sep=" "))
print(paste("LAN1 Rating:", lan1_rating,sep=" "))


###########################################################
# LAN2 - land Use Index
###########################################################
# a lot of this is based on http://zevross.com/blog/2015/03/30/map-and-analyze-raster-data-in-r/
arc.progress_label("Calculating LAN2...")
print("Calculating LAN2...")

ext <- raster::extract(landcover_crop, site_buffer, method='simple')

# Function to tabulate land use by region and return a data.frame
tabFunc<-function(indx, extracted, region, regname) {
  dat<-as.data.frame(table(extracted[[indx]]))
  dat$name<-region[[regname]][[indx]]
  return(dat)
}

# run through each region and compute a table of the count
# of raster cells by land use. Produces a list (see below)
tabs<-lapply(seq(ext), tabFunc, ext, site_buffer, "distance")
# assemble into one data frame
tabs<-do.call("rbind",tabs )

#colnames(tabs)[colnames(tabs)=="old_name"] <- "cellcount"
tabs1 <- tabs %>% group_by(name) %>% mutate(per=round(Freq/sum(Freq), 3)) 

# get lookup table for landuse index conversion
lu_lan2coeff <- arc.open(here(eia_gdb, "lu_NLCD2011_LAN2coefficients"))
lu_lan2coeff <- arc.select(lu_lan2coeff)
lu_lan2coeff$OBJECTID <- NULL

tabs_final <- merge(tabs1, lu_lan2coeff[c("Value","coefficient")], by.x="Var1", by.y="Value")
tabs_final$score <- tabs_final$per * tabs_final$coefficient

tabs_final1 <- tabs_final[c("name","Var1", "score")]

# use the spread function from tidyr to make nicer
tabs_final2 <- tabs_final1 %>%
  group_by(Var1) %>% # group by region
  spread(key=name, value=score, fill=0) # make wide format

lan2_100m_totscore <- sum(tabs_final2$`100`)
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
lan2_500m_totscore <- sum(tabs_final2$`500`)
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

lan2_combined_totscore <- round((lan2_100m_totscore*0.6) + (lan2_500m_totscore*0.4),3)
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

print(paste("LAN2 Score:",lan2_combined_totscore, sep=" "))
print(paste("LAN2 Rating:",lan2_combined_rating, sep=" "))


###########################################################
# BUF1 - Perimeter with Natural Buffer
###########################################################
arc.progress_label("Calculating BUF1...")
print("Calculating BUF1...")

sf_ln <- st_cast(site_sf, 'LINESTRING')

ext1 <- raster::extract(natcov, sf_ln, method='simple')

a <- lapply(seq(ext1), tabFunc, ext1, sf_ln, "ID")
a2 <- do.call("rbind", a)

buf1_score <- a2$Freq[which(a2$Var1==1)] /sum(a2$Freq)
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

print(paste("Buf1 Score:",buf1_score,sep=" "))
print(paste("Buf1 Rating:",buf1_rating,sep=" "))


###########################################################
# BUF2 - 
###########################################################
arc.progress_label("Calculating BUF2...")
print("Calculating BUF2...")
site_100m <- site_buffer[which(site_buffer$distance==100),]

natcov_100m <- st_intersection(site_100m, natcov_cont)
natcov_100m <- st_union(site_sf, natcov_100m)
natcov_100m <- st_union(natcov_100m, by_feature = FALSE)

# convert to points
natcov_100m_pts <- st_cast(natcov_100m, to="MULTIPOINT")
natcov_100m_pts <- st_cast(natcov_100m_pts, to="POINT")  # may be good to thin these

natcov_dist <- st_distance(natcov_100m_pts, site_sf)
natcov_dist <- as.numeric(natcov_dist)
buf2_score <- round(mean(natcov_dist),3)

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
print(paste("Buf2 Score:",buf2_score,sep=" "))
print(paste("Buf2 Rating:",buf2_rating,sep=" "))

###########################################################
# Assemble results into a csv table 
###########################################################

EIA_results <- data.frame("Site ID"=siteID, "LAN1_Score"=lan1_score, "LAN1_Rating"=lan1_rating, "LAN2_Rating"=lan2_combined_totscore, "LAN2_Score"=lan2_combined_rating, "BUF1_Score"=buf1_score, "BUF1_Rating"=buf1_rating, "BUF2_Score"=buf2_score, "BUF2_Rating"=buf2_rating) 

write.csv(EIA_results, paste("site",siteID,"_EIAresults.csv",sep=""), row.names=FALSE)

}