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

#tool_exec <- function(in_params)  #, out_params
#{

# check and load required libraries  
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
require(dplyr)
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
require(here)
if (!requireNamespace("arcgisbinding", quietly = TRUE)) install.packages("arcgisbinding")
require(arcgisbinding)
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
require(sf)
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
require(tidyr)
if (!requireNamespace("raster", quietly = TRUE)) install.packages("raster")
require(raster)


arc.check_product()

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

# open the NHA feature class and select and NHA
site <- arc.open(here(eia_gdb, "_TestSite"))
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

###########################################################
# LAN2 - 
###########################################################

ext <- extract(landcover_crop, site_buffer, method='simple')

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

# name the land uses
tabs$Var1<-factor(tabs$Var1, levels=c(0,1,9,13), labels=c("Water", "Green", "Shrubland", "Urban"))
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

lan2_combined_totscore <- (lan2_100m_totscore*0.6) + (lan2_500m_totscore*0.4)
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


###########################################################
# BUF1 - Perimeter with Natural Buffer
###########################################################



###########################################################
# BUF2 - 
###########################################################

#}