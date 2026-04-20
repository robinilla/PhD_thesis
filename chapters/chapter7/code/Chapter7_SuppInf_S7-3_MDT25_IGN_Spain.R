# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: Chapter7_SuppInf_S7-3_MDT25_IGN_Spain.R
# This script was used for joining MDT25 files in Spain & 
# generating elevation and slope covariate rasters
# 
# tiles from a same UTM zone has been previously joined with QGIS
# and saved with the same UTM zone.
# Note that there is information from two flights. For some parts of 
# Spain the second flight had empty information so the information 
# from the first flight was needed from these places 
# 
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------
library(terra)       # v1.8-50
library(tidyverse)   # v2.0.0

rm(list=ls())
# load the rasters from each UTM zone
r1 <- rast("~/MDT25_IGN_ES_1_2_Cob/MDT25_Espana_Final_H29.tif")
r2 <- rast("~/MDT25_IGN_ES_1_2_Cob/MDT25_Espana_Final_H30.tif")
r3 <- rast("~/MDT25_IGN_ES_1_2_Cob/MDT25_Espana_Final_H31.tif")
r1_ <- rast("~/MDT25_1Cob_Espana_Final_H29.tif")
r2_ <- rast("~/MDT25_1Cob_Espana_Final_H30.tif")

# defined their UTM zone. They are all already "EPSG:25830"
crs(r1) <- "EPSG:25830" # Zone 29
crs(r1_)<- "EPSG:25830" # Zone 29
crs(r2) <- "EPSG:25830" # Zone 30
crs(r2_)<- "EPSG:25830" # Zone 30
crs(r3) <- "EPSG:25830" # Zone 31

# join the raster layers in the same object
r_col <- sprc(list(r1, r1_, r2, r2_, r3))
r_final_mosaic <- mosaic(r_col, fun = "mean"); # rm(r1_proj, r1.1_proj, r2, r2_, r3_proj)

# save the MDT raster 
# writeRaster(r_final_mosaic, 
#             "MDT25_allMerged_Espana_Peninsular_H30.tif", 
#             overwrite = TRUE,
#             wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")))


# calculate slope 
slope <- terrain(r_final, v = "slope", unit = "radians")
slope_pct <- tan(slope) * 100 # transform to percentage degrees

# save slope raster
# writeRaster(slope_pct, 
#             "MDT25_ES_PerSlope.tif", 
#             overwrite = TRUE,
#             wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")))
