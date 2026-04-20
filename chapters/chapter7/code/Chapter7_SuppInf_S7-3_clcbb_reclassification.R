# -----------------------------------------
# From the work: Illanas, et al (in prep). 
# Large-scale modeling of wild ungulate relative abundance 
# from a national camera trap network
# -----------------------------------------
# script: Chapter7_SuppInf_S7-3_clcbb_reclassification.R
# This script was used for extract corine land cover plus classes
# author: Sonia Illanas
# v1.0.0
# date last modification: 2026/03/17
# -----------------------------------------

# clc Land cover
library(terra)       # v1.8-50
library(tidyverse)   # v2.0.0
library(sf)          # v1.0-21
rm(list=ls())

# load the raster file downloaded from CLCplus
clc <- rast("CLCplus_RASTER_2021_010m_03035.tif")

needle<-clc==2
plot(needle)
terra::writeRaster(needle, "clc_NeedleLeavedTrees.tif")

bl_decid<-clc==3
terra::writeRaster(bl_decid, "clc_broadLeavedDTrees.tif")

bl_evergreen<-clc==4
terra::writeRaster(bl_evergreen, "clc_broadLeavedEGTrees.tif")

shrubs<-clc==5
terra::writeRaster(shrubs, "clc_shrubs.tif")

forest<-needle + bl_decid + bl_evergreen
forest
# plot(forest)
terra::writeRaster(forest, "clc_Needle_broadLeavedDecAndEGTrees.tif")
