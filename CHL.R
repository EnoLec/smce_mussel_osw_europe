#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%  SMCE - CHL  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Enora LECORDIER
# Date: 11/02/2025

# Notes: 
# 1. Compute the average chl over the 2019-2023 timeseries (No need to run again)
# 2. Compute the fuzzy score for chl for each year (2019-2023) following Filgueira et al (2011, Equation 2b)
#    See Excel file Models_for_thresholds for mathematical function (functional response)
#    Also, this script averages the fuzzy score over the timeseries

# Filgueira et al 2011 - DOI: 110.1016/j.seares.2011.04.006


# 0. Initialisation ----
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/Users/enora/OneDrive - hull.ac.uk/001_Enora_PhD/Data/Suitable_WF_GIS_study/MCE/CMEMS_data/")


# 1. Average chl 2019-2023----
# CHL_files <- list.files(pattern = "^CHL_20") 
# 
# CHL_files <- lapply(CHL_files, function(filename) {
#   raster::brick(filename, varname = "CHL") # Open the data and convert kelvin in Celsius
# })
# 
# CHL_stack <- stack(CHL_files)
# 
# temp <- mean(CHL_stack,  na.rm=TRUE)
# temp <- mask(temp, coastline, inverse=TRUE)
# 
# writeRaster(temp,"CHL_averaged_2019-2023.tif", overwrite=T)


# 2. Compute fuzzy score ----
CHL_files <- list.files(pattern = "^CHL_20") # List chl files in folder

# Split the files into 5 groups for the 5 years
file_groups <- split(CHL_files, cut(seq_along(CHL_files), 5, labels = FALSE))

chl_fuzzy <- c()  
for(i in 1:length(file_groups)){
  
  date <- substr(file_groups[[i]][1], 5, 8) # Extract the year
  
  print(paste("Computing chl fuzzy method for:", date, paste0("(", i, "/", length(file_groups),")")))
  
  CHL_daily <- lapply(file_groups[[i]], function(filename) { # Open CHL files per year
    raster::brick(filename)
  })
  
  CHL_daily <- stack(CHL_daily) # One stack for every day of the year
  
  print("Raster stack generated. Now computing the fuzzy scores... This might take a while.")
  
  f <- CHL_daily / (CHL_daily + 1.06)  # Xk = 1.06 from  Filgueira et al (2011, Table 2)
  
  chl_score <- mean(f, na.rm=T)
  
  chl_fuzzy[[i]] <- chl_score     # Fill the result list
  names(chl_fuzzy[[i]]) <- date   # Name the list object with the year 
}


chl_fuzzy <- stack(chl_fuzzy) 
writeRaster(chl_fuzzy,"Thresholds/CHL/CHL_2019-2023_fuzzy.tif", overwrite=T)

