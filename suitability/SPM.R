#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%  SMCE - SPM  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Enora LECORDIER
# Date: 11/02/2025

# Notes: 
# 1. Compute the average spm over the 2019-2023 timeseries (No need to run if not needed)
# 2. Compute the fuzzy score for spm for each year (2019-2023) following Kiorboe et al (1980, DOI: 10.1080/00785326.1980.10425516). 
#    See Excel file Models_for_thresholds for mathematical function
#    Also, this script averages the fuzzy score over the timeseries

# 0. Initialisation ####
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/my/path/")


# 1. Average spm 2019-2023 ----
## DO NOT RUN AGAIN ## Especially if not needed, as this can take a while. It is mainly there for mapping purpose if you need a average SPM raster to present your study site
# coastline <- read_sf("C:/Users/enora/OneDrive - hull.ac.uk/001_Enora_PhD/Data/Suitable_WF_GIS_study/MCE/GIS/Europe_coastline_shapefile/Europe_coastline_poly_4326.shp")
# SPM_files <- c("SPM_2019.nc","SPM_2020.nc","SPM_2021.nc","SPM_2022.nc","SPM_2023.nc")
# 
# # Read every .nc as brick file
# SPM_stack_list <- lapply(SPM_files, function(filename) {
#   raster::brick(filename, varname = "SPM")
# })
# 
# SPM_stack <- stack(SPM_stack_list) # Stack the bricks together to get only one file
# 
# temp <- mean(SPM_stack,  na.rm=TRUE)
# temp <- mask(temp, coastline, inverse=TRUE)
# 
# writeRaster(temp,"SPM_averaged_2019-2023.tif", overwrite=T)


# 2. Compute fuzzy score ----
SPM_files <- list.files(pattern = "^SPM_20") # List spm files in folder

spm_fuzzy <- c()                          # Generate result file
for(i in 1:length(SPM_files)){            # To loop through each year
  
  date <- substr(SPM_files[[i]][1], 5, 8) # Extract the year
  
  print(paste("Computing spm fuzzy method for:", date, paste0("(", i, "/", length(SPM_files),")")))
  
  SPM_monthly <- raster::brick(SPM_files[i])
  
  spm_score <-  -0.0002 * (SPM_monthly^2) + 0.0042 * SPM_monthly + 0.8273 # (see Excel file Models_for_thresholds) 
  
  spm_score[spm_score<0] <- 0 # If negative values (SPM out of range of the model leading to negative suitability) <- 0
  
  print("Fuzzy scores computed")
  
  spm_score <- mean(spm_score, na.rm=TRUE) # make yearly average of the fuzzy score
  
  spm_fuzzy[[i]] <- spm_score     # Fill the result list
  names(spm_fuzzy[[i]]) <- date   # Name the list object with the year 
  
}

spm_fuzzy <- stack(spm_fuzzy)

writeRaster(spm_fuzzy,"SPM_2019-2023_fuzzy.tif", overwrite=T)

