#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%  SMCE - SPM  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Enora LECORDIER
# Date: 11/02/2025

# Notes: 
# Compute the fuzzy score for spm for each year (2019-2023) following Kiorboe et al (1980, DOI: 10.1080/00785326.1980.10425516). 
# See Excel file Models_for_thresholds for mathematical function
# Also, this script averages the fuzzy score over the timeseries

# Initialisation ####
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/my/path/")

# Compute fuzzy score ----
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


