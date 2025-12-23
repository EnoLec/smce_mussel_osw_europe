#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%  SMCE - SSS  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Enora LECORDIER
# Date: 11/02/2025

# Notes: 
# 1. Compute the average sss over the 2019-2023 timeseries (No need to run again)
# 2. Compute the fuzzy score for sss for each year (2019-2023) following Maar et al (2015) and MMO (2019) data. 
#    See Excel file Models_for_thresholds for mathematical function
#    Also, this script averages the fuzzy score over the timeseries

# Maar et al 2015 - DOI: 10.1016/j.jmarsys.2015.02.003
# MMO 2019 - URL: https://assets.publishing.service.gov.uk/media/5dfb8f9840f0b6665e801834/MMO1184_AquaPotential_forPub_191210.pdf


# 0. Initialisation ----
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/Users/my/path/")


# 1. Average sss 2019-2023 ----
## DO NOT RUN AGAIN ##
# SSS_files <- list.files(pattern = "^SSS_20") # "^SSS_20" because we don't want the "SSS_monthly_201901-202312.nc" file generated in the folder with this script
# 
# SSS_001_031 <- SSS_files[1:4] # Read every .nc as brick file
# 
# SSS_files_001_031 <- lapply(SSS_001_031, function(filename) {
#   raster::brick(filename, varname = "so_cglo") # Open the data and convert kelvin in Celsius
# })
# 
# SSS_001_031 <- stack(SSS_files_001_031)
# 
# SSS_001_024 <- brick(SSS_files[5], varname = "so")
# 
# SSS_001_024 <- resample(SSS_001_024,SSS_001_031)
# 
# # Stack the bricks together to get only one file
# SSS_stack <- stack(SSS_001_031,SSS_001_024)
# temp <- mean(SSS_stack,  na.rm=TRUE)
# 
# writeRaster(temp,"SSS_averaged_2019-2023.tif", overwrite=T)

# 2. Compute fuzzy score ----
SSS_files <- list.files(pattern = "^SSS_20") # List sss files in folder

resampler <- raster("SSS_2019.nc") # Load only the first layer of SSS_2019 to resample 2023 on a 0.25 resolution (currently 0.08)

sss_fuzzy <- c()
for(i in 1:length(SSS_files)){
  
  date <- substr(SSS_files[[i]][1], 5, 8) # Extract the year
  
  print(paste("Computing sss fuzzy method for:", date, paste0("(", i, "/", length(SSS_files),")")))
  
  # Open SSS files per year
  SSS_monthly <- raster::brick(SSS_files[i])
  
  sss_score <-  -0.0000085283 * SSS_monthly^4 + 0.0006914926 * SSS_monthly^3 - 0.0208964469 * SSS_monthly^2 + 0.2971865272 * SSS_monthly - 0.8758292861 # DW/WW Ratio (see Excel file Models_for_thresholds) 
  sss_score[sss_score<0] <- 0 # If negative values (SSS out of range of the model leading to negative suitability) <- 0
  sss_score[sss_score>1] <- 1 # mathematical function maximum is 1.00117. Flatten to 1 for practicality
  
  # Resample if needed (2023 only)
  if(ext(sss_score) != ext(resampler)){
    sss_score <- resample(sss_score,resampler)
  }
  
  print("Fuzzy scores computed")
  
  sss_score <- mean(sss_score, na.rm=TRUE) # make yearly average of the fuzzy score
  
  sss_fuzzy[[i]] <- sss_score     # Fill the result list
  names(sss_fuzzy[[i]]) <- date   # Name the list object with the year 
  
}

sss_fuzzy <- stack(sss_fuzzy)     
writeRaster(sss_fuzzy,"Thresholds/SSS/SSS_2019-2023_fuzzy.tif", overwrite=T)



