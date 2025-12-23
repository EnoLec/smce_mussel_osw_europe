#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%  SMCE - CHL  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Enora LECORDIER
# Date: 11/02/2025

# Notes: 
# Compute the fuzzy score for CHL for each year (2019-2023) following Filgueira et al (2011, DOI: 110.1016/j.seares.2011.04.006)
# See Excel file Models_for_thresholds for mathematical function (functional response)
# Also, this script averages the fuzzy score over the 5-year timeseries

# Initialisation ----
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/Users/my/path/")

# Compute fuzzy score ----
# List CHL files in folder
CHL_files <- list.files(pattern = "^CHL_20") 

# Split the files into 5 groups for the 5 years
file_groups <- split(CHL_files, cut(seq_along(CHL_files), 5, labels = FALSE))

chl_fuzzy <- c()  
for(i in 1:length(file_groups)){

  # Extract the year
  date <- substr(file_groups[[i]][1], 5, 8) 
  
  print(paste("Computing chl fuzzy method for:", date, paste0("(", i, "/", length(file_groups),")")))

  # Open CHL files per year
  CHL_daily <- lapply(file_groups[[i]], function(filename) { 
    raster::brick(filename)
  })
  
  CHL_daily <- stack(CHL_daily) # Stack daily data for one year a single object
  
  print("Raster stack generated. Now computing the fuzzy scores... This might take a while.")

  # Compute the function response on our dataset
  f <- CHL_daily / (CHL_daily + 1.06)  # Xk = 1.06 from  Filgueira et al (2011)

  # Average the result across the time serie
  chl_score <- mean(f, na.rm=T)
  
  chl_fuzzy[[i]] <- chl_score     # Fill the result list
  names(chl_fuzzy[[i]]) <- date   # Name the list object with the year 
}

chl_fuzzy <- stack(chl_fuzzy) 
writeRaster(chl_fuzzy,"CHL_2019-2023_fuzzy.tif", overwrite=T)


