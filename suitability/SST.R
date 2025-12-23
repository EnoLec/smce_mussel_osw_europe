#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%  SMCE - SST  %%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Enora LECORDIER
# Date: 11/02/2025

# Notes: 
# 1. Compute the fuzzy score for sst for each year (2019-2023) following Lauzon-Guay et al (2006, DOI: 10.3354/meps323171)
#    and its average over the timeseries
# 2. Compute the fuzzy score same way than 1. for Bio-Oracle data on 2010-2020 and 2040-2050 timeseries


# 0. Initialisation ----
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/my/path/")

# Variables from Lauzon-Guay et al (2006)
T_max <- 30     # Upper temperature at which growth stops
T_opt <- 15.8   # Optimal temperature for maximum growth
a <- 1.14       # Scale for growth
b <- -0.547     # Allometric exponent
c <- 0.393      # Constant related to temperature dependency

XTmax_exclusion = function(x, y) {    # Takes the same pixels of fT and XTmax
  x[y[]==1] <- 0                      # If a pixel in XTmax == 1, then the same pixel in fT <- 0
  return(x)
}

# 1. Compute fuzzy score ----
SST_files <- list.files(pattern = "^SST_20") # List sst files in folder

file_groups <- split(SST_files, cut(seq_along(SST_files), 5, labels = FALSE)) # Split the files into 5 groups for the 5 years

sst_fuzzy <- c()                            # Generate result file
for(i in 1:length(file_groups)){            # To loop through each year
  
  date <- substr(file_groups[[i]][1], 5, 8) # Extract the year
  
  print(paste("Computing sst fuzzy method for:", date, paste0("(", i, "/", length(file_groups),")")))
  
  SST_daily <- lapply(file_groups[[i]], function(filename) { # Open SST files per year
    raster::brick(filename, varname = "analysed_sst")-272.15 # Open the data and convert kelvin in Celsius
  })
  
  SST_daily <- stack(SST_daily) # One stack for every day of the year
  
  print("Raster stack generated. Now computing the fuzzy scores... This might take a while.")
  
  # Compute the upper section of Equation 2 in Lauzon-Guay (2006)
  A <- (T_max-SST_daily)/(T_max-T_opt)
  B <- c*(T_max-T_opt)
  C <- exp(c*(SST_daily-T_opt))
  
  fT <- (A^B)*C                # Compute fT whatever T_max
  XTmax <- SST_daily >= T_max  # Spot when T>=T_max
  XTmax[XTmax != 0]  <- 1      # If XTmax != 0 -> the threshold is crossed
  
  # When XTmax == 1, that means that the threshold of T_max was crossed for that day -> the fT value is probably already very low, 
  # but need to be absolute 0 to respect Lauzon-Guay's second part of the model.
  
  sst_score <- c()                   # Generate an empty score object
  for (j in 1:nlayers(SST_daily)){   # Loop through each day
    sst_score[[j]] <- overlay(fT[[j]], XTmax[[j]], fun = XTmax_exclusion)
  }
  
  print("Fuzzy scores computed")
  
  sst_score <- stack(sst_score)            # stack the days together for the entire year
  sst_score <- mean(sst_score, na.rm=TRUE) # make yearly average of the fuzzy score
  
  sst_fuzzy[[i]] <- sst_score     # Fill the result list
  names(sst_fuzzy[[i]]) <- date   # Name the list object with the year 

}

sst_fuzzy <- stack(sst_fuzzy)     # stack sst_fuzzy
writeRaster(sst_fuzzy,"SST_2019-2023_fuzzy.tif", overwrite=T)


# 2. Bio-oracle timeseries ----
## 2010-2020 ----
SST <- raster("Bio-oracle_2010-2020_meanSST_clip.tif") 
e <- as(extent(-12, 34, 33, 73), 'SpatialPolygons')
SST <-crop(SST,e)

# Fit the model
A <- (T_max-SST)/(T_max-T_opt)
B <- c*(T_max-T_opt)
C <- exp(c*(SST-T_opt))
fT <- (A^B)*C

XTmax <- SST >= T_max
XTmax[XTmax != 0]  <- 1 

sst_score <- c()
sst_score <- overlay(fT, XTmax, fun=XTmax_exclusion)

writeRaster(sst_score,"SST_2010-2020_fuzzy.tif", overwrite=T)


## 2040-2050 ----
SST <- raster("Bio-oracle_SSP85_2040-2050_meanSST.nc") 
e <- as(extent(-12, 34, 33, 73), 'SpatialPolygons')
SST <-crop(SST,e)

# Fit the model
A <- (T_max-SST)/(T_max-T_opt)
B <- c*(T_max-T_opt)
C <- exp(c*(SST-T_opt))
fT <- (A^B)*C

XTmax <- SST >= T_max
XTmax[XTmax != 0]  <- 1 

sst_score <- c()
sst_score <- overlay(fT, XTmax, fun=XTmax_exclusion)

writeRaster(sst_score,"SST/SST_2040-2050_fuzzy.tif", overwrite=T)
