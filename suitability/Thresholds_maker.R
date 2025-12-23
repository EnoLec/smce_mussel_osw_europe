# Initialisation ####
#library(dplyr)
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/Users/enora/OneDrive - hull.ac.uk/001_Enora_PhD/Data/Suitable_WF_GIS_study/MCE/CMEMS_data/")


# SMCE full Europe ----
sss_fuzzy <- raster("Thresholds/SSS/SSS_2019-2023_fuzzy.tif")
sst_fuzzy <- raster("Thresholds/SST/SST_2019-2023_fuzzy.tif")
chl_fuzzy <- raster("Thresholds/CHL/CHL_2019-2023_fuzzy.tif")
spm_fuzzy <- raster("Thresholds/SPM/SPM_2019-2023_fuzzy.tif")

res(sss_fuzzy)
res(sst_fuzzy)
res(chl_fuzzy)
res(spm_fuzzy)


sst_fuzzy <- resample(sst_fuzzy, sss_fuzzy)
chl_fuzzy <- resample(chl_fuzzy, sss_fuzzy)
spm_fuzzy <- resample(spm_fuzzy, sss_fuzzy)

sss <- mean(sss_fuzzy, na.rm=T)
sst <- mean(sst_fuzzy, na.rm=T)
chl <- mean(chl_fuzzy, na.rm=T)
spm <- mean(spm_fuzzy, na.rm=T)

SMCE_stack <- stack(sst, chl, sss, spm)

weights <- c(0.45,0.30,0.16,0.09)
SMCE_europe <- calc(SMCE_stack, fun = function(x) sum(x * weights))

writeRaster(SMCE_europe,"Thresholds/SMCE_full_europe_2019-2023.tif", overwrite=T)


# SMCE feasible Europe ----
feasibility <- raster("C:/Users/enora/OneDrive - hull.ac.uk/001_Enora_PhD/Data/Suitable_WF_GIS_study/MCE/feasibility/feasibility_mask.tif")
SMCE_europe <- raster("Thresholds/SMCE_full_europe_2019-2023.tif")
coast25NM <- read_sf("C:/Users/enora/OneDrive - hull.ac.uk/001_Enora_PhD/Data/Suitable_WF_GIS_study/MCE/GIS/Europe_coastline_shapefile/Europe_coastline_25NM_buffer.shp")

feasibility_res <- resample(feasibility, SMCE_europe, method="ngb")

SMCE_feasible <- SMCE_europe

SMCE_feasible[feasibility_res==0] <- NA

plot(SMCE_feasible)

writeRaster(SMCE_feasible,"Thresholds/SMCE_feasible_europe_2019-2023.tif", overwrite=T)


# SMCE future trend ----
## SMCE 2010-2020 ----
sst_present <- raster("Thresholds/SST/SST_2010-2020_fuzzy.tif")

sst_present <- resample(sst_present, sss_fuzzy)

SMCE_stack_present <- stack(sst_present, chl, sss, spm)
weights <- c(0.45,0.30,0.16,0.09)
SMCE_europe_present <- calc(SMCE_stack_present, fun = function(x) sum(x * weights))

## SMCE 2040-2050 ----
sst_future <- raster("Thresholds/SST/SST_2040-2050_fuzzy.tif")
sst_future <- resample(sst_future, sss_fuzzy)

SMCE_stack_future <- stack(sst_future, chl, sss, spm)
weights <- c(0.45,0.30,0.16,0.09)
SMCE_europe_future <- calc(SMCE_stack_future, fun = function(x) sum(x * weights))

## Compute trend ----
### Absolute future trend ----
SMCE_europe_trend <- SMCE_europe_future-SMCE_europe_present

# writeRaster(SMCE_europe_trend,"Thresholds/SMCE_full_europe_absolute_trend.tif", overwrite=T)

SMCE_feasible_trend <- SMCE_europe_trend

SMCE_feasible_trend[feasibility_res==0] <- NA

plot(SMCE_feasible_trend)
writeRaster(SMCE_feasible_trend,"Thresholds/SMCE_feasible_europe_absolute_trend.tif", overwrite=T)

### Relative future trend ----
SMCE_europe_trend <- 100*(SMCE_europe_future-SMCE_europe_present)/SMCE_europe_present

# writeRaster(SMCE_europe_trend,"Thresholds/SMCE_full_europe_absolute_trend.tif", overwrite=T)

SMCE_feasible_trend <- SMCE_europe_trend

SMCE_feasible_trend[feasibility_res==0] <- NA

plot(SMCE_feasible_trend)
writeRaster(SMCE_feasible_trend,"Thresholds/SMCE_feasible_europe_relative_trend.tif", overwrite=T)

# ### From palmer et el 2020 ####
# MCE_4.5 <- raster("Thresholds/SMCE_SSP4.5.tif")
# MCE_8.5 <- raster("Thresholds/SMCE_SSP8.5.tif")
# stability <- abs(MCE_8.5-MCE_4.5)/max(MCE_8.5,MCE_4.5)
# 
# writeRaster(stability,"Thresholds/stability_SSP.tif", overwrite=T)
# 
