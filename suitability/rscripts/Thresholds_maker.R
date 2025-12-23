#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%% FEASIBILITY ANALYSIS %%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Here we combine CHL, SST, SSS and SPM suitability scores together
# First, we are going to generate a suitability index of M. edulis growth for the whole of Europe (without feasibility mask applied), to allow a complete study of the region.
# Then, we will aplpy the Feasiblity mask to identify the best areas.
# Finally, we'll compute the suitability index (SI) trend by 2050.

# Initialisation ####
library(raster)
library(sp)
library(sf)
library(ncdf4)
library(terra)

setwd("C:/Users/my/path/")

# SMCE full Europe ----
# Load the suitability score rasters for each factor
sss_fuzzy <- raster("SSS_2019-2023_fuzzy.tif")
sst_fuzzy <- raster("SST_2019-2023_fuzzy.tif")
chl_fuzzy <- raster("CHL_2019-2023_fuzzy.tif")
spm_fuzzy <- raster("SPM_2019-2023_fuzzy.tif")

# Check their resolution
res(sss_fuzzy)
res(sst_fuzzy)
res(chl_fuzzy)
res(spm_fuzzy)

# Resample every layers on on the layer of your choice, here SSS 
sst_fuzzy <- resample(sst_fuzzy, sss_fuzzy)
chl_fuzzy <- resample(chl_fuzzy, sss_fuzzy)
spm_fuzzy <- resample(spm_fuzzy, sss_fuzzy)

# Compute the average suitablity score for each factor stack (until then we had one suitability score raster per year)
sss <- mean(sss_fuzzy, na.rm=T)
sst <- mean(sst_fuzzy, na.rm=T)
chl <- mean(chl_fuzzy, na.rm=T)
spm <- mean(spm_fuzzy, na.rm=T)

# Stack the result
SMCE_stack <- stack(sst, chl, sss, spm)

# Cmpute the suitability index (SI)
weights <- c(0.45,0.30,0.16,0.09)
SMCE_europe <- calc(SMCE_stack, fun = function(x) sum(x * weights))

writeRaster(SMCE_europe,"SMCE_full_europe_2019-2023.tif", overwrite=T)


# SMCE feasible Europe ----
# Now we apply the Feasiblity mask
feasibility <- raster("feasibility_mask.tif")
SMCE_europe <- raster("SMCE_full_europe_2019-2023.tif") # reload our full Europe SI if needed

# Rescale the feasibility mask on the SI
feasibility_res <- resample(feasibility, SMCE_europe, method="ngb")

# Duplicate the full Europe SI to crop it with the mask
SMCE_feasible <- SMCE_europe

# Mask out the full Europe SI with the feasibilty mask
SMCE_feasible[feasibility_res==0] <- NA

# Check your results
plot(SMCE_feasible)

writeRaster(SMCE_feasible,"SMCE_feasible_europe_2019-2023.tif", overwrite=T)


# SMCE future trend ----
# We compute the SI the same way we did previously but this time we use a different SST data (mean SST 2010-2020 and mean SST 2040-2050) 
## SMCE 2010-2020 ----
sst_present <- raster("SST_2010-2020_fuzzy.tif")

# Again, use SSS resolution
sst_present <- resample(sst_present, sss_fuzzy)

# STack the factors together and compute the SI
SMCE_stack_present <- stack(sst_present, chl, sss, spm)
weights <- c(0.45,0.30,0.16,0.09)
SMCE_europe_present <- calc(SMCE_stack_present, fun = function(x) sum(x * weights))

                            
## SMCE 2040-2050 ----
sst_future <- raster("SST_2040-2050_fuzzy.tif")
sst_future <- resample(sst_future, sss_fuzzy)

SMCE_stack_future <- stack(sst_future, chl, sss, spm)
weights <- c(0.45,0.30,0.16,0.09)
SMCE_europe_future <- calc(SMCE_stack_future, fun = function(x) sum(x * weights))

## Compute trend ----
SMCE_europe_trend <- 100*(SMCE_europe_future-SMCE_europe_present)/SMCE_europe_present

# Optional if your want the trend on the full extent of Europe
# writeRaster(SMCE_europe_trend,"Thresholds/SMCE_full_europe_absolute_trend.tif", overwrite=T)

# Generate an object to mask
SMCE_feasible_trend <- SMCE_europe_trend

# Mask the raster with the Feasibility mask
SMCE_feasible_trend[feasibility_res==0] <- NA

plot(SMCE_feasible_trend) # Check your results
writeRaster(SMCE_feasible_trend,"SMCE_feasible_europe_relative_trend.tif", overwrite=T)

