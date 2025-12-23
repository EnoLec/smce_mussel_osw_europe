#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Feasibility area %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(sf)
library(raster)


setwd("C:/Users/enora/OneDrive - hull.ac.uk/001_Enora_PhD/Data/Suitable_WF_GIS_study")

# Objectives: end up with a raster of 0 (not feasible) or 1 (feasible) for growing M. edulis
#   The output file could be a shapefile of feasible  areas too


# Bathymetry -----
# Bathy needs to be between 5 and 100 m
bathy <- raster("MCE/GIS/GEBCO_European_elevation/GEBCO_euro_bathy.tif")

bathy_f <- bathy
plot(bathy_f)

bathy_f[bathy_f < -100] <- NA
bathy_f[bathy_f > -5] <- NA

bathy_f[!is.na(bathy_f)] <- 1
bathy_f[is.na(bathy_f)] <- 0


save(bathy_f, file="MCE/feasibility/bathy_feasibility.RData")
#writeRaster(bathy_f, filename="MCE/feasibility/bathy_mask.tif", format="GTiff", overwrite=TRUE)


# Current speed -----
# if current speed reaches 1.5 m/s we drop this location
extreme_current <- function(values, threshold, num_layers) {
  # Count the number of values above the threshold
  count_extreme_current <- sum(values >= threshold, na.rm = TRUE)
  return(count_extreme_current)
}

## Load uo and vo ---
uo <- raster::brick("MCE/CMEMS_data/cmems_mod_glo_phy-all_my_0.25deg_P1D-m_uo-vo_2019-2023.nc", varname = "uo_oras")
vo <- raster::brick("MCE/CMEMS_data/cmems_mod_glo_phy-all_my_0.25deg_P1D-m_uo-vo_2019-2023.nc", varname = "vo_oras")


## Compute current speed ---
current_speed <- sqrt(uo^2 + vo^2)


# Free some space
rm(uo, vo)

# Compute the how many time the maximum current speed is exceeded
speed_max <- 1
r <- extreme_current(current_speed, speed_max, dim(current_speed)[[3]])
r
plot(r)

current_f <- r
current_f[current_f>=1] <- -666
current_f[current_f==0] <- 1
current_f[current_f==-666] <- 0
plot(current_f)

save(current_f, file="MCE/feasibility/current_feasibility.RData")
#writeRaster(current_f, filename="MCE/feasibility/current_speed_mask.tif", format="GTiff", overwrite=TRUE)

# Maximum wave height -----
# Add in supplementary material as it can be an issue only for certain types of aquaculture lines (close to surface). But for
# long lines below 5m depth, wave height is not an issue.
# So not used for feasibility but more as a complementary data for special cases.
maximum_waves <- function(values, threshold, num_layers) {
  # Count the number of values above the threshold
  count_max_waves <- sum(values >= threshold, na.rm = TRUE)
  return(count_max_waves)
}

waves_files <- list.files(path="MCE/CMEMS_data", pattern = "^waves_20")

# Patterns to split files into group
patterns <- c("waves_2019", "waves_2020", "waves_2021", "waves_2022", "waves_2023")

# Function to group raster names based on pattern
group_rasters <- function(pattern, raster_names) {
  grep(paste0("^", pattern), raster_names, value = TRUE)
}

# Group the .nc files by year with a pattern
file_groups <- lapply(patterns, group_rasters, waves_files)


# LOad only the first layer of waves_2019 to resample 2022 and 2023 on a 0.2 resolution (currently 0.08)
resampler <- raster("MCE/CMEMS_data/waves_2019.nc")
ext_resampler <- extent(resampler)


Hs_max <- 6

wave_stack <- c()

for(i in 1:length(file_groups)){
  
  print(paste0(i," out of ",length(file_groups)))
  
  # Open SST files per year
  Hs_hourly <- lapply(file_groups[[i]], function(filename) {
    raster::brick(paste0("MCE/CMEMS_data/", filename))
  })
  
  Hs_hourly <- stack(Hs_hourly) # 1 stack of the year's daily data
  
  
  # Resample if needed (2023 only)
  if(extent(Hs_hourly) != ext_resampler){
    Hs_hourly <- resample(Hs_hourly,resampler)
  }
  
  
  # Compute the number of times max wave height is crosses
  r <- maximum_waves(Hs_hourly, Hs_max, dim(Hs_hourly)[[3]])
  
  ## Yearly average
  names(r) <- substr(file_groups[i][[1]][1], 7, 10)
  
  
  # Loop to store output raster in a raster stack
  if(i==1){
    waves_stack <- r
  } else{
    waves_stack <- stack(waves_stack, r)
  }
  
  
}


waves_stack_mean <- sum(waves_stack, na.rm=TRUE)

waves_stack_mean[waves_stack_mean>=1] <- 1

writeRaster(waves_stack_mean, filename="MCE/CMEMS_data/Thresholds/Hs/max_waves.tif", format="GTiff", overwrite=TRUE)




# Significant wave height -----
# Function to calculate the percentage of values below the threshold for each pixel
accessibility_time <- function(values, threshold, num_layers) {
  # Count the number of values below the threshold
  count_below_threshold <- sum(values < threshold, na.rm = TRUE)
  # Calculate the percentage
  percentage_below_threshold <- (count_below_threshold / num_layers) * 100
  return(percentage_below_threshold)
}


waves_files <- list.files(path="MCE/CMEMS_data", pattern = "^waves_20")

# Patterns to split files into group
patterns <- c("waves_2019", "waves_2020", "waves_2021", "waves_2022", "waves_2023")

# Function to group raster names based on pattern
group_rasters <- function(pattern, raster_names) {
  grep(paste0("^", pattern), raster_names, value = TRUE)
}

# Group the .nc files by year with a pattern
file_groups <- lapply(patterns, group_rasters, waves_files)


# LOad only the first layer of waves_2019 to resample 2022 and 2023 on a 0.2 resolution (currently 0.08)
resampler <- raster("MCE/CMEMS_data/waves_2019.nc")
ext_resampler <- extent(resampler)


Hs_max <- 1.5

wave_stack <- c()

for(i in 1:length(file_groups)){
 
  print(paste0(i," out of ",length(file_groups)))
 
   # Open SST files per year
  Hs_hourly <- lapply(file_groups[[i]], function(filename) {
    raster::brick(paste0("MCE/CMEMS_data/", filename))
  })
  
  Hs_hourly <- stack(Hs_hourly) # 1 stack of the year's daily data

  
  # Resample if needed (2023 only)
  if(extent(Hs_hourly) != ext_resampler){
   Hs_hourly <- resample(Hs_hourly,resampler)
  }


  # Compute the % when Hs < Hs_max
  r <- accessibility_time(Hs_hourly, Hs_max, dim(Hs_hourly)[[3]])
  
  ## Yearly average
  names(r) <- substr(file_groups[i][[1]][1], 7, 10)
  
  
  # Loop to store output raster in a raster stack
  if(i==1){
    waves_stack <- r
  } else{
    waves_stack <- stack(waves_stack, r)
  }
  
  
}


waves_stack_mean <- mean(waves_stack, na.rm=TRUE)



#save(waves_f, file="MCE/feasibility/significant_wave_feasibility.RData")
writeRaster(waves_stack_mean, filename="MCE/CMEMS_data/Thresholds/Hs/accessibility.tiff", format="GTiff", overwrite=TRUE)


# Heat spike -----
# This function counts the total number of days within sequences of at least three consecutive days above Tmax
heat_spike_count <- function(arr, Tmax) {
  n_layers <- dim(arr)[3]  # Get the number of days we have to got through
  
  # Create an array to store the count of consecutive days
  count <- array(0, dim = dim(arr)[1:2])  # Initialize all values to zero
  
  # Loop through each pixel
  for (i in 1:dim(arr)[1]) {
    for (j in 1:dim(arr)[2]) {
      consecutive_days <- 0  # Initialize the consecutive days counter
      
      for (k in 1:n_layers) {
        if (!is.na(arr[i, j, k]) && arr[i, j, k] >= Tmax){
          consecutive_days <- consecutive_days + 1  # If the pixel is above Tmax we increment consecutive_days
        } else {
          # If the pixel is below the Tmax, we look at the consecutive_days valu: 
          if (consecutive_days >= 3) { # I consecutive_days >= 3, that means at least the previous 3 days can be considered as a heat spike. We add all those days to the count
            count[i, j] <- count[i, j] + consecutive_days
          }
          consecutive_days <- 0  # Then, we reset the counter if the value is below the threshold
        }
      }
      
      # Final check after loop in case the sequence ends at the last day
      if (consecutive_days >= 3) {
        count[i, j] <- count[i, j] + consecutive_days
      }
    }
  }
  
  
  return(count = count)
}

# # Example array with the given temperatures for one pixel
# arr <- array(c(20, 22, 23, 25, NA, 20, 27, 28, 29, 21, 30, 31, 30, 29, 25, 24, 19, 15, 14, 14), dim = c(1, 1, 20))
# Tmax <- 24
# 
# result <- heat_spike_count(arr, Tmax)
# print(result)

Tmax <- 25 # Set a maximum temperature fro heat spikes

SST_files <- list.files(path="MCE/CMEMS_data", pattern = "^SST_20") 


file_groups <- split(SST_files, cut(seq_along(SST_files), 5, labels = FALSE)) # Group the files per year (5 years in total)

heat_spike_stack <- c()
for(i in 1:length(file_groups)){
  
  print(paste0(i," out of ",length(file_groups)))
  
  # Open SST files per year
  SST_daily <- lapply(file_groups[[i]], function(filename) {
    raster::brick(paste0("MCE/CMEMS_data/", filename), varname = "analysed_sst")-272.15 # Open the data and convert kelvin in Celsius
  })
  
  SST_daily <- stack(SST_daily)                          # 1 stack of the year's daily data
  r <- SST_daily@layers[[1]]                             # Need a raster file with the right dimension to overwrite with the count of heat spikes
  brick_array <- as.array(SST_daily)                     # Convert raster brick to an array for applying the function
  above_threshold <- heat_spike_count(brick_array, Tmax) # Count how many time 3 consecutive days are above 25 degrees
  values(r) <- above_threshold                           # overwrite the raster with the count of heat spikes for the year
  names(r) <- substr(names(r), 2, 5)                     # overwrite the name with the year
  
  # Loop to store output raster in a raster stack
  if(i==1){
    heat_spike_stack <- r
  } else{
    heat_spike_stack <- stack(heat_spike_stack, r)
  }
  
}


heat_spike_stack_count <- heat_spike_stack
plot(heat_spike_stack_count)

heat_spike_stack_count[heat_spike_stack_count != 0] <- 1
heat_spike_stack_count[heat_spike_stack_count == 0] <- 0
heat_spike_stack_count <- sum(heat_spike_stack_count)
plot(heat_spike_stack_count)

heat_spike_f <- heat_spike_stack_count
heat_spike_f[heat_spike_f != 0] <- 1000
heat_spike_f[heat_spike_f == 0] <- 1
heat_spike_f[heat_spike_f == 1000] <- 0
plot(heat_spike_f)

save(heat_spike_f, file="MCE/feasibility/heat_spike_feasibility.RData")
#writeRaster(heat_spike_f, filename="MCE/feasibility/heat_spike_mask.tif", format="GTiff", overwrite=TRUE)


# Create feasibility mask ----
## Load feasible masks ----
load(file="MCE/feasibility/heat_spike_feasibility.RData") #heat_spike_f
load(file="MCE/feasibility/bathy_feasibility.RData") #bathy_f
load(file="MCE/feasibility/current_feasibility.RData") #current_f
#load(file="MCE/feasibility/significant_wave_feasibility.RData") #waves_f

res(heat_spike_f)
res(bathy_f)
res(current_f)
#res(waves_f)

extent(heat_spike_f)
extent(bathy_f)
extent(current_f)


bathy_f <- resample(bathy_f,heat_spike_f, method="ngb")

current_f <- resample(current_f,heat_spike_f,method="ngb")




test <- stack(heat_spike_f,bathy_f,current_f)
plot(test)

testsum <- sum(test)
testmsk <- testsum

testmsk[testmsk < dim(test)[3]] <- 0
testmsk[testmsk == dim(test)[3]] <- 1
plot(testmsk)



writeRaster(testmsk, filename="MCE/feasibility/feasibility_mask.tif", format="GTiff", overwrite=TRUE)



