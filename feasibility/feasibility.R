#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Feasibility analysis %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(sf)
library(raster)
library(ncdf4)

setwd("C:/Users/my/path")

# Objectives: end up with a raster of 0 (not feasible) or 1 (feasible) for growing M. edulis
#   The output file could be a shapefile of feasible  areas too

#%%%% Main Criteria %%%%%%

# Bathymetry -----
# Thresholds chosen for the bathymetry is 5 to 100 m
bathy <- raster("bathymmetry.tif")

bathy_f <- bathy # create a bathy_feasible object
#plot(bathy_f)

bathy_f[bathy_f < -100] <- NA # NA to every pixels deeper than 100 m
bathy_f[bathy_f > -5] <- NA # NA to every pixels above 5 m in depth

bathy_f[!is.na(bathy_f)] <- 1 # 1 to pixels we keep
bathy_f[is.na(bathy_f)] <- 0 # 0 to pixels we mask out

save(bathy_f, file="bathy_feasibility.RData")
#writeRaster(bathy_f, filename="bathy_mask.tif", format="GTiff", overwrite=TRUE)


# Current speed -----
# Thresholds chosen for the current speed is 1 m/s

# Function to count how many time current sped is above threshold
# This is kind of overkill as we will mask out every pixel where the threshold is crossed at least once, but it's an interesting info to keep in mind 
extreme_current <- function(values, threshold, num_layers) {
  count_extreme_current <- sum(values >= threshold, na.rm = TRUE)
  return(count_extreme_current)
}


## Load uo (Eastward) and vo (Northward) vectors ---
uo <- raster::brick("cmems_mod_glo_phy-all_my_0.25deg_P1D-m_uo-vo_2019-2023.nc", varname = "uo_oras")
vo <- raster::brick("cmems_mod_glo_phy-all_my_0.25deg_P1D-m_uo-vo_2019-2023.nc", varname = "vo_oras")

## Compute current speed ---
current_speed <- sqrt(uo^2 + vo^2)

# Free some space if needed
# rm(uo, vo)

# Create our threshold
speed_max <- 1

# Count how many times the threshold is crossed for each pixel
r <- extreme_current(current_speed, speed_max, dim(current_speed)[[3]])
# plot(r)

# Create a current speed feasible object
current_f <- r

current_f[current_f>=1] <- -666    # Extreme negative value when the threshold is crossed  
current_f[current_f==0] <- 1       # 1 where the threshold has neven been crossed
current_f[current_f==-666] <- 0    # 0 where the threshold was crossed at least once
plot(current_f)

save(current_f, file="current_feasibility.RData")
#writeRaster(current_f, filename="current_speed_mask.tif", format="GTiff", overwrite=TRUE)


# Heat spike -----
# This function counts the total number of days within sequences of at least three consecutive days above Tmax (25°C in our case)
heat_spike_count <- function(arr, Tmax) {
  n_layers <- dim(arr)[3]  # Get the number of days we have to got through
  
  # Create an array to store the count of consecutive days
  count <- array(0, dim = dim(arr)[1:2])  # Initialise count to 0
  
  # Loop through each pixel
  for (i in 1:dim(arr)[1]) {
    for (j in 1:dim(arr)[2]) {
      consecutive_days <- 0  # Initialise the consecutive days counter
      
      for (k in 1:n_layers) {
        # If the pixel is above Tmax we increment consecutive_days by 1
        if (!is.na(arr[i, j, k]) && arr[i, j, k] >= Tmax){
          consecutive_days <- consecutive_days + 1  
        } else {
          # If the pixel is below the Tmax, we look at the consecutive_days value: 
          if (consecutive_days >= 3) { # If consecutive_days >= 3, that means at least the previous 3 days can be considered as a heat spike. We add all those days to the count. If not, we do nothing.
            count[i, j] <- count[i, j] + consecutive_days
          }
          consecutive_days <- 0  # Then, we reset the counter to 0
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

# Set a maximum temperature for heat spikes
Tmax <- 25 

# List all SST files in your folder. In my case, I had one file per year (otherwise the files were too big for a 5-year time serie. Later (for max wave height, I'll show you how to do it if you ahve more than one file per year)
SST_files <- list.files(path="my/path", pattern = "^SST_20") 

# Group the files per year (5 years in total)
file_groups <- split(SST_files, cut(seq_along(SST_files), 5, labels = FALSE)) 

heat_spike_stack <- c()
for(i in 1:length(file_groups)){
  
  print(paste0(i," out of ",length(file_groups))) # To keep track of what your code is up to
  
  # Open SST files per year
  SST_daily <- lapply(file_groups[[i]], function(filename) {
    raster::brick(paste0("my/path/", filename), varname = "analysed_sst")-272.15 # Open the data and convert kelvin in Celsius
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

# Create a raster stack with your heat spike counts
heat_spike_stack_count <- heat_spike_stack
plot(heat_spike_stack_count)

heat_spike_stack_count[heat_spike_stack_count != 0] <- 1     # 1 to keep pixels without any heat spike within the 5-year time serie
heat_spike_stack_count[heat_spike_stack_count == 0] <- 0     # 0 to mask out pixels with heat spike(s) within the 5-year time serie
heat_spike_stack_count <- sum(heat_spike_stack_count)        # count your heat spikes
plot(heat_spike_stack_count)

# Create a heat spike feasibility mask
heat_spike_f <- heat_spike_stack_count

heat_spike_f[heat_spike_f != 0] <- 1000    # Extreme value where heat spikes occurred 
heat_spike_f[heat_spike_f == 0] <- 1       # 1 when we keep piwels ithout heat spikes
heat_spike_f[heat_spike_f == 1000] <- 0    # 0 when we mask out pixels with heat spikes
plot(heat_spike_f)

save(heat_spike_f, file="heat_spike_feasibility.RData")
#writeRaster(heat_spike_f, filename="heat_spike_mask.tif", format="GTiff", overwrite=TRUE)



# Create feasibility mask ----
## Load feasible masks ----
load(file="heat_spike_feasibility.RData")    #heat_spike_f
load(file="bathy_feasibility.RData")         #bathy_f
load(file="current_feasibility.RData")       #current_f

# Check their resolution and extent, it's unlikely they have the same
res(heat_spike_f)
res(bathy_f)
res(current_f)

extent(heat_spike_f)
extent(bathy_f)
extent(current_f)

# Here, we want to use heat spike resolution
# Reascale the other masks on the heat spike mask
bathy_f <- resample(bathy_f,heat_spike_f, method="ngb")
current_f <- resample(current_f,heat_spike_f,method="ngb")

# Stack them up
maskstack <- stack(heat_spike_f,bathy_f,current_f)

# Sum them up. If the result of a pixel is 3 (or alternatively the number of masks you used), that mean this pixel can be considered 'feasible'
maskstacksum <- sum(maskstack)

maskstacksum[maskstacksum < dim(maskstack)[3]] <- 0    # 0 if non-feasible
maskstacksum[maskstacksum == dim(maskstack)[3]] <- 1   # 1 if 'feasible'
plot(maskstacksum)

writeRaster(maskstacksum, filename="feasibility_mask.tif", format="GTiff", overwrite=TRUE)


#%%%% Addtional Criteria %%%%%%


# Maximum wave height -----
# Add in supplementary material, not used for feasibility mask but more as a complementary data for special cases.

# Count how many times the wave height is above threshold (6 m in our case)
maximum_waves <- function(values, threshold, num_layers) {
  # Count the number of values above the threshold
  count_max_waves <- sum(values >= threshold, na.rm = TRUE)
  return(count_max_waves)
}

# List waves data files in our folder
waves_files <- list.files(path="my/path", pattern = "^waves_20") 
# Important note: I have one file per year for 2019, 2020 and 2021. I have 4 files per year for 2022 and 2023 (better satial resolution means bigger files so I had to split them)

# Patterns to split files into group
patterns <- c("waves_2019", "waves_2020", "waves_2021", "waves_2022", "waves_2023")

# Function to group raster names based on pattern
group_rasters <- function(pattern, raster_names) {
  grep(paste0("^", pattern), raster_names, value = TRUE)
}

# Group the .nc files by year with a pattern
file_groups <- lapply(patterns, group_rasters, waves_files)


# Load only the first layer of waves_2019 to resample 2022 and 2023 on a 0.2 resolution (currently 0.08, this is why hourly data for a year is heavy!)
resampler <- raster("waves_2019.nc")
ext_resampler <- extent(resampler)

# Create our threshold
Hs_max <- 6

# Loop through ecah year
wave_stack <- c()
for(i in 1:length(file_groups)){
  
  print(paste0(i," out of ",length(file_groups)))
  
  # Open a file per year
  Hs_hourly <- lapply(file_groups[[i]], function(filename) {
    raster::brick(paste0("MCE/CMEMS_data/", filename))
  })
  
  Hs_hourly <- stack(Hs_hourly) # 1 stack of the year's hourly data
  
  
  # Resample when needed
  if(extent(Hs_hourly) != ext_resampler){
    Hs_hourly <- resample(Hs_hourly,resampler)
  }

  # Compute the number of times max wave height is crosses
  r <- maximum_waves(Hs_hourly, Hs_max, dim(Hs_hourly)[[3]])
  
  names(r) <- substr(file_groups[i][[1]][1], 7, 10)
  
  # Loop to store output raster in a raster stack
  if(i==1){
    waves_stack <- r
  } else{
    waves_stack <- stack(waves_stack, r)
  }
}

# Count how many times the threshold is crossed (just out of curiousity)
waves_stack_sum <- sum(waves_stack, na.rm=TRUE)
plot(waves_stack_sum)

# Keep a binary mask for the GIS project
waves_stack_sum[waves_stack_sum>=1] <- 1

# Save your file as a tif file for a map
writeRaster(waves_stack_sum, filename="max_waves.tif", format="GTiff", overwrite=TRUE)



# Significant wave height -----
# Here we want to compute an "Accessibility" metric: relative time a site (pixel) is accessible. This means that we will compute over the tiome serie, when the piwel present a wave height below the threshold
# In the end, an Accessibility of 100% means that the site is always accesible (no rough sea), 0% means the site is never accesible as the sea is too rough (significant wave height always above max threshold).

# Function to calculate the percentage of values below the threshold for each pixel throughout the time serie
accessibility_time <- function(values, threshold, num_layers) {
  # Count the number of values below the threshold
  count_below_threshold <- sum(values < threshold, na.rm = TRUE)
  # Calculate the percentage
  percentage_below_threshold <- (count_below_threshold / num_layers) * 100 # Here 'layers' mean time steps (hours across the time serie)
  return(percentage_below_threshold)
}

# List the waves data files in the folder
waves_files <- list.files(path="MCE/CMEMS_data", pattern = "^waves_20")

# Patterns to split files into yearly groups
patterns <- c("waves_2019", "waves_2020", "waves_2021", "waves_2022", "waves_2023")

# Function to group raster names based on pattern
group_rasters <- function(pattern, raster_names) {
  grep(paste0("^", pattern), raster_names, value = TRUE)
}

# Group the .nc files by year with a pattern
file_groups <- lapply(patterns, group_rasters, waves_files)


# Load only the first layer of waves_2019 to resample 2022 and 2023 on a 0.2 resolution (currently 0.08)
resampler <- raster("MCE/CMEMS_data/waves_2019.nc")
ext_resampler <- extent(resampler)

# Create our threshold
Hs_max <- 1.5
wave_stack <- c()

for(i in 1:length(file_groups)){
  
  print(paste0(i," out of ",length(file_groups)))
 
   # Open files per year
  Hs_hourly <- lapply(file_groups[[i]], function(filename) {
    raster::brick(paste0("MCE/CMEMS_data/", filename))
  })
  
  Hs_hourly <- stack(Hs_hourly) # 1 stack of the year's hourly data

  # Resample if needed
  if(extent(Hs_hourly) != ext_resampler){
   Hs_hourly <- resample(Hs_hourly,resampler)
  }

  # Compute the % when Hs < Hs_max
  r <- accessibility_time(Hs_hourly, Hs_max, dim(Hs_hourly)[[3]])
  
  names(r) <- substr(file_groups[i][[1]][1], 7, 10)
  
  # Loop to store output raster in a raster stack
  if(i==1){
    waves_stack <- r
  } else{
    waves_stack <- stack(waves_stack, r)
  }
   
}

# Compute the mean of each year's Accessibility
waves_stack_mean <- mean(waves_stack, na.rm=TRUE)

writeRaster(waves_stack_mean, filename="accessibility.tiff", format="GTiff", overwrite=TRUE)










