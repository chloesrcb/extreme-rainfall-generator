# Load required libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
source("./script/load_libraries.R")

# OMSEV rainfall 
filename_omsev <- file.path(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain <- read.csv(filename_omsev)
rain$dates <- as.POSIXct(rain$dates, tz = "UTC")

# OMSEV gauge locations 
filename_loc <- file.path(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem","mse","poly","um","cefe","cnrs","crbm",
                             "archiw","archie","um35","chu1","chu2","chu3",
                             "chu4","chu5","chu6","chu7","cines","brives","hydro")

# COMEPHORE rainfall and pixel locations
filename_com <- file.path(data_folder, "comephore/comephore_2008_2024_2km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
comephore_raw$date <- as.POSIXct(comephore_raw$date, tz = "UTC")

filename_loc_px <- file.path(data_folder, "comephore/loc_px_zoom_2km.csv")
loc_px <- read.csv(filename_loc_px, sep = ",")

# Reshape COMEPHORE data to long format and join with pixel coordinates
com_long <- comephore_raw %>%
  pivot_longer(-date, names_to = "Station", values_to = "rain") %>%
  rename(dates = date) %>%
  left_join(loc_px, by = c("Station" = "pixel_name"))

# Function to compute advection for a single episode
# t0: episode time (POSIXct)
# pre_episode: time before episode (difftime)
# post_episode: time after episode (difftime)
# data: data frame with columns 'dates', 'Longitude', 'Latitude', 'rain'
# rain_threshold: minimum rain to consider for barycentre computation
# meters_per_degree: conversion factor from degrees to meters
# cos_lat: cosine of mean latitude for longitudinal distance correction
compute_episode_advection <- function(t0, pre_episode, post_episode, data,
                                      rain_threshold, meters_per_degree, cos_lat) {
  # Extract time window around the episode
  t_start <- t0 - pre_episode
  t_end <- t0 + post_episode

  # Select all unique time steps in this window
  times_ep <- data %>%
    filter(dates >= t_start & dates <= t_end) %>%
    pull(dates) %>%
    unique() %>%
    sort()

  centers <- list()
  valid_times <- c()

  # Compute rain-weighted barycentre at each time step
  for (t in times_ep) {
    subset <- data %>% filter(dates == t)
    rain <- subset$rain
    lon <- subset$Longitude
    lat <- subset$Latitude
    
    # Skip time steps with insufficient rain
    if (length(rain) == 0 || max(rain, na.rm = TRUE) < rain_threshold) next
    
    # Compute barycentre coordinates
    x_c <- sum(rain * lon, na.rm = TRUE) / sum(rain, na.rm = TRUE)
    y_c <- sum(rain * lat, na.rm = TRUE) / sum(rain, na.rm = TRUE)
    if (is.na(x_c) || is.na(y_c)) next
    
    # Store the center and time
    centers[[length(centers) + 1]] <- c(x_c, y_c)
    valid_times <- c(valid_times, t)
  }

  # If less than 2 time steps, return zeros
  if (length(valid_times) < 2) {
    return(list(mean_dx_kmh = 0, mean_dy_kmh = 0))
  }

  centers_mat <- do.call(rbind, centers)
  duration_hours <- as.numeric(difftime(max(valid_times), min(valid_times), units = "hours"))

  # Compute displacement in km/h
  total_disp <- centers_mat[nrow(centers_mat), ] - centers_mat[1, ]
  dx_kmh <- (total_disp[1] * meters_per_degree * cos_lat) / 1000 / duration_hours
  dy_kmh <- (total_disp[2] * meters_per_degree) / 1000 / duration_hours

  list(mean_dx_kmh = dx_kmh, mean_dy_kmh = dy_kmh)
}


# Compute advections for different configurations
q_values <- c(95)
delta_values <- c(12)
min_spatial_dist_values <- c(1600)
beta_lag <- 0
if (beta_lag == 0) {
end_file <- ""
} else {
end_file <- paste0("_beta", beta_lag)
}

# Reshape rainfall data to long format
rain_long <- rain %>%
  pivot_longer(-dates, names_to = "Station", values_to = "rain") %>%
  left_join(location_gauges, by = "Station")

for (quantile in q_values) {
  for (delta in delta_values) {
    for (dmin in min_spatial_dist_values) {

      # Load episode start times (t0)
      filename_t0 <- sprintf("t0_5min_episodes_q%d_delta%d_dmin%d", quantile, delta, dmin)
      filename_t0 <- paste0(filename_t0, end_file, ".csv")
      t0_df <- read_csv(file.path(data_folder, "omsev", filename_t0), show_col_types = FALSE)
      colnames(t0_df)[1] <- "t0"
      t0_df$t0 <- as.POSIXct(ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}$", t0_df$t0),
                                    paste0(t0_df$t0, " 00:00:00"), t0_df$t0),
                             tz = "UTC")
      # Get minimum rain threshold above zero
      rain_threshold <- as.numeric(min(rain[rain > 0], na.rm = TRUE))
      # Constants for distance conversion
      meters_per_degree <- 111320
      cos_lat_omsev <- cos(mean(location_gauges$Latitude) * pi/180)
      cos_lat_com <- cos(mean(loc_px$Latitude) * pi/180)
      pre_episode_omsev <- minutes(30) # time before episode
      post_episode_omsev <- minutes(30) # time after episode
      pre_episode_com <- hours(1) # time before episode
      post_episode_com <- hours(1) # time after episode

      # OMSEV barycentric advection
      episode_advect_omsev <- lapply(seq_len(nrow(t0_df)), function(i) {
        t0 <- t0_df$t0[i]
        res <- compute_episode_advection(
          t0, pre_episode_omsev, post_episode_omsev,
          data = rain_long,
          rain_threshold = rain_threshold,
          meters_per_degree = meters_per_degree,
          cos_lat = cos_lat_omsev
        )
        data.frame(episode = i, t0_omsev = t0,
                   mean_dx_kmh_omsev = res$mean_dx_kmh,
                   mean_dy_kmh_omsev = res$mean_dy_kmh)
      }) %>% bind_rows()

      # COMEPHORE gridded advection
      episode_advect_com <- lapply(seq_len(nrow(t0_df)), function(i) {
        t0 <- ceiling_date(t0_df$t0[i], unit = "hour")  # align to hourly COMEPHORE timestep
        res <- compute_episode_advection(
          t0, pre_episode_com, post_episode_com,
          data = com_long,
          rain_threshold = 0.1,
          meters_per_degree = meters_per_degree,
          cos_lat = cos_lat_com
        )
        data.frame(episode = i, t0_comephore = t0,
                   mean_dx_kmh_comephore = res$mean_dx_kmh,
                   mean_dy_kmh_comephore = res$mean_dy_kmh)
      }) %>% bind_rows()

      # Combine the two sources
      episode_advect_df <- left_join(episode_advect_omsev, episode_advect_com, by = "episode")
      episode_advect_df$vx_final <- ifelse(episode_advect_df$mean_dx_kmh_comephore == 0
                                                  & episode_advect_df$mean_dy_kmh_comephore == 0,
                                                  episode_advect_df$mean_dx_kmh_omsev,
                                                  episode_advect_df$mean_dx_kmh_comephore)
      episode_advect_df$vy_final <- ifelse(episode_advect_df$mean_dx_kmh_comephore == 0
                                                  & episode_advect_df$mean_dy_kmh_comephore == 0,
                                                  episode_advect_df$mean_dy_kmh_omsev,
                                                  episode_advect_df$mean_dy_kmh_comephore)
      # if comephore advection is 0 put omsev advection
      # OMSEV, COMEPHORE, and COMBINED advections
      out_dir <- file.path(data_folder, "omsev/adv_estim")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      filename_out <- sprintf("episode_advection_q%d_delta%d_dmin%d",
                              quantile, delta, dmin)
      filename_out <- paste0(filename_out, end_file, ".csv")
      out_dir_omsev <- file.path(out_dir, "bary_omsev")
      write.csv(episode_advect_omsev,
                file.path(out_dir_omsev, filename_out), row.names = FALSE)

      out_dir_com <- file.path(out_dir, "bary_comephore")
      if (!dir.exists(out_dir_com)) dir.create(out_dir_com, recursive = TRUE)
      write.csv(episode_advect_com,
                file.path(out_dir_com, filename_out), row.names = FALSE)
      out_dir_comb <- file.path(out_dir, "combined_comephore_omsev")
      if (!dir.exists(out_dir_comb)) dir.create(out_dir_comb, recursive = TRUE)
      write.csv(episode_advect_df,
                file.path(out_dir_comb, filename_out), row.names = FALSE)

      cat("Saved OMSEV, COMEPHORE, and combined advections for q =",
          quantile, "delta =", delta, "dmin =", dmin, "\n")
    }
  }
}