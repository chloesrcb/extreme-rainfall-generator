rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(animation)
library(tidyr)
library(lubridate)
library(purrr)
library(viridis)
library(magick)
library(sf)
library(units)
library(ggspatial)  

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

################################################################################
# DATA AND COORDS
################################################################################

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)

# egpd fit
filename_egpd <- paste0(data_folder, "../thesis/resources/images/EGPD/OMSEV/2019_2024/egpd_results.csv")
egpd_params <- read.csv(filename_egpd)

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]
location_gauges <- location_gauges[location_gauges$Station != "cines" &
                                   location_gauges$Station != "hydro" &
                                   location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)


sites_names <- colnames(rain)
sites_coords <- location_gauges[, c("Longitude", "Latitude")]

rownames(sites_coords) <- location_gauges$Station

sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)

sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000

# get distance matrix
grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
colnames(grid_coords_km) <- c("Longitude", "Latitude")


dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
df_dist_km <- df_dist
df_dist_km$value <- df_dist$value / 1000

rownames(grid_coords_m) <- rownames(sites_coords)


################################################################################
## Marginal parameters
################################################################################

# keep rain rows with at least 50% non NA
# rain_complete <- rain[rowSums(is.na(rain)) <= (ncol(rain) / 2), ]
# head(rain_complete)
# rain <- rain_complete

# to compute p0 remove all cumul over 1h above 10 mm
window_time <- 12 # 1 hour = 12 * 5min
# rain_cumul_1h <- zoo::rollapply(rain, width = window_time, FUN = sum, align = "right", fill = NA)
# rain_cumul_1h[rain_cumul_1h > 22] <- NA

p0_values <- sapply(colnames(rain), function(s) {
  mean(rain[, s] == 0, na.rm = TRUE)
})

kappa_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "kappa"]
xi_vect    <- egpd_params[match(colnames(rain), egpd_params$Site), "xi"]
sigma_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "sigma"]

params_margins <- list(
  xi    = xi_vect,
  sigma = sigma_vect,
  kappa = kappa_vect,
  p0    = p0_values
)

# for all sites names in params_margins$p0, put p0 to 0.95
# for (i in seq_along(params_margins$p0)) {
#   site <- names(params_margins$p0)[i]
#   if (site %in% colnames(rain)) {
#     params_margins$p0[i] <- 0.95
#   }
# }

#################################################################################
# GET EPISODES
#################################################################################
# rain_complete <- rain[rowSums(is.na(rain)) <= (ncol(rain) / 2), ]
# head(rain_complete)
# rain <- rain_complete
q <- 0.95
delta <- 12
dmin <- 1200  # m
times <- 0:(delta - 1)

# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
                            min_spatial_dist = dmin,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)
selected_points <- s0t0_set
selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

n_episodes <- length(selected_points$s0)
table(selected_points$s0) # number of episodes per s0

# do a plot for the number of episodes per site
ggplot(selected_points, aes(x = s0)) +
  geom_bar(fill = btfgreen, alpha = 0.7) +
  labs(
    x = "Site (s0)",
    y = "Number of episodes"
  ) +
  theme_minimal() 

# save plot
filename_plot <- paste0(
  im_folder,
  "swg/omsev/number_episodes_per_site_q",
  q * 100,
  "dmin",
  dmin,
  "delta",
  delta,
  ".png"
)
ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)

# number of episodes par t0 date
n_s0_per_t0 <- table(selected_points$t0_date)

t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0

# get extreme episodes
list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes
tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_m)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  u <- u_list[i]
  excesses <- empirical_excesses_rpar(episode, threshold = u,
                                      df_lags = lags, t0 = ind_t0_ep)
  # tau is in 5 minutes
  lags$tau <- lags$tau * 5 / 60 # convert to hours
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")
df_lags <- list_lags[[10]] # km/h
df_excesses <- list_excesses[[13]]
sum(df_excesses$kij)

length(list_episodes)

###################################################################################
# VARIOS PARAMETERS FROM KM/H TO M/5MIN
###################################################################################
# From results
params_est <- c(1.2953654, 4.2212009, 0.2495586, 0.6657619, 3.8962660, 2.2208320)
# params_est <-  c(1.2743853, 4.2106519, 0.2175129, 0.6678552, 5.3556250, 2.1357170)
# param_est <- c(0.9740506, 4.5719732, 0.2321974, 0.7175306, 1.6210000, 5.2190000)
etas_estimates <- params_est[5:6]

params_kmh <- list(
  beta1 = params_est[1],
  beta2 = params_est[2],
  alpha1 = params_est[3],
  alpha2 = params_est[4]
)

# convert params from km/h to m/5min
c_x_m <- 1000    # for m
c_t_5min <- 12   # 1 hour = 12 * 5min

# convert params and ci to m/5min
params_m5min_beta <- convert_params(params_kmh$beta1, params_kmh$beta2, params_kmh$alpha1, params_kmh$alpha2,
                               c_x = c_x_m, c_t = c_t_5min)

params_m5min <- list(
  beta1 = params_m5min_beta$beta1,
  beta2 = params_m5min_beta$beta2,
  alpha1 =  params_est[3],
  alpha2 = params_est[4]
)

beta1 <- params_m5min$beta1
beta2 <- params_m5min$beta2
alpha1 <- params_m5min$alpha1
alpha2 <- params_m5min$alpha2

##################################################################################
# TRANSFORM ADVECTION SPEEDS
##################################################################################
adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df_raw <- read.csv(adv_filename, sep = ",")
head(adv_df_raw)
setDT(selected_points)
setDT(adv_df_raw)

# is there na in adv_df_raw?
anyNA(adv_df_raw)

adv_df_raw[, t0_omsev := as.POSIXct(t0_omsev, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
selected_points[, t0_date := as.POSIXct(t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")]

adv_df_t0 <- adv_df_raw[, .(
  vx_final = vx_final[1],
  vy_final = vy_final[1]
), by = t0_omsev]

# is there na in adv_df_t0?
anyNA(adv_df_t0)
setkey(adv_df_t0, t0_omsev)

# get advections for selected points
selected_episodes <- adv_df_t0[selected_points, on = .(t0_omsev = t0_date)]
setnames(selected_episodes, c("vx_final","vy_final"), c("adv_x","adv_y"))
anyNA(selected_episodes)
which(is.na(selected_episodes$adv_x))
selected_episodes[is.na(selected_episodes$adv_x), ]
# Advection
V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)
adv_df <- V_episodes
colnames(adv_df) <- c("vx_final", "vy_final")
adv_df_transfo <- adv_df

# Transform advections with etas estimates
adv_df_transfo$vnorm <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
adv_df_transfo$vnorm_t <- etas_estimates[1] * adv_df_transfo$vnorm^etas_estimates[2]

adv_df_transfo$vx_t <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vx_final / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

adv_df_transfo$vy_t <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vy_final / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

plot(adv_df$vx_final, adv_df$vy_final, pch = 19,
     xlab = "vx", ylab = "vy")

plot(adv_df_transfo$vx_t, adv_df_transfo$vy_t, pch = 19,
     xlab = "vx", ylab = "vy")

# if too big advections after transformation, cap them
# max_adv <- 100  # km/h
# adv_df_transfo$vx_t <- pmax(pmin(adv_df_transfo$vx_t, max_adv), -max_adv)
# adv_df_transfo$vy_t <- pmax(pmin(adv_df_transfo$vy_t, max_adv), -max_adv)

# Compute speed
speed_raw <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
summary(speed_raw)
speed_transfo <- sqrt(adv_df_transfo$vx_t^2 + adv_df_transfo$vy_t^2)
summary(speed_transfo)

hist(speed_raw, breaks = 20, main = "Histogram of raw advection speeds", xlab = "Speed (km/h)")
hist(speed_transfo, breaks = 20, main = "Histogram of transformed advection speeds", xlab = "Speed (km/h)")
# Compute direction
direction_raw <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
direction_transfo <- atan2(adv_df_transfo$vy_t, adv_df_transfo$vx_t) * (180 / pi)



################################################################################
# Simulations over all episodes
################################################################################
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
summary(adv_matrix)
speed <- sqrt(adv_matrix[, 1]^2 + adv_matrix[, 2]^2)
summary(speed)
speed[speed>150]
grid_omsev <- grid_coords_m

# convert advection from km/h to m/5min
adv_matrix <- adv_matrix * (1000 * 5 / 60)

# remove "invalid"
id_not_invalid <- which(adv_class$group != "invalid")
adv_matrix_filtered <- adv_matrix[id_not_invalid, , drop = FALSE]
list_episodes_filtered <- list_episodes[id_not_invalid]
s0_list_filtered <- s0_list[id_not_invalid]
u_list_filtered <- u_list[id_not_invalid]

M <- 100 # repetitions per episode
Nsim <- length(list_episodes_filtered)  # number of episodes to simulate
sims_by_ep <- vector("list", Nsim)
u_sim <- numeric(Nsim)
timesteps <- 0:11 #* 5 / 60 # every 5 minutes for 1 hour
set.seed(123)
for (j in seq_len(Nsim)) {
  ep_idx <- j
  s0_j  <- s0_list_filtered[ep_idx]
  u_j   <- u_list_filtered[ep_idx]
  adv_j <- adv_matrix_filtered[ep_idx, ]
  u_sim[j] <- u_j
  sims_by_ep[[j]] <- vector("list", M)

  for (m in seq_len(M)) {
  sims_by_ep[[j]][[m]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_j,
      params_vario = params_m5min,
      params_margins = params_margins,
      coords = grid_omsev,
      times = timesteps,
      adv = adv_j,
      t0 = 0,
      s0 = s0_j
    )[[1]]
  }
}

N <- Nsim
length(sims_by_ep) * length(sims_by_ep[[1]])

# sims_by_ep[[i]][[m]] -> sims_all[[k]]
sims_all <- unlist(sims_by_ep, recursive = FALSE)

site_names <- colnames(sims_all[[1]]$X)


library(ggplot2)
library(dplyr)



# OBS cumul par épisode
Cobs <- sapply(seq_len(N), function(j) sum(as.matrix(list_episodes_filtered[[j]]), na.rm=TRUE))

# SIM cumul par épisode (N x M)
Csim <- matrix(NA_real_, nrow=N, ncol=M)
for (j in seq_len(N)) {
  for (m in seq_len(M)) {
    Csim[j,m] <- sum(sims_by_ep[[j]][[m]]$X, na.rm=TRUE)
  }
}

summary(Cobs)
summary(as.vector(Csim))

Cobs_pos <- Cobs
Csim_pos <- as.vector(Csim)
Csim_pos <- Csim_pos

# Violin sur log1p
library(ggplot2)
df <- rbind(
  data.frame(type="Observations", value=Cobs_pos),
  data.frame(type="Simulations", value=Csim_pos)
)
df$y <- log1p(df$value)

ggplot(df, aes(x=type, y=y)) +
  geom_violin(trim=FALSE, alpha=0.7, color=NA, scale = "width", bounds = c(0, Inf), bw = 0.5, fill = btfgreen) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white") +
  theme_minimal() +
  labs(y="log(1 + Cumul)", x="") +
  btf_theme


q05 <- apply(Csim, 1, quantile, probs=0.05, na.rm=TRUE)
q95 <- apply(Csim, 1, quantile, probs=0.95, na.rm=TRUE)
covered90 <- (Cobs >= q05) & (Cobs <= q95)

mean(covered90, na.rm=TRUE)




q <- quantile(Cobs, probs = c(0.33, 0.66), na.rm=TRUE)

class_obs <- cut(
  Cobs,
  breaks = c(-Inf, q[1], q[2], Inf),
  labels = c("weak", "medium", "strong")
)


# function to discretize low values in sim to match obs
# low obs = 0, p, 2*p, etc
apply_discretization <- function(x, p) {
  x[x < 1e-2] <- 0
  x[x > 0 & x < p] <- p
  x[x > p & x < 2*p] <- 2*p
  # x[x > 2*p & x < 3*p] <- 3*p
  x
}


apply_measurement <- function(x, p) {
  x[x < 1e-2] <- 0
  x[x > 0 & x < p] <- p
  x
}


sims_all <- unlist(sims_by_ep, recursive = FALSE)

site_names <- c("cefe", "iem", "um", "cnrs")
site_names <- colnames(sims_all[[1]]$X)

for (site in site_names) {

  Xsim_site <- unlist(lapply(sims_all, function(sim) as.numeric(sim$X[, site])))

  Xobs_site <- unlist(lapply(list_episodes_filtered, function(ep) as.numeric(ep[, site])))
  Xobs_site <- Xobs_site[!is.na(Xobs_site)]

  Xsim_site <- Xsim_site[Xsim_site > 1e-2]

  # Xsim_meas <- apply_discretization(Xsim_site, p = 0.2152)
  Xsim_meas <- Xsim_site
  Xobs_meas <- Xobs_site

  Xsim_site <- Xsim_meas[Xsim_meas > 0]
  Xobs_site <- Xobs_meas[Xobs_meas > 0]

  config <- "above0"

  df <- rbind(
    data.frame(type = "Observations", value = Xobs_site),
    data.frame(type = "Simulations", value = Xsim_site)
  )

  bw_common <- bw.nrd0(df$value)

  df_log <- df %>% mutate(yplot = log1p(value))

  psite_log <- ggplot(df_log, aes(x = type, y = yplot)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = NA,
                bw = 0.2, scale = "width", bounds= c(0, Inf), fill = btfgreen) +
    labs(x = "", y = "log(1 + Rain)") +
    theme_minimal() + btf_theme +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(legend.position = "none")

  folder_site <- paste0(im_folder, "swg/omsev/margins/", config, "/")
  if (!dir.exists(folder_site)) dir.create(folder_site, recursive = TRUE)

  ggsave(paste0(folder_site, "violin_log_", config, "_", site, ".png"),
         psite_log, width = 6, height = 4)

  psite_dens <- ggplot(df, aes(x = value, fill = type)) +
            geom_density(alpha = 0.4, bw = 0.2) +

            # boxplot en bas
            geom_boxplot(
              aes(y = -0.1, group = type),
              width = 0.15,
              alpha = 0.6,
              outlier.size = 0.8
            ) +

            labs(x = "Rainfall", y = "Density", fill = "Source") +
            theme_minimal() + btf_theme +
            xlim(0, 10)

  ggsave(paste0(folder_site, "density_box_", config, "_", site, ".png"),
         psite_dens, width = 8, height = 6)
}


# CUMULATIVE RAINFALL PER EPISODE AND SITE

sites_names <- colnames(rain)

# gestion of NA values
# in data if NA in obs, set sim to NA at the same episode and site
for (j in seq_len(Nsim)) {
  for (site in sites_names) {
    if (all(is.na(list_episodes_filtered[[j]][,site]))) {
      for (m in seq_along(sims_by_ep[[j]])) {
        sims_by_ep[[j]][[m]]$X[, site] <- NA
      }
    }
  }
}

# apply apply_measurement to all sims
for (j in seq_len(Nsim)) {
  for (m in seq_along(sims_by_ep[[j]])) {
    sims_by_ep[[j]][[m]]$X <- apply_discretization(sims_by_ep[[j]][[m]]$X, p = 0.2152)
  }
}

cumul_ep <- c()
for (j in seq_len(Nsim)) {
  ep_mat <- list_episodes_filtered[[j]]
  cumul_ep[j] <- sum(ep_mat, na.rm = TRUE)
}

df_cumul_obs <- data.frame(
  ep_id = seq_len(Nsim),
  cumul = cumul_ep
)

cumul_sim <- matrix(NA_real_, nrow=N, ncol=M)
mean_ep_sim <- numeric(Nsim)
for (j in seq_len(Nsim)) {
  cumul_j <- numeric(M)
  for (m in seq_len(M)) {
    cumul_sim[j,m] <- sum(sims_by_ep[[j]][[m]]$X, na.rm = TRUE)
    if(cumul_sim[j,m] < 1) {
      print(paste("Episode", j, "Rep", m))
      print(paste("Cumul sim:", cumul_sim[j,m]))   
     }
    cumul_j[m] <- cumul_sim[j,m]
  }
  mean_ep_sim[j] <- mean(cumul_j, na.rm = TRUE)
}

# [1] "Episode 237 Rep 27"
# [1] "Cumul sim: 0.865534598725329"
# ep_id <- 237
# rep <- 27
# sims_by_ep[[ep_id]][[rep]]$X

min_cum_sim <- min(mean_ep_sim, na.rm = TRUE)

df_cumul_sim <- as.data.frame(cumul_sim) %>%
  pivot_longer(cols = everything(), names_to = "rep", values_to = "cumul") %>%
  mutate(rep = as.integer(gsub("V", "", rep)))

df_cumul <- rbind(
  data.frame(type = "Observations", cumul = cumul_ep),
  data.frame(type = "Simulations", cumul = as.vector(cumul_sim))
)

# # # plot density
# ggplot(df_cumul_sim, aes(x = cumul)) +
#   geom_density(alpha = 0.4, bw = 20) +
#   labs(x = "Cumulative rainfall", y = "Density") +
#   btf_theme + 
#   xlim(0, 300)

# ggplot(df_cumul_obs, aes(x = cumul)) +
#   geom_density(alpha = 0.4, bw = 20) +
#   labs(x = "Cumulative rainfall", y = "Density") +
#   btf_theme + 
#   xlim(0, 300)

min_cum_sim <- min(df_cumul_sim$cumul, na.rm = TRUE)

df_cumul <- df_cumul %>% filter(is.finite(cumul))

bw_common <- bw.nrd(df_cumul$cumul)
ggplot(df_cumul, aes(x = cumul, fill = type)) +
  geom_density(alpha = 0.4, bw = 40) +
  labs(x = "Cumulative rainfall", y = "Density", fill = "Source") +
  geom_boxplot(
    aes(y = -0.0015, group = type),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  btf_theme + 
  coord_cartesian(xlim=c(0,700))

# save plot
folder_site <- paste0(im_folder, "swg/omsev/cumuls/")
filename <- paste0(folder_site, "cumul_density_all.png")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}
ggsave(filename, width = 10, height = 6)  


# Cumulatives observed and simulated
cum_obs_mat <- t(vapply(
  seq_along(list_episodes_filtered),
  function(i) sum_over_time_by_site(list_episodes_filtered[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) |>
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) |>
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

cum_sim_long <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    simX <- sims_by_ep[[j]][[m]]$X
    simX <- apply_measurement(simX, p = 0.2152)
    v <- sum_over_time_by_site(simX, sites_names)
    tibble(
      ep_id = j,
      rep   = m,
      site  = names(v),
      cum_sim = as.numeric(v)
    )
  })
})

selected_points_filtered <- selected_points[selected_points$speed_class != "invalid", ]

meta_ep <- tibble(
  ep_id = seq_len(length(list_episodes_filtered)),
  speed_class =   selected_points_filtered$speed_class
)

# Join OBS + SIM ep site
df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))
table(df_cum$speed_class)

# Cumuls observed and simulated per site
df_cum_site <- df_cum |>
  group_by(ep_id, site) |>
  summarise(
    cum_obs = first(cum_obs),
    cum_sim = mean(cum_sim, na.rm = FALSE),
    .groups = "drop"
  )


df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))
unique(df_cum$speed_class)
max(df_cum$cum_obs, na.rm=TRUE)
max(df_cum$cum_sim, na.rm=TRUE)
df_cum[is.na(df_cum$cum_obs), ]$cum_sim <- NA

head(sort(unique(df_cum$cum_obs)))
head(sort(unique(df_cum$cum_sim)))
summary(df_cum$cum_obs)
summary(df_cum$cum_sim)
obs_ep <- cum_obs_long |>
  group_by(ep_id) |>
  summarise(cumul_obs = sum(cum_obs, na.rm=TRUE), .groups="drop")

sim_eprep <- cum_sim_long |>
  group_by(ep_id, rep) |>
  summarise(cumul_sim = sum(cum_sim, na.rm=TRUE), .groups="drop")

df_all <- bind_rows(
  obs_ep |> transmute(ep_id, source="obs", cumul=cumul_obs),
  sim_eprep |> transmute(ep_id, source="sim", cumul=cumul_sim)
)

df_plot <- df_all |>
  left_join(meta_ep, by = "ep_id")
# density
ggplot(df_plot, aes(x = cumul, fill = source)) +
  geom_density(alpha = 0.4, bw = 40) +
  labs(x = "Episode cumul", y = "Density", fill = "Source") +
  btf_theme +
  geom_boxplot(
    aes(y = -0.0015, group = source),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  coord_cartesian(xlim=c(0,700))

# save plot
folder_site <- paste0(im_folder, "swg/omsev/cumuls/")
filename <- paste0(folder_site, "cumul_density_episode_all_disc.png")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}
ggsave(filename, width = 10, height = 6)

df_cumul <- df_cumul %>% filter(is.finite(cumul))



cumul_episode <- function(
  ep_mat,
  sites,
  s0,
  radius,
  time_idx = NULL,
  transfo = FALSE
) {

  x <- as.matrix(ep_mat)
  if (transfo) {
     x <- apply_measurement(x, p = 0.2152)
  }

  is_time_by_site <- !is.null(colnames(x)) && all(colnames(x) %in% sites)
  if (is_time_by_site) {
    x <- x[, sites, drop = FALSE]
    nT <- nrow(x)
  } else {
    x <- x[sites, , drop = FALSE]
    nT <- ncol(x)
  }

  if (!is.null(time_idx)) {
    time_idx2 <- time_idx[time_idx >= 1 & time_idx <= nT]
    if (length(time_idx2) == 0) return(NA_real_)

    if (is_time_by_site) {
      x <- x[time_idx2, , drop = FALSE]
    } else {
      x <- x[, time_idx2, drop = FALSE]
    }
  }

  sum(x, na.rm = TRUE)
}




# ============================================================================
# CUMULS OBSERVED VS SIMULATED WITH NA-MATCHING
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)

sites_names <- colnames(rain)

Nep <- length(list_episodes_filtered)
M   <- length(sims_by_ep[[1]])

for (j in seq_along(list_episodes_filtered)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])

  for (m in seq_along(sims_by_ep[[j]])) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])

    sim_mat[is.na(obs_mat)] <- NA

    sims_by_ep[[j]][[m]]$X[, sites_names] <- sim_mat
  }
}

for (j in seq_along(sims_by_ep)) {
  for (m in seq_along(sims_by_ep[[j]])) {
    sims_by_ep[[j]][[m]]$X <- apply_discretization(
      sims_by_ep[[j]][[m]]$X,
      p = 0.2152
    )
  }
}

cumul_obs_raw  <- numeric(Nep)
cumul_obs_norm <- numeric(Nep)
n_avail_obs    <- integer(Nep)

cumul_sim_raw  <- matrix(NA_real_, nrow = Nep, ncol = M)
cumul_sim_norm <- matrix(NA_real_, nrow = Nep, ncol = M)

for (j in seq_along(list_episodes_filtered)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])

  n_avail_obs[j] <- sum(!is.na(obs_mat))
  cumul_obs_raw[j] <- sum(obs_mat, na.rm = TRUE)

  if (n_avail_obs[j] > 0) {
    cumul_obs_norm[j] <- cumul_obs_raw[j] / n_avail_obs[j]
  } else {
    cumul_obs_norm[j] <- NA_real_
  }

  for (m in seq_along(sims_by_ep[[j]])) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])

    cumul_sim_raw[j, m] <- sum(sim_mat, na.rm = TRUE)

    if (n_avail_obs[j] > 0) {
      cumul_sim_norm[j, m] <- cumul_sim_raw[j, m] / n_avail_obs[j]
    } else {
      cumul_sim_norm[j, m] <- NA_real_
    }
  }
}

df_cumul_obs_raw <- tibble(
  ep_id = seq_len(Nep),
  cumul = cumul_obs_raw,
  source = "Observations"
)

df_cumul_obs_norm <- tibble(
  ep_id = seq_len(Nep),
  cumul = cumul_obs_norm,
  source = "Observations"
)

df_cumul_sim_raw <- as.data.frame(cumul_sim_raw) %>%
  mutate(ep_id = seq_len(n())) %>%
  pivot_longer(
    cols = -ep_id,
    names_to = "rep",
    values_to = "cumul"
  ) %>%
  mutate(
    rep = as.integer(gsub("V", "", rep)),
    source = "Simulations"
  )

df_cumul_sim_norm <- as.data.frame(cumul_sim_norm) %>%
  mutate(ep_id = seq_len(n())) %>%
  pivot_longer(
    cols = -ep_id,
    names_to = "rep",
    values_to = "cumul"
  ) %>%
  mutate(
    rep = as.integer(gsub("V", "", rep)),
    source = "Simulations"
  )

df_plot_raw <- bind_rows(
  df_cumul_obs_raw,
  df_cumul_sim_raw %>% select(ep_id, cumul, source)
) %>%
  filter(is.finite(cumul))

df_plot_norm <- bind_rows(
  df_cumul_obs_norm,
  df_cumul_sim_norm %>% select(ep_id, cumul, source)
) %>%
  filter(is.finite(cumul))

cat("Raw cumulative rainfall summaries:\n")
print(df_plot_raw %>% group_by(source) %>% summarise(
  n = n(),
  min = min(cumul, na.rm = TRUE),
  q25 = quantile(cumul, 0.25, na.rm = TRUE),
  median = median(cumul, na.rm = TRUE),
  mean = mean(cumul, na.rm = TRUE),
  q75 = quantile(cumul, 0.75, na.rm = TRUE),
  max = max(cumul, na.rm = TRUE)
))

cat("\nNormalized cumulative rainfall summaries:\n")
print(df_plot_norm %>% group_by(source) %>% summarise(
  n = n(),
  min = min(cumul, na.rm = TRUE),
  q25 = quantile(cumul, 0.25, na.rm = TRUE),
  median = median(cumul, na.rm = TRUE),
  mean = mean(cumul, na.rm = TRUE),
  q75 = quantile(cumul, 0.75, na.rm = TRUE),
  max = max(cumul, na.rm = TRUE)
))

p_raw <- ggplot(df_plot_raw, aes(x = cumul, fill = source)) +
  geom_density(alpha = 0.4, bw = 40) +
  labs(
    x = "Episode cumulative rainfall",
    y = "Density",
    fill = "Source"
  ) +
  geom_boxplot(
    aes(y = -0.0015, group = source),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  btf_theme +
  coord_cartesian(xlim = c(0, 700))

print(p_raw)

p_norm <- ggplot(df_plot_norm, aes(x = cumul, fill = source)) +
  geom_density(alpha = 0.4, bw = 0.8) +
  labs(
    x = "Normalized episode cumulative rainfall",
    y = "Density",
    fill = "Source"
  ) +
  geom_boxplot(
    aes(y = -0.05, group = source),
    width = 0.1,
    alpha = 0.5,
    outlier.size = 0.1
  ) +
  btf_theme +
  coord_cartesian(xlim = c(0, 10))

print(p_norm)


folder_site <- paste0(im_folder, "swg/omsev/cumuls/")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}

filename_raw <- paste0(folder_site, "cumul_density_episode_all_disc_masked.png")
ggsave(filename_raw, p_raw, width = 10, height = 6)

filename_norm <- paste0(folder_site, "cumul_density_episode_all_disc_masked_normalized.png")
ggsave(filename_norm, p_norm, width = 10, height = 6)

cat("Saved:\n")
cat(filename_raw, "\n")
cat(filename_norm, "\n")


mean_ep_sim_raw <- apply(cumul_sim_raw, 1, mean, na.rm = TRUE)
mean_ep_sim_norm <- apply(cumul_sim_norm, 1, mean, na.rm = TRUE)

df_compare_means <- tibble(
  ep_id = seq_len(Nep),
  obs_raw = cumul_obs_raw,
  sim_mean_raw = mean_ep_sim_raw,
  obs_norm = cumul_obs_norm,
  sim_mean_norm = mean_ep_sim_norm,
  n_avail_obs = n_avail_obs
)

head(df_compare_means)
summary(df_compare_means)

# ----------------------------------------------------------------------------
# 10) Optional scatterplots obs vs mean simulated
# ----------------------------------------------------------------------------
p_scatter_raw <- ggplot(df_compare_means, aes(x = obs_raw, y = sim_mean_raw)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  labs(
    x = "Observed episode cumulative rainfall",
    y = "Mean simulated episode cumulative rainfall"
  ) +
  btf_theme +
  coord_cartesian(xlim = c(0, 700), ylim = c(0, 700))

print(p_scatter_raw)

filename_scatter_raw <- paste0(folder_site, "cumul_scatter_obs_vs_sim_mean_raw.png")
ggsave(filename_scatter_raw, p_scatter_raw, width = 7, height = 6)

p_scatter_norm <- ggplot(df_compare_means, aes(x = obs_norm, y = sim_mean_norm)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  labs(
    x = "Observed normalized cumulative rainfall",
    y = "Mean simulated normalized cumulative rainfall"
  ) +
  btf_theme

print(p_scatter_norm)

filename_scatter_norm <- paste0(folder_site, "cumul_scatter_obs_vs_sim_mean_normalized.png")
ggsave(filename_scatter_norm, p_scatter_norm, width = 7, height = 6)

cat(filename_scatter_raw, "\n")
cat(filename_scatter_norm, "\n")


df_plot_raw %>%
  group_by(source) %>%
  summarise(
    q50 = quantile(cumul, 0.50, na.rm = TRUE),
    q75 = quantile(cumul, 0.75, na.rm = TRUE),
    q90 = quantile(cumul, 0.90, na.rm = TRUE),
    q95 = quantile(cumul, 0.95, na.rm = TRUE),
    q99 = quantile(cumul, 0.99, na.rm = TRUE),
    max = max(cumul, na.rm = TRUE)
  )


summary(df_compare_means$sim_mean_raw - df_compare_means$obs_raw)
cor(df_compare_means$obs_raw, df_compare_means$sim_mean_raw, use = "complete.obs")

df_compare_means %>%
  mutate(ratio = sim_mean_raw / obs_raw) %>%
  ggplot(aes(x = ratio)) +
  geom_histogram(bins = 50) +
  btf_theme

ggplot(df_compare_means, aes(x = obs_raw, y = sim_mean_raw)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
  btf_theme



################################################################################
# Plot function for split violin plots
################################################################################

# Split-violin geom
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, panel_scales, coord, draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- data
    newdata$x <- if (grp %% 2 == 1) data$xminv else data$xmaxv
    newdata <- newdata[order(newdata$y), ]
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- newdata[1, "x"]
    GeomPolygon$draw_panel(newdata, panel_scales, coord)
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., trim = TRUE, scale = "area",
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, na.rm = na.rm, ...)
  )
}


################################################################################
# Simulations
################################################################################
grid_omsev <- grid_coords_m
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
Nsim <- length(list_episodes)
s0_sim <- integer(Nsim)
sims_all <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  # idx <- sample(seq_along(list_episodes), 1)
  idx <- i
  s0_i <- s0_list[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list[idx]
  adv_i <- adv_matrix[idx, ]
  adv_sim[i, ] <- adv_i
  sims_all[[i]] <- simulate_many_episodes(
    N = 1,
    u_emp = u_i,
    params_vario = params_m5min,
    params_margins = params_margins,
    coords = grid_omsev, # km
    times = times,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}




Xsim_all <- do.call(cbind, lapply(sims_all, function(sim) sim$X))
# q <- c(0.1,0.25,0.5,0.75,0.9,0.95,0.99)
# apply_measurement <- function(x, p) {
#   x[x < 1e-4] <- 0
#   x[x > 0 & x < p] <- p
#   x
# }

# Qobs <- quantile(Xobs_site[Xobs_site>0], q, na.rm=TRUE)
# Qsim <- quantile(Xsim_site[Xsim_site>0.06], q, na.rm=TRUE)
# Xsim_meas <- apply_measurement(Xsim_site, p = 0.2152)
# Xsim_pos <-  Xsim_meas[Xsim_meas > 0.001]
# Qsim_m <- quantile(Xsim_pos, q, na.rm=TRUE)

# cbind(q=q, obs=Qobs, sim=Qsim, sim_measured=Qsim_m,
#       diff_before=Qsim-Qobs, diff_after=Qsim_m-Qobs)
apply_measurement <- function(x, p) {
  # x[x < 1e-4] <- 0
  x[x > 0 & x < p] <- p
  x
}

apply_disc <- function(x, p) {
  x[x < 1e-4] <- 0
  x[x > 0 & x < p] <- p
  x[x > p & x < 2*p] <- 2*p
  x
}

site_names <- colnames(sims_all[[1]]$X)

for (site in site_names) {
  # sim <- sims_all[[1]]
  Xsim_site <- unlist(lapply(sims_all, function(sim) as.numeric(sim$X[, site])))
  Xobs_site <- unlist(lapply(list_episodes, function(ep) as.numeric(ep[, site])))
  Xobs_site <- Xobs_site[!is.na(Xobs_site)]
  xr <- range(c(Xobs_site, Xsim_site), finite = TRUE)
  dens     <- density(Xobs_site, from = xr[1], to = xr[2], na.rm = TRUE)
  dens_sim <- density(Xsim_site, from = xr[1], to = xr[2], na.rm = TRUE)

  Xobs_site <- Xobs_site[Xobs_site > 0]
  p <- min(Xobs_site[Xobs_site > 0])
  # # put all values between 0 and 0.001 to 0 
  # Xsim_site[Xsim_site < 0.1] <- 0
  # and all values between 0.001 and 0.2 to 0.2
  # Xsim_site[Xsim_site >= 0.001 & Xsim_site < p] <- p
  Xsim_site <- Xsim_site[Xsim_site >  1e-3]

  Xsim_meas <- apply_disc(Xsim_site, p = 0.2152)
  Xsim_meas <- Xsim_site
  Xobs_meas <- Xobs_site

  # éventuellement retirer les 0 pour comparer les positives
  Xsim_site <- Xsim_meas[Xsim_meas > 0]
  Xobs_site <- Xobs_meas[Xobs_meas > 0]


  config <- "above0_simdisc_p02152"
  df <- rbind(
    data.frame(type = "Observed episodes", value = Xobs_site),
    data.frame(type = "Simulated episodes", value = Xsim_site)
  )
  bw_common <- bw.nrd0(df$value)

  # psite <- ggplot(df, aes(x = type, y = value, fill = type)) +
  #   geom_violin(alpha = 0.7, trim = F, scale = "width",
  #               bw = bw_common) +
  #   labs(x = "", y = "Rainfall (mm/5min)") +
  #   theme_minimal() + btf_theme +
  #   theme(legend.position = "none") 

  # # save plot
  # ggsave(paste0(im_folder, "swg/omsev/margins/violin_all", site, ".png"), psite, width = 6, height = 4)
  df_log <- df %>% mutate(yplot = log1p(value))
  # df_log <- df_log[df_log$yplot > 0,]
  psite_log <- ggplot(df_log, aes(x = type, y = yplot, fill = type)) +
    geom_violin(alpha = 0.7, trim = F, color = NA, bw = bw_common, scale = "width", bounds= c(0, Inf)) +
    labs(x = "", y = "log(1 + X)") +
    # scale_y_log10() +
    theme_minimal() + btf_theme +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(legend.position = "none")
  
  # save plot
  folder_site <- paste0(im_folder, "swg/omsev/margins/", config, "/")
  if (!dir.exists(folder_site)) {
    dir.create(folder_site, recursive = TRUE)
  }
  ggsave(paste0(folder_site, "violin_log_", config, "_", site, ".png"), psite_log, width = 6, height = 4)
}



# remove episode with more than 80% of active sites ie less than 20% of NA
list_episodes_filtered <- list_episodes[sapply(list_episodes, function(ep) mean(is.na(ep)) < 0.8)]
length(list_episodes_filtered) # 17 episodes left


cum_obs <- sapply(list_episodes_filtered, function(ep) sum(ep, na.rm = TRUE))
cum_sim <- sapply(sims_all, function(sim) sum(sim$X, na.rm = TRUE))

dens_cum_obs <- density(cum_obs, from = 0, to = max(cum_obs), na.rm = TRUE)
dens_cum_sim <- density(cum_sim, from = 0, to = max(cum_sim), na.rm = TRUE)
plot(dens_cum_obs, main = "Cumulative rainfall distribution", xlab = "Cumulative rainfall (mm)", ylab = "Density", col = "blue")
lines(dens_cum_sim, add = TRUE, col = "red")


dens_cum_obs$y <- dens_cum_obs$y / max(dens_cum_obs$y)
dens_cum_sim$y <- dens_cum_sim$y / max(dens_cum_sim$y)

plot(dens_cum_obs, col = "blue",
     main = "Normalized cumulative rainfall distribution",
     xlab = "Cumulative rainfall (mm)", ylab = "Relative density")
lines(dens_cum_sim, col = "red")
legend("topright", legend = c("Observed", "Simulated"),
       col = c("blue", "red"), lwd = 2)

# plot cumulative distribution
df_cum <- rbind(
  data.frame(type = "Observed episodes", cumul = cum_obs),
  data.frame(type = "Simulated episodes", cumul = cum_sim)
)

ggplot(df_cum, aes(x = type, y = cumul, fill = type)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  labs(x = "", y = "Total cumulative rainfall (mm)") +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")
