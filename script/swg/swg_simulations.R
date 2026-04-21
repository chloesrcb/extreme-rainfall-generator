library(sf)
library(dplyr)
library(ggplot2)
library(magick)
library(grid)
library(generain)

# Simulate a spatio-temporal r-Pareto field on a 5-minute grid
sim_rpareto_coords_m5 <- function(coords,
                                  steps,
                                  beta1, beta2, alpha1, alpha2,
                                  adv_m5 = c(0, 0),
                                  threshold = 1,
                                  s0_index = 1,
                                  t0_index = 1,
                                  seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  RandomFields::RFoptions(
    spConform = FALSE,
    allow_duplicated_locations = TRUE,
    install = "no"
  )

  # coords in Lambert-93 meters
  coords <- as.data.frame(coords)
  n_sites <- nrow(coords)
  n_steps <- length(steps)

  x_coords <- coords$x
  y_coords <- coords$y

  # Build full spatio-temporal grid
  step_grid <- expand.grid(
    site = seq_len(n_sites),
    it   = seq_len(n_steps)
  )

  x <- x_coords[step_grid$site]
  y <- y_coords[step_grid$site]
  step_t <- steps[step_grid$it]

  # Apply advection shift
  # adv_m5 is in meters per 5 min
  x_shift <- x - step_t * adv_m5[1]
  y_shift <- y - step_t * adv_m5[2]

  # get conditioning point (s0, t0)
  ind0 <- which(step_grid$site == s0_index & step_grid$it == t0_index)

  x0_shift <- x_shift[ind0]
  y0_shift <- y_shift[ind0]
  t0_step  <- step_t[ind0]

  # Define spatial and temporal fractional Brownian motion models
  model_space <- RandomFields::RMfbm(
    alpha = alpha1,
    var   = 2 * beta1,
    scale = 1
  )

  model_time <- RandomFields::RMfbm(
    alpha = alpha2,
    var   = 2 * beta2,
    scale = 1
  )

  # Compute gamma0(s, t) relative to the conditioning point
  # gamma0 = gamma_space + gamma_time
  gamma_space_vec <- RandomFields::RFvariogram(
    model_space,
    x = x_shift - x0_shift,
    y = y_shift - y0_shift
  )

  gamma_time_vec <- RandomFields::RFvariogram(
    model_time,
    x = step_t - t0_step
  )

  gamma0_vec <- gamma_space_vec + gamma_time_vec
  gamma0 <- matrix(gamma0_vec, nrow = n_sites, ncol = n_steps, byrow = FALSE)

  # Simulate spatial Gaussian component on unique shifted coordinates
  coord_key <- paste0(
    sprintf("%.8f", x_shift), "_",
    sprintf("%.8f", y_shift)
  )

  unique_keys <- unique(coord_key)
  map_idx <- match(coord_key, unique_keys)

  xy_unique <- do.call(rbind, strsplit(unique_keys, "_", fixed = TRUE))
  x_unique <- as.numeric(xy_unique[, 1])
  y_unique <- as.numeric(xy_unique[, 2])

  W_space_unique <- RandomFields::RFsimulate(
    model_space,
    x = x_unique,
    y = y_unique,
    grid = FALSE
  )

  W_space <- W_space_unique[map_idx]

  W_time <- RandomFields::RFsimulate(
    model_time,
    x = steps,
    grid = TRUE
  )
  # Combine spatial and temporal components
  W <- W_space + W_time[step_grid$it]
  W_mat <- matrix(W, nrow = n_sites, ncol = n_steps, byrow = FALSE)

  # r-Pareto transformation
  W0 <- W_mat[s0_index, t0_index]
  Y  <- exp(W_mat - W0 - gamma0)

  # simple Pareto
  R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)

  # Final field on standardized scale
  Z <- threshold * R * Y
  rownames(Z) <- rownames(coords)

  return(list(
    Z = Z,
    W = W_mat,
    gamma0 = gamma0,
    R = R
  ))
}

# Simulate one rainfall episode on the original rainfall scale
sim_episode_grid_m5 <- function(params_vario,
                                params_margins_common,
                                coords,
                                steps,
                                adv_m5,
                                t0,
                                s0_pixel_id,
                                u_emp) {

  # get conditioning pixel index
  s0_index <- which(rownames(coords) == s0_pixel_id)

  # Convert from rainfall scale to standard Pareto scale
  x_s0 <- pEGPD_full(
    u_emp,
    p0    = params_margins_common$p0,
    xi    = params_margins_common$xi,
    sigma = params_margins_common$sigma,
    kappa = params_margins_common$kappa
  )

  # G inverse CDF to standardize
  u <- G_std_inv(x_s0, p0 = params_margins_common$p0)

  # Simulate r-Pareto episode
  sim <- sim_rpareto_coords_m5(
    coords     = coords,
    steps      = steps,
    beta1      = params_vario$beta1,
    beta2      = params_vario$beta2,
    alpha1     = params_vario$alpha1,
    alpha2     = params_vario$alpha2,
    adv_m5     = adv_m5,
    s0_index   = s0_index,
    t0_index   = t0 + 1,
    threshold  = u
  )

  Z <- sim$Z
  n_sites <- nrow(coords)
  n_steps <- length(steps)

  # Transform back to rainfall scale
  X <- matrix(NA_real_, n_sites, n_steps, dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, n_sites, n_steps, dimnames = list(rownames(coords), NULL))

  for (k in seq_len(n_sites)) {
    Zk <- Z[k, ]

    V[k, ] <- G_std(Zk, p0 = params_margins_common$p0)

    X[k, ] <- qEGPD_full(
      V[k, ],
      p0    = params_margins_common$p0,
      xi    = params_margins_common$xi,
      sigma = params_margins_common$sigma,
      kappa = params_margins_common$kappa
    )
  }

  return(X)
}

s0_pixel_id <- "pixel_200"
u_emp <- 1 # 95% quantile of all observations ie = 1 mm/5min

# Number of 5 min time steps
nT <- 12
steps <- 0:(nT - 1)

# Example advection vector m/5min
adv_m5 <- c(50, -30)

# Convert from lon/lat to Lambert-93
grid_pts_l93 <- st_as_sf(
  grid_df,
  coords = c("Longitude", "Latitude"),
  crs = 4326
) %>%
  st_transform(2154)

grid_pts_xy <- st_coordinates(grid_pts_l93)

grid_df_l93 <- data.frame(
  x = grid_pts_xy[, 1],
  y = grid_pts_xy[, 2]
)
rownames(grid_df_l93) <- rownames(grid_df)

grid_poly_l93 <- st_transform(grid_latlon_poly, 2154)
sites_l93 <- st_transform(sites_sf, 2154)
s0_poly_l93 <- grid_poly_l93 %>%
  filter(pixel_id == s0_pixel_id)

# Simulate one episode on the rainfall scale
sim_episode <- sim_episode_grid_m5(
  params_vario           = params_m5min,
  params_margins_common  = params_margins_common,
  coords                 = grid_df_l93,
  steps                  = steps,
  adv_m5                 = adv_m5,
  t0                     = 0,
  s0_pixel_id            = s0_pixel_id,
  u_emp                  = u_emp
)

threshold <- u_emp
positive_vals <- sim_episode[is.finite(sim_episode) & sim_episode > 0]

fill_max <- if (length(positive_vals) > 0) {
  quantile(positive_vals, 0.995, na.rm = TRUE)
} else {
  threshold
}

fill_max <- max(fill_max, 2 * threshold)
fill_limits <- c(0, fill_max)

dir_frames <- file.path(
  im_folder,
  paste0("swg/omsev/frames_episode")
)

dir.create(dir_frames, recursive = TRUE, showWarnings = FALSE)

# polygon centroids in Lambert-93
grid_cent_l93 <- st_centroid(grid_poly_l93)
cent_xy <- st_coordinates(grid_cent_l93)

cent_df <- data.frame(
  pixel_id = grid_poly_l93$pixel_id,
  x = cent_xy[, 1],
  y = cent_xy[, 2]
)

dx_m <- adv_m5[1]
dy_m <- adv_m5[2]

# one frame per time step (5min)
for (tt in seq_len(nT)) {

  # Rainfall at current time step
  df_t <- data.frame(
    pixel_id = rownames(sim_episode),
    rain = sim_episode[, tt]
  )

  # add rainfall values to polygons
  map_t <- grid_poly_l93 %>%
    left_join(df_t, by = "pixel_id") %>%
    mutate(
      extreme   = rain > threshold,
      rain_plot = rain
    )

  # rainfall-weighted barycenter
  tmp <- df_t %>%
    left_join(cent_df, by = "pixel_id") %>%
    mutate(w = pmax(rain, 0))

  has_rain <- any(tmp$w > threshold / 4, na.rm = TRUE)

  # grid map
  p <- ggplot() +
    geom_sf(data = map_t, aes(fill = rain_plot), color = NA) +
    geom_sf(data = s0_poly_l93, fill = NA, color = "#ee8686", linewidth = 1.2) +
    geom_sf(data = sites_l93, color = "white", size = 3) +
    geom_sf(data = sites_l93, color = "grey40", size = 2) +
    scale_fill_gradientn(
      colours = c("white", "#dbe9f6", "#9ecae1", "#4a90d9", "#08519c"),
      values = scales::rescale(c(0, 0.2, threshold, 2 * threshold, fill_limits[2])),
      limits = fill_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = paste0("Rainfall (mm/5min)\nThreshold u = ", round(threshold, 2))
    ) +
    labs(
      title = paste0(
        "t = ", tt,
        "   |   advection per 5 min = (", dx_m, ", ", dy_m, ") m"
      )
    ) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(
      plot.title         = element_text(size = 18, hjust = 0.5),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      legend.background  = element_rect(fill = "white", color = NA),
      legend.key         = element_rect(fill = "white", color = NA),
      legend.text        = element_text(size = 14)
    )

  # rainfall-weighted barycenter and advection vector arrow
  if (has_rain) {

    bx <- sum(tmp$x * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)
    by <- sum(tmp$y * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)

    bary_pt_l93 <- st_sfc(st_point(c(bx, by)), crs = 2154) %>%
      st_sf()

    arrow_line_l93 <- st_sfc(
      st_linestring(rbind(
        c(bx, by),
        c(bx + dx_m, by + dy_m)
      )),
      crs = 2154
    ) %>%
      st_sf()

    p <- p +
      geom_sf(data = bary_pt_l93, color = "white", size = 2) +
      geom_sf(data = bary_pt_l93, color = "red", size = 1) +
      geom_sf(
        data = arrow_line_l93,
        arrow = arrow(length = unit(0.35, "cm")),
        color = "white",
        linewidth = 2
      ) +
      geom_sf(
        data = arrow_line_l93,
        arrow = arrow(length = unit(0.3, "cm")),
        color = "red",
        linewidth = 1
      )
  }

  # Save it
  out_png <- file.path(dir_frames, sprintf("frame_%03d.png", tt))

  ggsave(
    filename = out_png,
    plot = p,
    width = 7,
    height = 6,
    dpi = 150,
    bg = "transparent"
  )
}

# GIF
options(str = NULL)

files <- list.files(
  dir_frames,
  pattern = "frame_\\d+\\.png$",
  full.names = TRUE
)
files <- sort(files)

img <- image_read(files)
gif <- image_animate(img, fps = 1)

out_gif <- file.path(
  im_folder,
  paste0("swg/omsev/simulated_episode_", episode_idx, ".gif")
)

image_write(gif, path = out_gif)

out_gif