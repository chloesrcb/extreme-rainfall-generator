adv <- c(0.5, 0.2) # advection
params <- c(0.4, 0.8, 0.2, 0.7)
ngrid <- 7
temp <- 0:11
s0 <- c(1, 1)
t0 <- 0
random_s0 <- TRUE # random initial location
s0_radius <- Inf # radius between episodes
M <- 50 # number of simulations
m <- 500 # number of extreme episodes
distance_type <- "euclidean" # "euclidean" or "lalpha"
eta1 <- 0.5 # transformation parameter for advection
eta2 <- 1.5
fixed_eta1 <- eta1
fixed_eta2 <- eta2
use_wind_data <- TRUE
