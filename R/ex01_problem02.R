rm(list = ls())
library(geoR)     # Analysis of Geostatistical Data
library(akima)    # Interpolation of Irregularly and Regularly Spaced Data
library(fields)   # Tools for Spatial Data
library(MASS)     # Support Functions and Datasets for Venables and Ripley's MASS
library(viridis)  # Default Color Maps from 'matplotlib'
library(tidyverse)

# Read data
data_path <- "./Data/"
fig_path <- "./Figures/"
data <- read.table(file = paste0(data_path, "topo.dat"), header = TRUE)

n <- 315
spl <- with(data,
            interp(x, y, z,
                   xo     = seq(0, n, by = 5),
                   yo     = seq(0, n, by = 5),
                   linear = FALSE,
                   extrap = TRUE
            )
)

# Plot
# pdf(
#   file = paste0(fig_path, "interp_obs.pdf"),
#   width = 6,
#   height = 5.75
# )
image.plot(spl, col = terrain.colors(n = 256), asp = 1)
contour(spl, nlevels = 11, add = TRUE)
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

# dev.off()

## c)
m = 100
gridline <- seq(1, n, length.out = m)
locs <- expand.grid(x = gridline, y = gridline)
geodata <- as.geodata(data)
ord_kc <- krige.conv(
  geodata   = geodata,
  locations = locs,
  krige     = krige.control(
    type.krige = "ok",       # Ordinary Kriging
    trend.d    = "cte",      # Constant mean of data values
    trend.l    = "cte",      # Constant mean at prediction locations
    cov.model  = "powered.exponential",
    cov.pars   = c(2500, 100),  # (sigma_r^2, xi)
    kappa      = 1.5            # nu
  )
)

ok <- list(
  x   = gridline,
  y   = gridline,
  z   = matrix(ord_kc$predict, nrow = m),
  var = matrix(ord_kc$krige.var, nrow = m)
)
image.plot(ok, col = terrain.colors(n = 256), asp = 1)
contour(ok, nlevels = 11, add = TRUE)
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

with(ok, image.plot(x, y, var, col = inferno(n = 256), asp = 1))
with(ok, contour(x, y, var, nlevels = 10, add = TRUE))
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))


## d)
uni_kc <- krige.conv(
  geodata   = geodata,
  locations = locs,
  krige     = krige.control(
    type.krige = "ok",  # Universal Kriging, specified with trend below
    trend.d = trend.spatial(~ x + y + I(x * y) + I(x^2) + I(y^2), geodata),
    trend.l = with(locs, trend.spatial(~ x + y + I(x * y) + I(x^2) + I(y^2))),
    cov.model  = "powered.exponential",
    cov.pars   = c(2500, 100),  # (sigma_r^2, xi)
    kappa      = 1.5            # nu
  )
)

uk <- list(
  x   = gridline,
  y   = gridline,
  z   = matrix(uni_kc$predict, nrow = m),
  var = matrix(uni_kc$krige.var, nrow = m)
)

image.plot(uk, col = terrain.colors(n = 256), asp = 1)
contour(uk, nlevels = 11, add = TRUE)
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

with(uk, image.plot(x, y, var, col = inferno(n = 256), asp = 1))
with(uk, contour(x, y, var, nlevels = 10, add = TRUE))
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))



## e)
ord_kc2 <- krige.conv(
  geodata   = geodata,
  locations = expand.grid(1:n, 1:n),
  krige     = krige.control(
    type.krige = "ok",       # Ordinary Kriging
    trend.d    = "cte",      # Constant mean of data values
    trend.l    = "cte",      # Constant mean at prediction locations
    cov.model  = "powered.exponential",
    cov.pars   = c(2500, 100),  # (sigma_r^2, xi)
    kappa      = 1.5            # nu
  )
)

ok2 <- list(
  x   = gridline,
  y   = gridline,
  z   = matrix(ord_kc2$predict, nrow = n),
  var = matrix(ord_kc2$krige.var, nrow = n)
)
x0 <- 100
y0 <- 100
r0_hat <- ok2$z[x0, y0]          # Kriging predictor at (x0, y0) = (100, 100)
sigma_r_hat2 <- ok2$var[x0, y0]  # Prediction variance at (x0, y0) = (100, 100)

pnorm(850, mean = r0_hat, sd = sqrt(sigma_r_hat2), lower.tail = FALSE)  # Pr(r0 > 850 | d)
qnorm(0.9, mean = r0_hat, sd = sqrt(sigma_r_hat2))  # z such that Pr(r0 < z | d) = 0.90
