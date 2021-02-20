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

lin <- with(data, interp(x, y, z, nx = 70, ny = 70))

# Set-up for color in `persp()` plot
nbcol <- 256
zfacet <- with(lin, z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])
facetcol <- cut(zfacet, nbcol)


# Plot
# pdf(
#   file = paste0(fig_path, "interp_obs.pdf"),
#   height = 4,
#   width = 7
# )

# Plot
mar_default <- c(5, 4, 4, 2) + 0.1
# par(mfrow = c(1, 2))

image.plot(lin, col = terrain.colors(nbcol), asp = 1)
contour(lin, nlevels = 11, add = TRUE)
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))
persp(lin, scale = FALSE,
  theta = 135, phi = 30,
  col = terrain.colors(nbcol)[facetcol],
  lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8),
  xlab = "X", zlab = "Elevation",
  xaxs = "i", yaxs = "i"
)


# dev.off()

## c)
m = 70
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
  x    = gridline,
  y    = gridline,
  z    = matrix(ord_kc$predict, nrow = m),
  var  = matrix(ord_kc$krige.var, nrow = m),
  beta = matrix(ord_kc$beta.est, nrow = m, ncol = m)
)

# Set-up for color in `persp()` plot
zfacet <- with(ok, z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])
facetcol <- cut(zfacet, nbcol)

## Plot
# pdf(
#   file = paste0(fig_path, "ord_krig.pdf"),
#   height = 1.4 * 4,
#   width = 7
# )

# par(mfrow = c(2, 2), mar = mar_default - c(3, 0, 3, 0))

# Kriging prediction
image.plot(ok, col = terrain.colors(nbcol), asp = 1)
contour(ok, nlevels = 11, add = TRUE)
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

persp(ok,
  theta = 135, phi = 30,
  col = terrain.colors(nbcol)[facetcol],
  lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8),
  scale = FALSE,
  xlab = "X", zlab = "Elevation",
  xaxs = "i", yaxs = "i"
)

# Prediction standard deviation
with(ok, image.plot(x, y, sqrt(var), col = inferno(nbcol), asp = 1, ylab = ""))
with(ok, contour(x, y, sqrt(var), nlevels = 5, add = TRUE))
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

with(ok, persp(x, y, sqrt(var),
           theta = 135, phi = 30,
           col = inferno(nbcol)[facetcol],
           lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8),
           scale = FALSE,
           xlab = "X", zlab = "Elevation",
           xaxs = "i", yaxs = "i"
         )
)

# dev.off()

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
  var = matrix(uni_kc$krige.var, nrow = m),
  zprior = outer(
    X   = gridline,
    Y   = gridline,
    FUN = function(X, Y) {
      uni_kc$beta.est[1] + X * uni_kc$beta.est[2] + Y * uni_kc$beta.est[3] +
        X * Y * uni_kc$beta.est[4] + X^2 * uni_kc$beta.est[5] + Y^2 * uni_kc$beta.est[6]
    }
  )
)

# Set-up for color in `persp()` plot
zfacet <- with(uk, z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])
facetcol <- cut(zfacet, nbcol)

## Plot
# pdf(
#   file = paste0(fig_path, "uni_krig.pdf"),
#   height = 2.1 * 4,
#   width = 7
# )

# par(mfrow = c(3, 2), mar = mar_default - c(3, 0, 3, 0))

# Expected prior
with(uk, image.plot(x, y, zprior, col = terrain.colors(nbcol), asp = 1, ylab = ""))
with(uk, contour(x, y, zprior, nlevels = 10, add = TRUE))

with(uk, persp(x, y, zprior,
           theta = 135, phi = 30,
           col = terrain.colors(nbcol)[facetcol],
           lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8),
           scale = FALSE,
           xlab = "X", zlab = "Elevation",
           xaxs = "i", yaxs = "i"
         )
)

# Kriging prediction
image.plot(uk, col = terrain.colors(nbcol), asp = 1)
contour(uk, nlevels = 11, add = TRUE)
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

persp(uk,
      theta = 135, phi = 30,
      col = terrain.colors(nbcol)[facetcol],
      lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8),
      scale = FALSE,
      xlab = "X", zlab = "Elevation",
      xaxs = "i", yaxs = "i"
)

# Prediction standard deviation
with(uk, image.plot(x, y, sqrt(var), col = inferno(nbcol), asp = 1, ylab = ""))
with(uk, contour(x, y, sqrt(var), nlevels = 5, add = TRUE))
with(data, points(x, y, pch = 4, col = "white"))
with(data, points(x, y, pch = 3, col = "black"))

with(uk, persp(x, y, sqrt(var),
           theta = 135, phi = 30,
           col = inferno(nbcol)[facetcol],
           lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8),
           scale = FALSE,
           xlab = "X", zlab = "Elevation",
           xaxs = "i", yaxs = "i"
         )
)

# dev.off()

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