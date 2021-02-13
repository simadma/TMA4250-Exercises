rm(list = ls())
library(geoR)    # Analysis of Geostatistical Data
library(akima)   # Interpolation of Irregularly and Regularly Spaced Data
library(fields)  # Tools for Spatial Data
library(MASS)    # Support Functions and Datasets for Venables and Ripley's MASS
library(tidyverse)

# Read data
data_path <- "./Data/"
fig_path <- "./Figures/"
data <- read.table(
          file = paste0(data_path, "topo.dat"),
          header = TRUE
        )

n <- 315

# Plot data
res <- akima::interp(
         x = data$x, y = data$y, z = data$z,
         xo = seq(0, n, by = 5), yo = seq(0, n, by = 5),
         linear = FALSE,
         extrap = TRUE
       )


# Plot
pdf(
  file = paste0(fig_path, "interp_obs.pdf"),
  width = 6,
  height = 5.75
)

fields::image.plot(res, col = terrain.colors(n = 256), asp = 1)
contour(res, nlevels = 11, add = TRUE)
points(x = data$x, y = data$y, pch = 4, col = "white")
points(x = data$x, y = data$y, pch = 3, col = "black")

dev.off()

?image.plot

viridis::ma