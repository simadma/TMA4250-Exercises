#Project 3 in TMA4250
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#PROBLEM 1

# Load libraries
library(tidyverse)
library(spatial)
library(fields)

data_path <- "./Exercise 3/Data/"
fig_path <- "./Exercise 3/Figures/"

## Read data
# Seismic data observations
seismic <- scan(paste0(data_path, "seismic.dat"))
# Lithology distribution, with observation code {0, 1} for {sand, shale}
complit <- as.matrix(read.table(paste0(data_path, "complit.dat")))
dimnames(complit) <- NULL
n_seis <- as.integer(sqrt(length(seismic)))

# mar_default <- c(5, 4, 4, 2) + 0.1
# Plot
# pdf(
#   file = paste0(fig_path, "obs.pdf"),
#   height = 4,
#   width = 5.2
# )
# par(mar = c(2, 2, 0, 2) + 0.5)
image.plot(x = 1:n_seis, y = 1:n_seis, z = matrix(seismic, nrow=n_seis),
           col=gray.colors(256), zlim = c(-0.08, 1.08), asp = 1, xlab = "", ylab = "")
# dev.off()


image(complit, col=gray.colors(256), asp=1)