library(tidyverse)
library(spatial)

data_path <- "./Exercise 3/Data/"
fig_path <- "./Exercise 3/Figures/"

## Read data
# Seismic data observations
seismic <- scan(paste0(data_path, "seismic.dat"))
# Lithology distribution, with observation code {0, 1} for {sand, shale}
complit <- as.matrix(read.table(paste0(data_path, "complit.dat")))
dimnames(complit) <- NULL

