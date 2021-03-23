#Problem 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-23-21

#libraries
library(tidyverse)
library(spatial)

#################################################
#           PROBLEM 4                           #
#################################################

#paths
data_path <- "./Exercise 2/Data/"
fig_path <- "./Exercise 2/Figures/"

## Figure size
w <- 10  # cm
h <- 6  # cm


#Get data
# Biological Cells Point Pattern
cells <- read_table(paste0(data_path, "cells.dat"), col_names = c("x", "y"))

##########################################################################################
#                                 Subproblem a)                                          #
##########################################################################################






