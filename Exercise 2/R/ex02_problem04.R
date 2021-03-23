#Exercise 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-23-21

#libraries
source("./Exercise 2/R/utils.R")
##########################################################################################
#                                     PROBLEM 4                                          #
##########################################################################################


## Figure size
w <- 10  # cm
h <- 6  # cm


## Read data
# Biological Cells Point Pattern
cells <- read_table(paste0(data_path, "cells.dat"), col_names = c("x", "y"))

##########################################################################################
#                                 Subproblem a)                                          #
##########################################################################################


