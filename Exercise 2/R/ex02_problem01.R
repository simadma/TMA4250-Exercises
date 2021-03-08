library(tidyverse)
data_path <- "./Exercise 2/Data/"
fig_path <- "./Exercise 2/Figures/"

##########################################################################################
#                                 Subproblem a)                                          #
##########################################################################################

## Read data

# Biological Cells Point Pattern
cells <- read_table(paste0(data_path, "cells.dat"), col_names = c("x", "y"))
# California Redwoods Point Pattern
redwood <- read_table(paste0(data_path, "redwood.dat"), col_names = c("x", "y"))
# Japanese Pines Point Pattern
pines <- read_table(paste0(data_path, "pines.dat"), col_names = c("x", "y"), skip = 3)


## Plot data

p1 <- ggplot(cells, aes(x, y)) +
  geom_point() +
  coord_fixed() +
  labs(title = "Biological Cells Point Pattern") +
  theme_bw()

p2 <- ggplot(redwood, aes(x, y)) +
  geom_point() +
  coord_fixed() +
  labs(title = "California Redwoods Point Pattern") +
  theme_bw()

p3 <- ggplot(pines, aes(x, y)) +
  geom_point() +
  coord_fixed() +
  labs(title = "Japanese Pines Point Pattern") +
  theme_bw()


## Save figures

# ggsave("cells.pdf",
#   plot = p1, path = fig_path,
#   width = 9, height = 9, units = "cm"
# )
# ggsave("redwood.pdf",
#        plot = p2, path = fig_path,
#        width = 9, height = 9, units = "cm"
# )
# ggsave("pines.pdf",
#        plot = p3, path = fig_path,
#        width = 9, height = 9, units = "cm"
# )
