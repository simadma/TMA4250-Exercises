library(tidyverse)
library(spatial)
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
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

p2 <- ggplot(redwood, aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

p3 <- ggplot(pines, aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

p1
p2
p3
## Save figures

# w <- 7  # cm
# h <- 7  # cm
# ggsave("cells.pdf",
#   plot = p1, path = fig_path,
#   width = w, height = h, units = "cm"
# )
# ggsave("redwood.pdf",
#        plot = p2, path = fig_path,
#        width = w, height = h, units = "cm"
# )
# ggsave("pines.pdf",
#        plot = p3, path = fig_path,
#        width = w, height = h, units = "cm"
# )


##########################################################################################
#                                 Subproblem b)                                          #
##########################################################################################

ppregion(xl = 0, xu = 1, yl = 0, yu = 1)  # Set up domain D (This is default)

# Computes the expected number of points within distance t of a point of the pattern
L <- function(data) {
  K <- with(data,
            Kfn(
              pp = list(x = x, y = y),
              fs = 1,
              k  = 100
            )
       )
  with(K, tibble(x = x, y = y))
}

# Plots the empirical and theoretical L(t)
plot_Lf <- function(data) {
  ggplot(L(data), aes(x, y)) +
    geom_step(aes(color = "Empirical L-function")) +
    geom_function(fun=function(x) x, aes(color="L(t) = t")) +
    coord_fixed() +
    labs(x = "t", y = "L(t)") +
    theme_bw() +
    theme(legend.position = c(.20, .91), legend.title = element_blank())
}

## Plot
p4 <- plot_Lf(cells)
p5 <- plot_Lf(redwood)  
p6 <- plot_Lf(pines)

p4
p5
p6
## Save figures

w <- 7  # cm
h <- 7  # cm
ggsave("L_cells.pdf",
  plot = p4, path = fig_path,
  width = w, height = h, units = "cm"
)
ggsave("L_redwood.pdf",
       plot = p5, path = fig_path,
       width = w, height = h, units = "cm"
)
ggsave("l_pines.pdf",
       plot = p6, path = fig_path,
       width = w, height = h, units = "cm"
)