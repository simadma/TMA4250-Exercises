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
L <- function(data, fs = .5, k = 100, dmin = FALSE) {
  K <- with(data, Kfn(pp = list(x = x, y = y), fs, k))
  if (dmin) {
    print(paste("Minimum distance:", K$dmin))
  }
  with(K, tibble(x = x, y = y))
}

# Plots the empirical and theoretical L(t)
plot_Lf <- function(data) {
  lab <- c("L(t) = t", "Empirical L-function")
  ggplot(L(data), aes(x, y)) +
    geom_function(fun=function(x) x, aes(color = lab[1], linetype = lab[1])) +
    geom_line(aes(color = lab[2], linetype = lab[2])) +
    scale_color_manual(values = c("blue", "black")) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_fixed() +
    labs(x = "t", y = "L(t)", color = "", linetype = "") +
    theme_bw() +
    theme(
      legend.position   = c(.35, .87),
      legend.background = element_blank()
    )
}

## Plot
p4 <- plot_Lf(cells)
p5 <- plot_Lf(redwood)  
p6 <- plot_Lf(pines)

p4
p5
p6

## Save figures

# w <- 7  # cm
# h <- 7  # cm
# ggsave("L_cells.pdf",
#   plot = p4, path = fig_path,
#   width = w, height = h, units = "cm"
# )
# ggsave("L_redwood.pdf",
#        plot = p5, path = fig_path,
#        width = w, height = h, units = "cm"
# )
# ggsave("L_pines.pdf",
#        plot = p6, path = fig_path,
#        width = w, height = h, units = "cm"
# )


##########################################################################################
#                                 Subproblem c)                                          #
##########################################################################################
sim_binom_process <- function(ncount, nsim = 100) {
  result <- matrix(0, nrow = nsim, ncol = 100)
  for (i in 1:nsim) {
    result[i, ] <- L(Psim(ncount))$y
  }
  result
}

pred_interval <- function(simulations, level = 0.10) {
  mu <- apply(simulations, 2, mean)
  sigma <- apply(simulations, 2, sd)
  crit_val <- qt(level/2, nrow(simulations), lower.tail = FALSE)
  list(
    lower = mu - crit_val*sigma,
    upper = mu + crit_val*sigma
  )
}

# Plots the empirical and theoretical L(t)
plot_Lf_with_interval <- function(data, nsim = 100, level = 0.10) {
  simulations <- sim_binom_process(nrow(data), nsim)
  predint <- pred_interval(simulations, level)
  lab <- sprintf("%.2f-interval", 1 - level)
  tib <- with(predint, cbind(L(data), lower = lower, upper = upper))
  ggplot(tib, aes(x, y)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = lab), alpha = 0.2) +
    geom_line(aes(color = "Empirical L-function")) +
    coord_fixed() +
    labs(x = "t", y = "L(t)", color = "", fill = "") +
    scale_color_manual(values = "blue") +
    scale_fill_manual(values = "black") +
    guides(fill  = guide_legend(order = 2),
           color = guide_legend(order = 1)) +
    theme_bw() +
    theme(
      legend.position   = c(.38, .87),
      legend.margin = margin(),
      legend.spacing.y = unit(0, "mm"),
      legend.background = element_blank()
    )
}

set.seed(71)
p7 <- plot_Lf_with_interval(cells)
p8 <- plot_Lf_with_interval(redwood)
p9 <- plot_Lf_with_interval(pines)

p7
p8
p9

## Save figures

# w <- 7  # cm
# h <- 7  # cm
# ggsave("sim_L_cells.pdf",
#   plot = p7, path = fig_path,
#   width = w, height = h, units = "cm"
# )
# ggsave("sim_L_redwood.pdf",
#        plot = p8, path = fig_path,
#        width = w, height = h, units = "cm"
# )
# ggsave("sim_L_pines.pdf",
#        plot = p9, path = fig_path,
#        width = w, height = h, units = "cm"
# )
