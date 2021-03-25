library(MASS)
library(tidyverse)
library(spatial)
data_path <- "./Exercise 2/Data/"
fig_path <- "./Exercise 2/Figures/"

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
plot_Lf <- function(data, fs = .5, k = 100, dmin = FALSE) {
  lab <- c("L(t) = t", "Empirical L-function")
  ggplot(L(data, fs, k, dmin), aes(x, y)) +
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


pred_interval <- function(simulations, level = 0.10) {
  mu <- apply(simulations, 2, mean)
  sigma <- apply(simulations, 2, sd)
  crit_val <- qt(level/2, nrow(simulations), lower.tail = FALSE)
  list(
    lower = mu - crit_val*sigma,
    upper = mu + crit_val*sigma
  )
}
