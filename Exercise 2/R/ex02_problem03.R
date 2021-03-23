source("./Exercise 2/R/utils.R")

## Read data
# California Redwoods Point Pattern
redwood <- read_table(paste0(data_path, "redwood.dat"), col_names = c("x", "y"))

## Simulate Mother-Children model
sequential_sim_NS <- function(lambda_M, lambda_C, sigma_c) {
  k_D_Mothers <- rpois(1, lambda = lambda_M)
  X_D <- tibble(x = numeric(), y = numeric(), family = factor())
  for (j in 1:k_D_Mothers) {
    x_Mother <- runif(2)
    k_j_Children <- rpois(1, lambda = lambda_C)
    if (k_j_Children > 0) {
      x_ji_ps <- mvrnorm(k_j_Children, mu = x_Mother, Sigma = diag(sigma_c^2, 2))
      x_ji_s <- x_ji_ps %% 1  # Torus wrap
      if (k_j_Children == 1) {
        x_ji_s <- t(as.matrix(x_ji_s))  # Reshape to row vector
      }
      X_D <- add_row(X_D,
               x      = x_ji_s[, 1],
               y      = x_ji_s[, 2],
               family = as.factor(j)
             )
    }
  }
  X_D
}

sim_multiple_NS <- function(lambda_M, lambda_C, sigma_c, nsim = 100) {
  result <- matrix(0, nrow = nsim, ncol = 100)
  for (i in 1:nsim) {
    result[i, ] <- L(sequential_sim_NS(lambda_M, lambda_C, sigma_c))$y
  }
  result
}

# 6.55, .09
lambda_M = 6.55
lambda_C = nrow(redwood)/lambda_M
sigma_c = .09
set.seed(41)
L_sim <- sim_multiple_NS(lambda_M, lambda_C, sigma_c)
level <- 0.1
predint <- pred_interval(L_sim, level)
lab <- sprintf("%.2f-interval", 1 - level)
tib <- with(predint, cbind(L(redwood), lower = lower, upper = upper))
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


## Plot
# Simulated point pattern
lambda_M = 6.5
lambda_C = nrow(redwood)/lambda_M
sigma_c = .08
set.seed(3)
ggplot(sequential_sim_NS(lambda_M, lambda_C, sigma_c), aes(x, y, color = family)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

# Redwoods Point pattern
ggplot(redwood, aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

plot_Lf(redwood)
plot_Lf(X_D)
