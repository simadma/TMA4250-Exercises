#Exercise 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-23-21

#libraries
source("./Exercise 2/R/utils.R")
##########################################################################################
#                                     PROBLEM 3                                          #
##########################################################################################

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
p1 <- ggplot(tib, aes(x, y)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = lab), alpha = 0.2) +
  geom_line(aes(color = "Empirical L-function")) +
  coord_fixed(
    xlim=c(0, .5), ylim=c(0, .5)
  ) +
  labs(x = "t", y = "L(t)", color = "", fill = "") +
  scale_color_manual(values = "blue") +
  scale_fill_manual(values = "black") +
  guides(fill  = guide_legend(order = 2),
         color = guide_legend(order = 1)) +
  theme_bw() +
  theme(
    legend.margin = margin(),
    legend.spacing.y = unit(0, "mm"),
    legend.background = element_blank()
  )

# w <- 7 + 4  # cm
# h <- 7  # cm
# ggsave("Neuman_Scott_L.pdf",
#   plot = p1, path = fig_path,
#   width = w, height = h, units = "cm"
# )


## Plot
# Simulated point pattern
set.seed(53)
p2 <- ggplot(sequential_sim_NS(lambda_M, lambda_C, sigma_c), aes(x, y, color = family)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.margin = margin(),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2))


p3 <- ggplot(sequential_sim_NS(lambda_M, lambda_C, sigma_c), aes(x, y, color = family)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.margin = margin(),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2))

p4 <- ggplot(sequential_sim_NS(lambda_M, lambda_C, sigma_c), aes(x, y, color = family)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.margin = margin(),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2))

p2
p3
p4

# w <- 7  # cm
# h <- 7 + 1  # cm
# ggsave("Neuman_Scott_sim1.pdf",
#   plot = p2, path = fig_path,
#   width = w, height = h, units = "cm"
# )
# 
# ggsave("Neuman_Scott_sim2.pdf",
#        plot = p3, path = fig_path,
#        width = w, height = h, units = "cm"
# )
# 
# ggsave("Neuman_Scott_sim3.pdf",
#        plot = p4, path = fig_path,
#        width = w, height = h, units = "cm"
# )



# Redwoods Point pattern
ggplot(redwood, aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()
