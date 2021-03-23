#Exercise 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-23-21

#libraries
source("./Exercise 2/R/utils.R")
##########################################################################################
#                                     PROBLEM 4                                          #
##########################################################################################




## Read data
# Biological Cells Point Pattern
cells <- read_table(paste0(data_path, "cells.dat"), col_names = c("x", "y"))

##########################################################################################
#                                 Subproblem a)                                          #
##########################################################################################

# Biological Cells point pattern
cells %>% ggplot(aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

plot_Lf(cells)
set.seed(10)
g <- Strauss(42,c=0.04, r = 0.10)
plot_Lf(g)
data.frame(x=g$x, y=g$y) %>% 
  ggplot(aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()


phi <- function(tau, tau_0, phi_0, phi_1) {
  tau <- norm(tau, type = "2")
  if_else(tau <= tau_0,
    true  = phi_0,
    false = phi_0*exp(-phi_1*(tau - tau_0))
  )
}

iterative_sim_Strauss <- function(k, tau_0, phi_0, phi_1, n_iter = 1000) {
  X_D <- matrix(runif(2*k), nrow = k, dimnames = list(NULL, c("x", "y")))
  for (i in 1:n_iter) {
    u <- sample(k, 1) # sample index
    x_p <- runif(2)   # proposal
    exponent <- 0
    for (j in 1:k) {
      exponent <- exponent + (phi(x_p - X_D[j, ], tau_0, phi_0, phi_1)
                              - phi(X_D[u, ] - X_D[j, ], tau_0, phi_0, phi_1))
    }
    alpha <- min(1, exp(-exponent))
    if (runif(1) < alpha) {
      X_D[u, ] <- x_p
    }
  }
  as_tibble(X_D)
}


test <- iterative_sim_Strauss(nrow(cells), tau_0 = 0.08, phi_0 = 5, phi_1 = 6)

test %>% 
  ggplot(aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()
