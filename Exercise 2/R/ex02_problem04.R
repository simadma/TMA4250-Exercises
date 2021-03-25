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


test <- iterative_sim_Strauss(nrow(cells), tau_0 = 0.07, phi_0 = 7000, phi_1 = 50)

test %>% 
  ggplot(aes(x, y)) +
  geom_point() +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_bw()

#Monte carlo test on the parameters with the L-function
sim_multiple_Strauss <- function(k, tau_0, phi_0, phi_1, nsim=100) {
  result <- matrix(0, nrow = nsim, ncol = 100) #hvorfor 100 kolonner?
  for (i in 1:nsim){
    result[i, ] <- L(iterative_sim_Strauss(k, tau_0, phi_0, phi_1))$y
  }
  result
}

#Monte carlo test on the parameters
#disse er ikke rett:/ er vanskelig å få den sånn flat på starten...
tau_0 <- 0.07
phi_0 <- 7000
phi_1 <- 90
set.seed(99)
L_sim <- sim_multiple_Strauss(nrow(cells), tau_0 = tau_0, phi_0 = phi_0, phi_1 = phi_1, nsim=10) #100 tar kjeeempelang tid
level <- 0.1
predint <- pred_interval(L_sim, level)
lab <- sprintf("%.2f-interval", 1 - level)
tib <- with(predint, cbind(L(cells), lower = lower, upper = upper))
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
p1
