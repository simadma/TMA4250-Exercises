rm(list = ls())
library(geoR)    # Analysis of Geostatistical Data
library(akima)   # Interpolation of Irregularly and Regularly Spaced Data
library(fields)  # Tools for Spatial Data
library(MASS)    # Support Functions and Datasets for Venables and Ripley's MASS
library(tidyverse)

## Problem 1: Gaussian RF - model characteristics

# Dataframe of all combinations of parameters
params <- rbind(
  expand.grid(
    cov.model = "powered.exponential",
    nu        = c(1, 1.9),
    sigma2    = c(1, 5),
    stringsAsFactors = FALSE
  ),
  expand.grid(
    cov.model = "matern",
    nu        = c(1, 3),
    sigma2    = c(1, 5),
    stringsAsFactors = FALSE
  )
)

# Initialize dataframe to store results from `geoR::cov.spatial()`
df <- data.frame(
  tau       = numeric(),
  cov       = numeric(),
  cov.model = character(),
  nu        = numeric(),
  sigma2    = numeric()
)

n <- 50  # Number of grid points
distances <- 1:n - 1
phi <- 10
for (i in 1:nrow(params)) {
  cov.model <- params[i, ]$cov.model
  nu        <- params[i, ]$nu
  sigma2    <- params[i, ]$sigma2
  covs <- cov.spatial(
    obj       = distances,
    cov.model = cov.model,
    cov.pars  = c(sigma2, phi),
    kappa     = nu
  )
  df <- rbind(  # Add result to dataframe
    df,
    data.frame(
      tau       = distances,
      cov       = covs,
      cov.model = cov.model,
      nu        = nu,
      sigma2    = sigma2
    )
  )
}

# Plot
df %>%
  filter(sigma2 == 1) %>%
  ggplot(mapping = aes(x = tau, y = cov, color = factor(nu), linetype = cov.model)) +
  geom_line(size = 1) +
  theme_minimal()

# Computes the covariance matrix
# cov: vector of covariances of pairs
cov_vec2cov_mat <- function(cov) {
  n <- length(cov)
  Sigma <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      h <- abs(i - j)
      Sigma[i, j] <- cov[h + 1]
    }
  }
  Sigma
}


autocov <- df %>%
  filter(
    cov.model == "powered.exponential",
    nu        == 1,
    sigma2    == 1
  ) %>% 
  pull(cov)



Sigma <- cov_vec2cov_mat(autocov)

test <- mvrnorm(n = 4, mu = rep(0, n), Sigma = Sigma)

test <- as.data.frame(t(test), row.names = c("h", "hd", "hf2", "hfs"))


paste0("realization", 1:4)
