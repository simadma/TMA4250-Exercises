#Project 3 in TMA4250
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#PROBLEM 1

# Load libraries
library(tidyverse)
library(spatial)
library(fields)

data_path <- "./Exercise 3/Data/"
fig_path <- "./Exercise 3/Figures/"

## Read data
# Seismic data observations
seismic <- scan(paste0(data_path, "seismic.dat"))
# Lithology distribution, with observation code {0, 1} for {sand, shale}
complit <- as.matrix(read.table(paste0(data_path, "complit.dat")))
dimnames(complit) <- NULL
n_seis <- as.integer(sqrt(length(seismic)))

# black and white color palette
bw <- gray.colors(256, start = 0, end = 1)
zlim <- c(-0.08, 1.08)
# mar_default <- c(5, 4, 4, 2) + 0.1
# Plot
# pdf(
#   file = paste0(fig_path, "obs.pdf"),
#   height = 4,
#   width = 5.2
# )
par(mar = c(2, 2, 0, 2) + 0.5)
image.plot(x = 1:n_seis, y = 1:n_seis, z = matrix(seismic, nrow=n_seis),
           col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# contour(x = 1:n_seis, y = 1:n_seis, z = matrix(seismic, nrow=n_seis),
#         nlevels = 3, drawlabels = FALSE, add = TRUE)
# dev.off()

image(complit, col=bw, asp=1)
##########################################################################################
##                                                                                      ##
##                               UNIFORM PRIOR MODEL                                    ##
##                                                                                      ##
##########################################################################################

##########################################################################################
##                 Six realizations of posterior with uniform prior                     ##
##########################################################################################
set.seed(81)
mu_0 <- 0.02
mu_1 <- 0.08
sigma_sq <- 0.06^2

# "success" probability p(l_i = 1 | d), i = 1,...,n
theta <- 1 / (1 + exp(((seismic - mu_1)^2 - (seismic - mu_0)^2) / (2*sigma_sq)))
nsim <- 6
posterior_unif_prior <- matrix(rbinom(n = nsim*length(theta), size = 1, prob = theta),
                               nrow = nsim, byrow = TRUE)
# Plot
# pdf(
#   file = paste0(fig_path, "post_real_unif_prior.pdf"),
#   height = 4,
#   width = 6.0
# )
par(mfrow = c(2, 3), mar = c(2, 2, 0, 0) + 0.5)
for (i in 1:nsim) {
  image(x = 1:n_seis, y = 1:n_seis, z = matrix(posterior_unif_prior[i, ], nrow=n_seis),
        col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
}
# dev.off()


##########################################################################################
##             Expectation and Variance of posterior with uniform prior                 ##
##########################################################################################

# Plot
# pdf(
#   file = paste0(fig_path, "post_exp_unif_prior.pdf"),
#   height = 4,
#   width = 5.2
# )
par(mar = c(2, 2, 0, 2) + 0.5)
image.plot(x = 1:n_seis, y = 1:n_seis, z = matrix(theta, nrow=n_seis),
           col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# contour(x = 1:n_seis, y = 1:n_seis, z = matrix(theta, nrow=n_seis),
#         nlevels = 3, drawlabels = FALSE, add = TRUE)
# dev.off()

# pdf(
#   file = paste0(fig_path, "post_var_unif_prior.pdf"),
#   height = 4,
#   width = 5.2
# )
par(mar = c(2, 2, 0, 2) + 0.5)
image.plot(x = 1:n_seis, y = 1:n_seis, z = matrix(theta*(1 - theta), nrow=n_seis),
           col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# contour(x = 1:n_seis, y = 1:n_seis, z = matrix(theta*(1 - theta), nrow=n_seis),
#         nlevels = 3, drawlabels = FALSE, add = TRUE)
# dev.off()

##########################################################################################
##                         MMAP of posterior with uniform prior                         ##
##########################################################################################

mu_critical <- (mu_0 + mu_1)/2
mmap <- 1*(seismic > mu_critical)

# Plot
# pdf(
#   file = paste0(fig_path, "mmap_unif_prior.pdf"),
#   height = 4,
#   width = 4.4
# )
par(mar = c(2, 2, 0, 2) + 0.5)
image(x = 1:n_seis, y = 1:n_seis, z = matrix(mmap, nrow=n_seis),
      col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


##########################################################################################
##                                                                                      ##
##                              MARKOV RF PRIOR MODEL                                   ##
##                                                                                      ##
##########################################################################################

