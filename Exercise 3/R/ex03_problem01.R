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

##########################################################################################
##                               Display observations                                   ##
##########################################################################################

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
#        nlevels = 3, drawlabels = FALSE, add = TRUE)
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
posterior_unif_prior <- matrix(1*(runif(nsim*length(theta)) < theta),
                               nrow = nsim, byrow = TRUE)

# Plot
pdf(
  file = paste0(fig_path, "post_real_unif_prior.pdf"),
  height = 4,
  width = 6.0
)
par(mfrow = c(2, 3), mar = c(2, 2, 0, 0) + 0.5)
for (i in 1:nsim) {
  image(x = 1:n_seis, y = 1:n_seis, z = matrix(posterior_unif_prior[i, ], nrow=n_seis),
        col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
}
dev.off()


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

##########################################################################################
##                          Observations from train data                                ##
##########################################################################################

n_comp <- nrow(complit)
# Plot
# pdf(
#   file = paste0(fig_path, "obs_train.pdf"),
#   height = 3.592,
#   width = 3.992
# )
par(mar = c(2, 2, 0, 2) + 0.5)
image(x = 1:n_comp, y = 1:n_comp, z = complit,
      col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# dev.off()

##########################################################################################
##            Compute Maximum Marginal pseudo-Likelihood (MMpL) estimate                ##
##########################################################################################
# MMpL <- function(beta, l) {
#   m <- nrow(l)
#   value <- 0
#   for (k_row in 1:m) {
#     for (k_col in 1:m) {
#       l_k <- l[k_row, k_col]
#       # first term
#       first_term <- 0
#       if (k_row > 1) first_term <- first_term + (l[k_row - 1, k_col] == l_k)  # v clique
#       if (k_row < m) first_term <- first_term + (l[k_row + 1, k_col] == l_k)  # v clique
#       if (k_col > 1) first_term <- first_term + (l[k_row, k_col - 1] == l_k)  # h clique
#       if (k_col < m) first_term <- first_term + (l[k_row, k_col + 1] == l_k)  # h clique
#       value <- value + log(beta)*first_term
#       
#       # second term
#       second_term <- 0
#       for (l_prime in 0:1) {
#         temp <- 0
#         if (k_row > 1) temp <- temp + (l[k_row - 1, k_col] == l_prime)  # v clique
#         if (k_row < m) temp <- temp + (l[k_row + 1, k_col] == l_prime)  # v clique
#         if (k_col > 1) temp <- temp + (l[k_row, k_col - 1] == l_prime)  # h clique
#         if (k_col < m) temp <- temp + (l[k_row, k_col + 1] == l_prime)  # h clique
#         second_term <- second_term + beta^temp
#       }
#       value <- value - log(second_term)
#     }
#   }
#   value
# }
MMpL <- function(beta, l) {
  m <- nrow(l)
  value <- 0
  for (k_row in 1:m) {
    for (k_col in 1:m) {
      l_k <- l[k_row, k_col]
      # first term
      first_term <- ifelse(
        test = k_row > 1,
        yes  = (l[k_row - 1, k_col] == l_k),
        no   = (l[m, k_col] == l_k)
      ) +
        ifelse(
          test = k_row < m,
          yes  = (l[k_row + 1, k_col] == l_k),
          no   = (l[1, k_col] == l_k)
        ) +
        ifelse(
          test = k_col > 1,
          yes  = (l[k_row, k_col - 1] == l_k),
          no   = (l[k_row, m] == l_k)
        ) +
        ifelse(
          test = k_col < m,
          yes  = (l[k_row, k_col + 1] == l_k),
          no   = (l[k_row, 1] == l_k)
        )
      value <- value + log(beta)*first_term
      
      # second term
      second_term <- 0
      for (l_prime in 0:1) {
        temp <- ifelse(
          test = k_row > 1,
          yes  = (l[k_row - 1, k_col] == l_prime),
          no   = (l[m, k_col] == l_prime)
        ) +
          ifelse(
            test = k_row < m,
            yes  = (l[k_row + 1, k_col] == l_prime),
            no   = (l[1, k_col] == l_prime)
          ) +
          ifelse(
            test = k_col > 1,
            yes  = (l[k_row, k_col - 1] == l_prime),
            no   = (l[k_row, m] == l_prime)
          ) +
          ifelse(
            test = k_col < m,
            yes  = (l[k_row, k_col + 1] == l_prime),
            no   = (l[k_row, 1] == l_prime)
          )
        second_term <- second_term + beta^temp
      }
      value <- value - log(second_term)
    }
  }
  value
}

beta_hat <- optimize(f = MMpL, interval = c(1, 100), l = complit, maximum = TRUE)$maximum
beta_hat

##########################################################################################
##                                Realizations                                          ##
##########################################################################################

#function to check proportion of sand
prop_sand <- function(l) {
  sum(l==0)/length(l)
}

#function to find the sum of I(l_i=1) for i in the neighborhood of k
neighborhood <- function(k, l, n, ncols){
  #if l is a vector
  
  #check if k is such that west point is outside (left column in matrix)
  if (k %in% seq(1, n, ncols)) l_W<-l[k + (ncols - 1)]
  else l_W <- l[k-1]
  #check if k is such that east point is outside (right column in matrix)
  if (k %in% seq(ncols,n, ncols)) l_E<-l[k - (ncols - 1)]
  else l_E <- l[k+1]
  #check if k is such that north point is outside (upper row in matrix)
  if (k %in% 1:ncols) l_N <- l[k + (ncols-1)*ncols]
  else l_N <- l[k-ncols]
  #check if k is such that south point is outside (lower row in matrix)
  if (k %in% (n-(ncols-1)):n) l_S <- l[k - (ncols-1)*ncols]
  else l_S <- l[k+ncols]
  
  return(l_E+l_W+l_N+l_S)
}


#function to compute lambda_k (successprobability) for the bernoulli dist of l_i|d
successprob <- function(neighbors, beta, d_k){
  return ((dnorm(d_k, mean=mu_1, sd=sqrt(sigma_sq)) *beta^(neighbors)) /
            (dnorm(d_k, mean=mu_0, sd=sqrt(sigma_sq))* beta^(4-neighbors) + 
               dnorm(d_k, mean=mu_1, sd=sqrt(sigma_sq))* beta^(neighbors)))
}


#proposals from the bernoulli distribution with probability of success = lambda_k
g <- function(l, d, beta){
  i <- sample(1:length(l), 1)
  neighbors <- neighborhood(k=i, l, length(l), as.integer(sqrt(length(l))))
  # cat(neighbors)
  lambda_i <- successprob(neighbors = neighbors, beta=beta, d[i])
  #cat(lambda_i)
  l_i <- rbinom(1,1,lambda_i)
  l[i] <- l_i
  return (l)
}

#actual algorithm to perform gibbs sampling
gibbsposterior <- function(l_0, d, beta, max_iter){
  n <- length(l_0)
  nsamples<- as.integer(max_iter/n)
  samples <- matrix(NA, nrow=nsamples, ncol=n)
  convergence <- matrix(NA, nrow=nsamples)
  l <- l_0
  sample_iter <- 1
  for (iter in 1:max_iter){
    #l^j drawn from g(l|l^{j-1}) (block gibbs)
    l <- g(l, d, beta)
    #save every n samples
    if ( iter %% n ==0) {
      samples[sample_iter,] <- l
      #check convergence with proportion of sand
      sandprob <- prop_sand(l)
      convergence[sample_iter] <- sandprob
      sample_iter <- sample_iter+1
    }
  }
  return (list(samples, convergence)) #return all samples of l and the sand proportion
}

##########################################################################################
##                                Initial guess: All sand (all 0s)                      ##
##########################################################################################

nsweeps <- 500

l_00 <- rep(0, length.out = length(seismic))
result0 <- gibbsposterior(l_00, seismic, beta_hat, length(seismic)*nsweeps)

proportion0 <- result0[[2]][,1]
distribution0 <- result0[[1]]


##########################################################################################
##                                Initial guess: All shale (all 1s)                     ##
##########################################################################################

l_01 <- rep(1, length.out = length(seismic))
result1 <- gibbsposterior(l_01, seismic, beta_hat, length(seismic)*nsweeps)

proportion1 <- result1[[2]][,1]
distribution1 <- result1[[1]]



##########################################################################################
##                                Initial guess: Uniform iid                            ##
##########################################################################################

l_0r <- sample(0:1, size = length(seismic), replace = TRUE) 
resultr <- gibbsposterior(l_0r, seismic, beta_hat, length(seismic)*nsweeps)

proportionr <- resultr[[2]][,1]
distributionr <- resultr[[1]]



#####################################
##              PLOTTING           ##
#####################################
# Plot
# pdf(
#   file = paste0(fig_path, "convergenceposterior.pdf"),
#   height = 4.592,
#   width = 4.992
# )
plot(1:nsweeps, proportion0, type="l", col = "blue", ylim = c(0,1), ylab = "Sand proportion", xlab = "Sweep number")
lines(proportion1, col = "red")
lines(proportionr, col = "black")
legend("bottomright", legend = c("All sand", "All shale", "Uniform iid"), col = c("blue", "red", "black"), lty=c(1,1,1))
# dev.off()


##########################################################################################
##                 Six realizations of posterior with uniform prior                     ##
##########################################################################################

realizationid <- sample(250:nsweeps, 2)



sample1 <- distribution0[realizationid[1],]
sample2 <- distribution0[realizationid[2],]

sample3 <- distribution1[realizationid[1],]
sample4 <- distribution1[realizationid[2],]

sample5 <- distributionr[realizationid[1],]
sample6 <- distributionr[realizationid[2],]

# Plot
# pdf(
#   file = paste0(fig_path, "post_real_markov_prior.pdf"),
#   height = 4,
#   width = 6.0
# )
# par(mfrow = c(2, 3), mar = c(2, 2, 0, 0) + 0.5)
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(sample1, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(sample2, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(sample3, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(sample4, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(sample5, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(sample6, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# dev.off()


##########################################################################################
##             Expectation and Variance of posterior with Markov prior                  ##
##########################################################################################

samples_combined <- rbind(distribution0[250:nsweeps,], rbind(distribution1[150:nsweeps,], distributionr[50:nsweeps,]))

mean_post <- apply(samples_combined, MARGIN = 2, FUN = mean)
var_post <- apply(samples_combined, MARGIN = 2, FUN = var)

# Plot
# pdf(
#   file = paste0(fig_path, "post_exp_markov_prior.pdf"),
#   height = 4,
#   width = 5.2
# )
# par(mar = c(2, 2, 0, 2) + 0.5)
# image.plot(x = 1:n_seis, y = 1:n_seis, z = matrix(mean_post, nrow=n_seis),
#            col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# # contour(x = 1:n_seis, y = 1:n_seis, z = matrix(mean_post, nrow=n_seis),
# #         nlevels = 3, drawlabels = FALSE, add = TRUE)
# dev.off()
# 
# pdf(
#   file = paste0(fig_path, "post_var_markov_prior.pdf"),
#   height = 4,
#   width = 5.2
# )
# par(mar = c(2, 2, 0, 2) + 0.5)
# image.plot(x = 1:n_seis, y = 1:n_seis, z = matrix(var_post, nrow=n_seis),
#            col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# # contour(x = 1:n_seis, y = 1:n_seis, z = matrix(var_post, nrow=n_seis),
# #         nlevels = 3, drawlabels = FALSE, add = TRUE)
# dev.off()

##########################################################################################
##                         MMAP of posterior with Markov prior                          ##
##########################################################################################

mu_critical2 <- 0.5
mmap2 <- 1*(mean_post > mu_critical2)

# Plot
# pdf(
#   file = paste0(fig_path, "mmap_markov_prior.pdf"),
#   height = 4,
#   width = 4.4
# )
# par(mar = c(2, 2, 0, 2) + 0.5)
# image(x = 1:n_seis, y = 1:n_seis, z = matrix(mmap2, nrow=n_seis),
#       col=bw, zlim = zlim, asp = 1, xlab = "", ylab = "")
# dev.off()

