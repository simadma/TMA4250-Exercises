#Problem 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-08-21

#fig_path
fig_path <- "C:/Users/karin/OneDrive - NTNU/8. semester/Romlig statistikk/TMA4250-Exercises/Exercise 2/Figures/"
## Figure size
w <- 10  # cm
h <- 6  # cm
#packages
library(ggplot2)
library(dplyr)

# plot2d <- function(treecount,filename="NONE", save=FALSE){
#   df2d <- data.frame(x=d$x, y=d$y, N_obs = treecount)
#   fig2di <- ggplot(data=df2d, aes(x=x, y=y, color=as.factor(N_obs))) +
#     geom_point(size=0.75) +
#     scale_colour_manual(values = cbPalette) +
#     labs(color="Tree count")+
#     theme_minimal()
#   if (save==TRUE){
#     ggsave(filename=filename,
#            plot = fig2di, path = fig_path,
#            width = w, height = h, units = "cm")
#   }  
# }

eventloc_gridnode <- function(x,y,N_obs) {
  x_loc <- runif(N_obs, min=x-5, max=x+5)
  y_loc <- runif(N_obs, min=y-5, max=y+5)
  return(data.frame(x=x_loc, y=y_loc))
}




eventloc <- function(df_event_count) {
  eventlocs<-data.frame(x=double(), y=double())
  for (i in 1:nrow(df_event_count)){
    if (df_event_count$N_obs[i]>0){
      temp <- with(df_event_count, eventloc_gridnode(x[i], y[i], N_obs[i]))
      eventlocs <- rbind(eventlocs, temp)
    }
  }
  return(eventlocs)
}

plot_eventloc1 <- function(eventlocs) {
  ggplot(data=eventlocs, aes(x=x, y=y, color="Observed")) +
    geom_vline(xintercept = seq(0,300, 10), color="black", alpha = 0.1, size=0.1) +
    geom_hline(yintercept = seq(0,300, 10), color="black", alpha = 0.1, size=0.1) +
    geom_point(size=0.5) +
    coord_fixed() +
    theme_minimal() +
    labs(color="") + 
    theme(legend.position = "none")
  
}

plot_eventloc <- function(eventlocs) {
  ggplot(data=eventlocs, aes(x=x, y=y)) +
    geom_vline(xintercept = seq(0,300, 10), color="black", alpha = 0.1, size=0.1) +
    geom_hline(yintercept = seq(0,300, 10), color="black", alpha = 0.1, size=0.1) +
    geom_point(size=0.5, aes(color=point)) +
    coord_fixed() +
    theme_minimal() +
    labs(color="") + 
    theme(legend.position = "top")

}



#PATHS FOR DATASET
obspath <- "C:/Users/karin/OneDrive - NTNU/8. semester/Romlig statistikk/TMA4250-Exercises/Exercise 2/Data/obspines.txt"
obsprobpath <- "C:/Users/karin/OneDrive - NTNU/8. semester/Romlig statistikk/TMA4250-Exercises/Exercise 2/Data/obsprob.txt"
d <- read.delim(obspath, header = TRUE, sep=" ")
alpha <- read.delim(obsprobpath, header = TRUE, sep=" ")


#Event location plot for observations
d_eventloc <- eventloc(d)
fig1 <- plot_eventloc1(d_eventloc)
fig1

#Plot of the observation probability
fig2 <- ggplot(data=alpha) +
  geom_tile(aes(x=x, y=y, fill=alpha))+
  theme_minimal()
#  theme(legend.position = "top")

fig2

#Saving figures for problem 1a
ggsave("obspinetrees.pdf",
  plot = fig1, path = fig_path,
  width = w, height = h, units = "cm"
)


ggsave("probpinetrees.pdf",
       plot = fig2, path = fig_path,
       width = w, height = h, units = "cm"
)

#c) Estimation of lambda
area <- 100 #m^2

lambda_hat <- sum(d$N_obs)/(area*sum(alpha$alpha))
cat("Estimation of lambda_k is ", lambda_hat, "per square meter")

# Generate 6 realizations from the prior
nsim <- 6 #number of simulations
n <- 900 #number of grid points
priors <- matrix(0, nrow=nsim, ncol=n) #for storage
for (i in 1:nsim){
  priors[i,] <- rpois(n=n, lambda=lambda_hat*area)
}


#plot all realizations
prior_dfs <- list()
for (i in 1:nsim){
  df_i <- data.frame(x=d$x, y=d$y, N_obs = priors[i,])
  eventlocs_i <- eventloc(df_i)
  prior_dfs[[i]] <- eventlocs_i
}

a <- bind_rows(prior_dfs, .id="realization")
fig2c <- plot_eventloc1(a) + facet_wrap(~realization, nrow=2)
fig2c
w2 <- 6*2  # cm
h2 <- 6*2  # cm
ggsave(filename="02priors.pdf", plot=fig2c, path=fig_path, height = h2, width=w2, units = "cm")


## 2 d) 
#make 6 realizations from posterior
#loop through the grid to draw samples
posterior <- matrix(NA, nrow=n, ncol=nsim)
for (i in 1:n) {
  unobs <- rpois(n=nsim, lambda = (1-alpha$alpha[i])*lambda_hat*area)
  grid_i_count <- unobs#+d$N_obs[i]
  posterior[i,]<- grid_i_count
}

#event locations for all points
posterior_dfs <- list()
for (i in 1:nsim){
  df_i <- data.frame(x=d$x, y=d$y, N_obs = posterior[,i])
  eventlocs_i <- eventloc(df_i)
  tot_eventloc <- rbind(cbind(eventlocs_i,point="Predicted"), cbind(d_eventloc,point="Observed"))
  posterior_dfs[[i]]<-tot_eventloc
}


b <- bind_rows(posterior_dfs, .id="realization")
fig2d <- plot_eventloc(b) + facet_wrap(~realization, nrow=2)

fig2d
## Figure size
w2 <- 6*2  # cm
h2 <- 6*2  # cm
ggsave(filename="02posteriors.pdf", plot=fig2d, path=fig_path, height = h2, width=w2, units = "cm")






#2 e
#Simulate 100 realizations of the discretized event-count model, both for
#the prior and the posterior models

#Prior
prior100<-matrix(NA, nrow=100, ncol=n)
for (i in 1:100){
  prior100[i,] <- rpois(n=n, lambda=lambda_hat*area)
}
mean_prior <- apply(prior100, MARGIN = 2, mean)





#posterior
posterior100<- matrix(NA, nrow=100, ncol=n)
for (i in 1:n) {
  unobs <- rpois(n=100, lambda = (1-alpha$alpha[i])*lambda_hat*area)
  grid_i_count <- unobs+d$N_obs[i]
  posterior100[,i]<- grid_i_count
}
#Find the mean
mean_posterior <- apply(posterior100, MARGIN = 2, mean)



#For plotting
rng = range(c((mean_prior), (mean_posterior))) #a range to have the same min and max for both plots

#dataframes for plotting
df.prior <- data.frame(x=d$x, y=d$y, mean=mean_prior)
df.posterior <- data.frame(x=d$x, y=d$y, mean=mean_posterior)


mprior <- ggplot(data=df.prior) +
  geom_tile(aes(x=x, y=y, fill=mean))+
  scale_fill_gradient2(low="blue", mid="cyan", high="lightblue", #colors in the scale
                       midpoint=mean(rng),    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,4), #breaks in the scale bar
                       limits=c(floor(rng[1]), ceiling(rng[2]))) + #same limits for plots
  labs(title="Mean tree count for the prior")+
  theme_minimal()

mprior


mposterior <- ggplot(data=df.posterior) +
  geom_tile(aes(x=x, y=y, fill=mean))+
  scale_fill_gradient2(low="blue", mid="cyan", high="lightblue", #colors in the scale
                       midpoint=mean(rng),    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,4), #breaks in the scale bar
                       limits=c(floor(rng[1]), ceiling(rng[2])))+ #same limits for plots
  labs(title="Mean tree count for the posterior")+
  theme_minimal()

mposterior


ggsave("02eprior100.pdf",
       plot = mprior, path = fig_path,
       width = w, height = h, units = "cm"
)
ggsave("02eposterior100.pdf",
       plot = mposterior, path = fig_path,
       width = w, height = h, units = "cm"
)
