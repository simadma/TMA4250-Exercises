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

plot2d <- function(treecount,filename="NONE", save=FALSE){
  df2d <- data.frame(x=d$x, y=d$y, N_obs = treecount)
  fig2di <- ggplot(data=df2d, aes(x=x, y=y, color=as.factor(N_obs))) +
    geom_point(size=0.75) +
    scale_colour_manual(values = cbPalette) +
    labs(color="Tree count")+
    theme_minimal()
  if (save==TRUE){
    ggsave(filename=filename,
           plot = fig2di, path = fig_path,
           width = w, height = h, units = "cm")
  }  
}

obspath <- "C:/Users/karin/OneDrive - NTNU/8. semester/Romlig statistikk/TMA4250-Exercises/Exercise 2/Data/obspines.txt"
obsprobpath <- "C:/Users/karin/OneDrive - NTNU/8. semester/Romlig statistikk/TMA4250-Exercises/Exercise 2/Data/obsprob.txt"
d <- read.delim(obspath, header = TRUE, sep=" ")
alpha <- read.delim(obsprobpath, header = TRUE, sep=" ")

# fig1 <- ggplot(data=d) +
#   geom_tile(aes(x=x, y=y, fill=N_obs))
# for (i in 1:29){
#   fig1 + geom_vline(aes(xintercept=i*10))
#   fig1 + geom_hline(aes(yintercept=i*10))  
# }  
# fig1
cbPalette <- c("white", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
fig1 <- ggplot(data=d, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point(size=0.9) +
  scale_colour_manual(values = cbPalette) +
  labs(color="Tree count")+
  labs(title="Observed number of pine trees") +
  theme_minimal()
#  theme(legend.position = "top")

fig1

fig2 <- ggplot(data=alpha) +
  geom_tile(aes(x=x, y=y, fill=alpha))+
  labs(title="The probability of observing pine trees")+
  theme_minimal()
#  theme(legend.position = "top") 

fig2


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

for (i in 1:nsim){
  plot2d(priors[i,], filename=paste0("02realization", i, ".pdf"), save=TRUE)
}



## 2 d) 
#make 6 realizations from posterior
#loop through the grid to draw samples
result2d <- matrix(NA, nrow=n, ncol=nsim)
for (i in 1:n) {
  unobs <- rpois(n=nsim, lambda = (1-alpha$alpha[i])*lambda_hat*area)
  grid_i_count <- unobs+d$N_obs[i]
  result2d[i,]<- grid_i_count
}


for (i in 1:nsim){
  plot2d(result2d[,i], filename=paste0("posterior",i,".pdf"), save=TRUE)
}

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

mean_posterior <- apply(posterior100, MARGIN = 2, mean)



#For plotting
rng = range(c((mean_prior), (mean_posterior))) #a range to have the same min and max for both plots

#dataframes for plotting
df.prior <- data.frame(x=d$x, y=d$y, mean=mean_prior)
df.posterior <- data.frame(x=d$x, y=d$y, mean=mean_posterior)


mprior <- ggplot(data=df.prior) +
  geom_tile(aes(x=x, y=y, fill=mean))+
  scale_fill_gradient2(low="blue", mid="cyan", high="purple", #colors in the scale
                       midpoint=mean(rng),    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,4), #breaks in the scale bar
                       limits=c(floor(rng[1]), ceiling(rng[2]))) + #same limits for plots

  labs(title="Mean tree count for the prior")+
  theme_minimal()

mprior


mposterior <- ggplot(data=df.posterior) +
  geom_tile(aes(x=x, y=y, fill=mean))+
  scale_fill_gradient2(low="blue", mid="cyan", high="purple", #colors in the scale
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
