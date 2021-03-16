#Problem 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-08-21

#fig_path
fig_path <- "C:/Users/karin/OneDrive - NTNU/8. semester/Romlig statistikk/TMA4250-Exercises/Exercise 2/Figures/"
#packages
library(ggplot2)

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

## Save figures
w <- 10  # cm
h <- 6  # cm
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
realizations <- matrix(0, nrow=nsim, ncol=n) #for storage
for (i in 1:nsim){
  realizations[i,] <- rpois(n=n, lambda=lambda_hat)
}

for (i in 1:nsim){
  sim_i <- data.frame(x=d$x, y=d$y, N_obs=realizations[i,])
  figc_i<-ggplot(data=sim_i, aes(x=x, y=y, color=as.factor(N_obs))) +
    geom_point() +
    labs(color="Tree count")+
    scale_colour_manual(values = cbPalette) +
    theme_minimal()
  ggsave(paste0("02realization", i, ".pdf"),
         plot = figc_i, path = fig_path,
         width = w, height = h, units = "cm"
  )
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


plot2d <- function(observations,filename="NONE", save=FALSE){
  df2d <- data.frame(x=d$x, y=d$y, N_obs = observations)
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

for (i in 1:nsim){
  plot2d(result2d[,i], filename=paste0("posterior",i,".pdf"), save=TRUE)
}

