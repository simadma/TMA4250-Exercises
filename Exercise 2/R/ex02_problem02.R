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

fig1 <- ggplot(data=d, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of observations")+
  labs(title="Observed number of pine trees")+
  theme_minimal()

fig1

fig2 <- ggplot(data=alpha) +
  geom_tile(aes(x=x, y=y, fill=alpha))+
  labs(title="The probability of observing pine trees")+
  theme_minimal()

fig2


#c) Estimation of lambda
area <- 100 #m^2

lambda_hat <- sum(d$N_obs)/(area*sum(alpha$alpha))
cat("Estimation of lambda_k is ", lambda_hat, "per square meter")

# Generate 6 realizations from the prior
nsim <- 6
n <- 900
realizations <- matrix(0, nrow=nsim, ncol=n)
for (i in 1:nsim){
  realizations[i,] <- rpois(n=n, lambda=lambda_hat)
}


sim1 <- data.frame(x=d$x, y=d$y, N_obs=realizations[1,])
figc1<-ggplot(data=sim1, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of trees")+
  labs(title="Number of pine trees")+
  theme_minimal()
sim2 <- data.frame(x=d$x, y=d$y, N_obs=realizations[2,])
figc2<-ggplot(data=sim2, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of trees")+
  labs(title="Number of pine trees")+
  theme_minimal()
sim3 <- data.frame(x=d$x, y=d$y, N_obs=realizations[3,])
figc3<-ggplot(data=sim3, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of trees")+
  labs(title="Number of pine trees")+
  theme_minimal()
sim4 <- data.frame(x=d$x, y=d$y, N_obs=realizations[4,])
figc4<-ggplot(data=sim4, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of trees")+
  labs(title="Number of pine trees")+
  theme_minimal()
sim5 <- data.frame(x=d$x, y=d$y, N_obs=realizations[5,])
figc5<-ggplot(data=sim5, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of trees")+
  labs(title="Number of pine trees")+
  theme_minimal()
sim6 <- data.frame(x=d$x, y=d$y, N_obs=realizations[6,])
figc6<-ggplot(data=sim6, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Number of trees")+
  labs(title="Number of pine trees")+
  theme_minimal()

pdf(file=paste0(fig_path,"priors6.pdf"), width=6, height=4)
op <- par(c(2,3))
figc1
figc2
figc3
figc4
figc5
figc6
dev.off()
