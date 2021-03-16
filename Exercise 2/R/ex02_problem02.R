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
  geom_point() +
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
realizations <- matrix(0, nrow=nsim, ncol=n)
for (i in 1:nsim){
  realizations[i,] <- rpois(n=n, lambda=lambda_hat)
}


sim1 <- data.frame(x=d$x, y=d$y, N_obs=realizations[1,])
figc1<-ggplot(data=sim1, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Tree count")+
  scale_colour_manual(values = cbPalette) +
  theme_minimal()
ggsave("02realization1.pdf",
       plot = figc1, path = fig_path,
       width = w, height = h, units = "cm"
)
sim2 <- data.frame(x=d$x, y=d$y, N_obs=realizations[2,])
figc2<-ggplot(data=sim2, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Tree count")+
  scale_colour_manual(values = cbPalette) +
  theme_minimal()
ggsave("02realization2.pdf",
       plot = figc2, path = fig_path,
       width = w, height = h, units = "cm"
)
sim3 <- data.frame(x=d$x, y=d$y, N_obs=realizations[3,])
figc3<-ggplot(data=sim3, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Tree count")+
  scale_colour_manual(values = cbPalette) +
  theme_minimal()
ggsave("02realization3.pdf",
       plot = figc3, path = fig_path,
       width = w, height = h, units = "cm"
)
sim4 <- data.frame(x=d$x, y=d$y, N_obs=realizations[4,])
figc4<-ggplot(data=sim4, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Tree count")+
  scale_colour_manual(values = cbPalette) +
  theme_minimal()
ggsave("02realization4.pdf",
       plot = figc4, path = fig_path,
       width = w, height = h, units = "cm"
)
sim5 <- data.frame(x=d$x, y=d$y, N_obs=realizations[5,])
figc5<-ggplot(data=sim5, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Tree count")+
  scale_colour_manual(values = cbPalette) +
  theme_minimal()
ggsave("02realization5.pdf",
       plot = figc5, path = fig_path,
       width = w, height = h, units = "cm"
)
sim6 <- data.frame(x=d$x, y=d$y, N_obs=realizations[6,])
figc6<-ggplot(data=sim6, aes(x=x, y=y, color=as.factor(N_obs))) +
  geom_point() +
  labs(color="Tree count")+
  scale_colour_manual(values = cbPalette) +
  theme_minimal()
ggsave("02realization6.pdf",
       plot = figc6, path = fig_path,
       width = w, height = h, units = "cm"
)




