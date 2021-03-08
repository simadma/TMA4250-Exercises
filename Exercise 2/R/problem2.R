#Problem 2
#Authors: Mads Adrian Simonsen and Karina Lilleborge
#Date: 03-08-21

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
