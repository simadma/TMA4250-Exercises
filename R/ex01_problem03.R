# TMA4245 - Spatial statistics
# Project 1
# Authors: Mads Adrian Simonsen and Karina Lilleborge
# Exercise 1
# Problem 3

#libraries
library(MASS)
library(fields)
library(akima)
library(geoR)
library(gridExtra)
library(latex2exp)
library(tidyverse)
library(reshape2)

fig_path <- "./Figures/" #For saving figures
##a)
#distances from different points in the grid
distance <- 0:43
#regular grid 
grid <- expand.grid(seq(1,30,length.out = 100), seq(1,30,length.out=100)) #2D list of coordinates
names(grid)=c('x','y')
#gaussian random field, parameters
mu_r <- 0
sigma_r <- 2
xi_r <- 3



#random field #1
set.seed(1998)
field1 <- grf(n=1, grid=grid, cov.model = "exponential", cov.pars = c(sigma_r, xi_r))
#empirical variogram based on realisation
empirical_variogram1 <- variog(field1)


#random field #2
set.seed(1997)
field2 <- grf(n=1, grid=grid, cov.model = "exponential", cov.pars = c(sigma_r, xi_r))
#empirical variogram based on realisation
empirical_variogram2 <- variog(field2)


#random field #3
set.seed(1996)
field3 <- grf(n=1, grid=grid, cov.model = "exponential", cov.pars = c(sigma_r, xi_r))
#empirical variogram based on realisation
empirical_variogram3 <- variog(field3)

image(field1, col=grey.colors(100))
#plots for the gaussian random fields
pdf(paste0(fig_path,"grfs.pdf"), width = 6, height = 2) #open pdf
op <- par(mfrow=c(1,3)) #initiate subplotting
image(field1, col=grey.colors(100,start=0,end=1)) #,  plot.title = {par(cex.main=1);title(main = "Realization 1",xlab = "x", ylab = "y")})
image(field2, col=grey.colors(100,start=0,end=1)) #, plot.title = {par(cex.main=1);title(main = "Realization 2",xlab = "x", ylab = "y")})
image(field3, col=grey.colors(100,start=0,end=1)) #, plot.title = {par(cex.main=1);title(main = "Realization 3",xlab = "x", ylab = "y")})
dev.off() #close pdf


#plots for variograms
pdf(paste0(fig_path,"empiricalvariograms3.pdf"), width = 6, height=2) #open pdf
op <- par(mfrow=c(1,3)) #initiate subplotting

plot(empirical_variogram1, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)

plot(empirical_variogram2, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)

plot(empirical_variogram3, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)
dev.off() #close pdf

## d) Generate 36 random locations on the grid
#Use a integer grid 
grid <- expand.grid(1:30,1:30)
set.seed(1996)
field <- grf(n=1, grid=grid, cov.model = "exponential", cov.pars = c(sigma_r, xi_r)) #new realization with new grid
set.seed(2021)
rpoints36 <- sample(1:900, 36)

#get these gridpoints 
rlocations36 <- grid[rpoints36,]


#get the data
data36 <- field$data[rpoints36]

#Empirical variance from the 36 randomly located grid points
emp_variogram36 <- variog(coords=rlocations36, data=data36)

pdf(paste0(fig_path,"empiricalvariogram36points.pdf"), width = 4*2, height = 4*1.5)
plot(emp_variogram36, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)
dev.off()

## e) Repeat d) with 9, 64 and 100 points
#first get random indexes from the grid
set.seed(2022)
rpoints9 <- sample(1:900, 9)
set.seed(2023)
rpoints64 <- sample(1:900, 64)
set.seed(2024)
rpoints100 <- sample(1:900, 100)


#get the grid points from these indexes
rlocations9 <- grid[rpoints9,]
rlocations64 <- grid[rpoints64,]
rlocations100 <- grid[rpoints100,]

#get the data for the random points
data9 <- field$data[rpoints9]
data64 <- field$data[rpoints64]
data100 <- field$data[rpoints100]

#Empirical variogram from the 9, 64 and 100 randomly located grid points
emp_variogram9 <- variog(coords=rlocations9, data=data9)
emp_variogram64 <- variog(coords=rlocations64, data=data64)
emp_variogram100 <- variog(coords=rlocations100, data=data100)

pdf(paste0(fig_path,"empiricalvariograms9-64-100.pdf"), width=6, height=2) #open pdf
op <- par(mfrow=c(1,3)) #start subplotting
plot(emp_variogram9, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)

plot(emp_variogram64, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)

plot(emp_variogram100, type="l")
lines(distance, sigma_r*(1-exp(-1*distance/xi_r)), lty=2)
dev.off() #close pdf

#loglikelihood estimation on sigma and xi
#empirical variogram based on realization
empirical_variogram1 <- variog(field)
#do estimation on sigma and xi based on the 36 points
ml_red36 <- likfit(coords=rlocations36, data=data36, ini=c(1,1.1), cov.model = "exponential")

#based on the 9 points
ml_red9 <- likfit(coords=rlocations9, data=data9, ini=c(1,1.1), cov.model = "exponential")

#based on the 64 points
ml_red64 <- likfit(coords=rlocations64, data=data64, ini=c(1,1.1), cov.model = "exponential")


#based on the 100 points
ml_red100 <- likfit(coords=rlocations100, data=data100, ini=c(1,1.1), cov.model = "exponential")


df.ml <- data.frame(xi.hat = c(xi_r, ml_red9$phi, ml_red36$phi, ml_red64$phi, ml_red100$phi), sigmasq.hat = c(sigma_r,ml_red9$sigmasq, ml_red36$sigmasq, ml_red64$sigmasq, ml_red100$sigmasq))


#Compare estimated variogram to exact
pow.exp.variogram <- function(distance, sigma, xi){
  return (sigma*(1-exp(-1*distance/xi)))
}
exact.vario = pow.exp.variogram(distance, sigma = df.ml$sigmasq.hat[1], xi=df.ml$xi.hat[1])
est.vario9 = pow.exp.variogram(distance, sigma = df.ml$sigmasq.hat[2],xi=df.ml$xi.hat[2])
est.vario36 = pow.exp.variogram(distance, sigma = df.ml$sigmasq.hat[3],xi=df.ml$xi.hat[3])
est.vario64 = pow.exp.variogram(distance, sigma = df.ml$sigmasq.hat[4],xi=df.ml$xi.hat[4])
est.vario100 = pow.exp.variogram(distance, sigma = df.ml$sigmasq.hat[5],xi=df.ml$xi.hat[5])

df.est.variogram <- data.frame(distance=distance, exact=exact.vario, est9=est.vario9, est36=est.vario36, est64=est.vario64, est100=est.vario100)

ggplot(data=df.est.variogram) +
  geom_line(aes(x=distance, y=exact, color = "Exact function")) +
  xlab(TeX("Distance, $\\tau$")) + 
  ylab(TeX("Variogram, $\\gamma$"))+
  geom_line(aes(x=distance, y=est9, color="9 points"))+
  geom_line(aes(x=distance, y=est36, color="36 points"))+
  geom_line(aes(x=distance, y=est64, color="64 points"))+
  geom_line(aes(x=distance, y=est100, color="100 points"))+
  labs(color="") + 
  theme_bw() +
  theme(legend.position = "right")
  
pdf(paste0(fig_path,"parameterestimation.pdf"), width=6, height=4)
dev.off()

