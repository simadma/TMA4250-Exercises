library(MASS)
library(fields)
library(akima)
library(geoR)


set.seed(1995)
# library

# loading data file

a <- read.table('topo')

# cov.spatial() geoR

distance <- 1:10
grid <- expand.grid(seq(0, n, by = 1), seq(0, 20, by = 1))
response <- rnorm(6636, 5, 0.2)
cov.spatial(distance, cov.model = "powered.exponential", cov.pars=c(1, 5), kappa = 1.5) # power, matern etx cov.spars
#sigmasq=1 matern gir samme som over

# interp() - gridded
data(akima)
# interp()
# asp - shrink
image.plot(interp(akima$x, akima$y, akima$z), asp = 2)
contour


# kriging Oppgave2
temp <- as.geodata(a)
#trend.d trend.l må vere like
krige.conv()# spesifiser sigmasq, og phi
krige.control()



#likfit
#variog bruk coords og data i stede for geodata(kanskje)

# uttrykkene står i dokumentasjonen
?cov.spatial

#plotting
