#this runs all the simulations and produces outputs for all species
source("Functions/runPopSims.R")         # To run the simulations for all species
#define the species we run simulations for

#define the species we run simulations for
Sps <- c("Ttru","Bwsp", "Fatt", "Ggri", "Gmac", "Kosp", "Pele", "Pmac", "Satt", "Sbre", "Scly", "Scoe", "Sfro", "Slon", "Ttro", "Ttrs")

#define the number of iterations to run
ns <- 9999
# define the number of years each iteration is run for
ny <- 75
for (i in Sps){
  runPopSims(Sp = i, nsims = ns, nyears = ny, seed = 7134672)
}
