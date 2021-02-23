#this runs all the simulations and produces outputs for all species
source("Functions/reqfuns.R")         # most functions are here
source("Functions/SilerFuns.R")       # Siler model functions are here
source("Functions/runSims.R")         # To run the simulations for all species
#define the species we run simulations for
Sps <- c("Bwsp", "Fatt", "Ggri", "Gmac", "Kosp", "Pele", "Pmac", "Satt", "Sbre", "Scly", "Scoe", "Sfro", "Slon", "Ttro", "Ttrs","Ttru")
#define the number of iterations to run
ns <- 2
# define the number of years each iteration is run for
ny <- 150
for (i in Sps){
runSims(Sp = i, nsims = ns, nyears = ny)
}