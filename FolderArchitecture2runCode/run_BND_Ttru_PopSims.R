#this runs all the simulations and produces outputs for Ttru BND species
source("Functions/runPopSims.R")        
#define the number of iterations to run
ns <- 9999
# define the number of years each iteration is run for
ny <- 75
runPopSims(Sp = "Ttru", nsims = ns, nyears = ny, seed = 7134672)

