#This runs the BB BND low salinity simulations
#Saves results to out.Rdata, with objects res.pop (population size per year)
# and res.mort (yearly total mortality)

source("Functions/runPopSims.R")         # To run the simulations for all species

#define strata
nstr <- 4
strata <- c("Island", "Southeast", "Central", "West")

#define the number of iterations to run
ns <- 1000
# define the number of years each iteration is run for
ny <- 75 #67 #2027 + 50 - 1 = 2076.  2076 - 2010 + 1 = 67

seed <- 3299728 #7134672

#create arrays to hold results
res.pop <- array(NA, dim = c(ny, ns, nstr, 3), 
             dimnames = list(year = 1:ny, sim = 1:ns, 
                             stratum = strata, reality = c("diversion", "no_diversion", "no_dwh")))
res.mort <- array(NA, dim = c(ny - 1, ns, nstr, 3), 
             dimnames = list(year = 1:(ny-1), sim = 1:ns, 
                             stratum = strata, reality = c("diversion", "no_diversion", "no_dwh")))
#vector of ages/stages used to calculate annual mortality
mort.vec <- c(1:60, 1:60 + 61 * 1, 1:60 + 61 * 2, 1:60 + 61 * 3, 1:60 + 61 * 4, 1:60 + 61 * 5)
for(i in 1:nstr){
  
  #Low salinity scenario
  runPopSims(Sp = "LowSal", stratum = i, nsims = ns, nyears = ny, seed = seed)
  load(paste0("InOutBySp/LowSal/LowSalsimres",ns,"Sim.RData"))
  #Work out total population size
  res.pop[, , i, 1] <- apply(simres[, , , 1], c(2, 3), sum)
  #Work out number of deaths
  tmp.mort <- simres[mort.vec, 1:(ny - 1), , 1] - simres[mort.vec + 1, 2:ny, , 1]
  res.mort[, , i, 1] <- apply(tmp.mort, c(2, 3), sum)
  
  #Baseline scenario
  runPopSims(Sp = "Ttru", stratum = i, nsims = ns, nyears = ny, seed = seed)
  load(paste0("InOutBySp/Bottlenose_dolphin_BB/Ttrusimres",ns,"Sim.RData"))
  res.pop[, , i, 2] <- apply(simres[, , , 1], c(2, 3), sum)
  tmp.mort <- simres[mort.vec, 1:(ny - 1), , 1] - simres[mort.vec + 1, 2:ny, , 1]
  res.mort[, , i, 2] <- apply(tmp.mort, c(2, 3), sum)
  #Also record numbers with no DWH 
  res.pop[, , i, 3] <- apply(simres[, , , 2], c(2, 3), sum)
  tmp.mort <- simres[mort.vec, 1:(ny - 1), , 2] - simres[mort.vec + 1, 2:ny, , 2]
  res.mort[, , i, 3] <- apply(tmp.mort, c(2, 3), sum)
}
save(res.pop, res.mort, file = "InOutBySp/LowSal/out.Rdata")
