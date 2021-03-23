#this runs all the simulations and produces outputs for all species
source("Functions/reqfuns.R")         # most functions are here
source("Functions/SilerFuns.R")       # Siler model functions are here
source("Functions/runSimsSensitivity.R")         # To run the simulations for all species
#define the species we run simulations for
Sps <- c("Ttru","Bwsp", "Fatt", "Ggri", "Gmac", "Kosp", "Pele", "Pmac", "Satt", "Sbre", "Scly", "Scoe", "Sfro", "Slon", "Ttro", "Ttrs")
# define the parameters that we evaluate sensitivity for
# "a1r" - age at first reproduction - uses mean of Gamma if par not "a1r" in runSimsSensitivity
# "N0" - initial population size - uses mean of realizations if par not "N0" in runSimsSensitivity
# "pe" - proportion exposed - uses mean of realizations if par not "N0" in runSimsSensitivity
# "per" - proportion exposed that recovers - uses mean of realizations if par not "per" in runSimsSensitivity
# "Fmax" - Maximum fecundity - uses mean of realizations if par not "Fmax" in runSimsSensitivity
# "Fnom" - Nominal Fecundity - uses mean of realizations if par not "Fnom" in runSimsSensitivity
# "rho" - density dependence parameter - uses mean of realizations if par not "rho" in runSimsSensitivity
# "br" - baseline reproduction - uses mean of realizations if par not "br" in runSimsSensitivity
# "por" - post spill reproductive success rate - uses mean of realizations if par not "por" in runSimsSensitivity
# "spos" - post oil spill survival  - uses mean of realizations if par not "spos" in runSimsSensitivity
# Next pars are ONLY RELEVANT FOR NON T TRU results - see if condition below
# "BrTt"- Baseline reproductive success rate for Ttru
# "PorTt" - post spill reproductive success rate for Ttru
# what about sensitivity to N0BB ???
parS <- c("a1r","N0","pe","per","Fmax","Fnom","rho","br","por","BrTt","PorTt")
#define the number of iterations to run
ns <- 3
# define the number of years each iteration is run for
ny <- 75
p<-1
for (j in parS) {
for (i in Sps){
  if(!(j %in% c("BrTt","PorTt")) | i !="Ttru") runSimsSensitivity(Sp = i, nsims = ns, nyears = ny, par =j)
print(p)
p<-p+1
  }
}