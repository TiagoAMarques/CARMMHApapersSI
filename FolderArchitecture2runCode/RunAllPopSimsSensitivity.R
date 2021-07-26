#this runs all the simulations and produces outputs for all species
source("Functions/reqfuns.R")         # auxiliary functions are here
source("Functions/SilerFuns.R")       # Siler model functions are here
source("Functions/runPopSims.R")      # the main simulation function

#define the species we run simulations for
Sps <- c("Bwsp", "Fatt", "Ggri", "Gmac", "Kosp", "Pele", "Pmac", "Satt", "Sbre", "Scly", "Scoe", "Sfro", "Slon", "Ttro", "Ttrs")

# ------------
# define the parameters that we evaluate sensitivity for
# "ASM" - age of sexual maturity - uses mean of Gamma if par not "a1r" in runSimsSensitivity
# "N0" - initial population size - uses mean of realizations if par not "N0" in runSimsSensitivity
# "per" - proportion exposed that recovers - uses mean of realizations if par not "per" in runSimsSensitivity
# "Fmax" - Maximum fecundity - uses mean of realizations if par not "Fmax" in runSimsSensitivity
# "Fnom" - Nominal Fecundity - uses mean of realizations if par not "Fnom" in runSimsSensitivity
# "rho" - density dependence parameter - uses mean of realizations if par not "rho" in runSimsSensitivity
# "br" - baseline reproduction - uses mean of realizations if par not "br" in runSimsSensitivity
# "por" - post spill reproductive success rate - uses mean of realizations if par not "por" in runSimsSensitivity
# "spos" - post oil spill survival  - uses mean of realizations if par not "spos" in runSimsSensitivity (SR for non Ttru)
# "ascS" -  age sex class survival -  uses mean of realizations of survival curves if par not "ascS" in runSimsSensitivity
# "PM" - proportion marked  -  uses mean of realizations of P(marked|age) if par not "PM" in runSimsSensitivity
# "pe" - proportion exposed - uses mean of realizations if par not "pe" in runSimsSensitivity
# ---------------
# "BrTt" - Baseline reproductive success rate for Ttru
# "PorTt" - post spill reproductive success rate for Ttru
# "SR" - survival reduction -  uses mean of realizations if par not "SR" in runSimsSensitivity
# ---------------

#list of parameters
parS <- c("ASM","N0","pe","per","Fmax","Fnom","rho","br","por","spos","ascS","PM","SR","BrTt","PorTt")


#define the number of iterations to run
ns <- 2
# define the number of years each iteration is run for
ny <- 75
# For each species, note Ttru comes first
for (i in Sps){
  # For each parameter we want to evaluate sensitivity for, see list above  
  for (j in parS) {
      runPopSims(Sp = i, type = "Sens", nsims = ns, nyears = ny, par =j)
    }
}