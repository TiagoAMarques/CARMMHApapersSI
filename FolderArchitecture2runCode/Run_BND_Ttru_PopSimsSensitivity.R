#this runs all the simulations and produces outputs for BND Ttru
source("Functions/reqfuns.R")         # auxiliary functions are here
source("Functions/SilerFuns.R")       # Siler model functions are here
source("Functions/runPopSims.R")      # the main simulation function

# ------------
# define the parameters that we evaluate sensitivity for
# "a1r" - age at first reproduction - uses mean of Gamma if par not "a1r" in runSimsSensitivity
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
# ---------------
# While pe is always 1 for Ttru, it is nonetheless evaluated, as a control: no variability is expected
# "pe" - proportion exposed - uses mean of realizations if par not "N0" in runSimsSensitivity
# ---------------

#list of parameters
parS <- c("a1r","N0","pe","per","Fmax","Fnom","rho","br","por","spos","ascS","PM")

#define the number of iterations to run
ns <- 500
# define the number of years each iteration is run for
ny <- 75
# For each parameter we want to evaluate sensitivity for, see list above  

for (j in parS) {
  runPopSims(Sp = "Ttru", type = "Sens", nsims = ns, nyears = ny, par =j)
}