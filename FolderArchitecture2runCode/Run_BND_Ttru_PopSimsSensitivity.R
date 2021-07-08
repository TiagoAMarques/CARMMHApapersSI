#this runs all the simulations and produces outputs for BND Ttru
source("Functions/reqfuns.R")         # auxiliary functions are here
source("Functions/SilerFuns.R")       # Siler model functions are here
source("Functions/runPopSims.R")      # the main simulation function

# ------------
# define the parameters that we evaluate sensitivity for
# "ASM" - age at sexual maturity - uses mean of Gamma if par not "ASM" in runPopSims
# "N0" - initial population size - uses mean of realizations if par not "N0" in runPopSims
# "per" - proportion exposed that recovers - uses mean of realizations if par not "per" in runPopSims
# "Fmax" - Maximum fecundity - uses mean of distribution if par not "Fmax" in runPopSims
# "Fnom" - Nominal Fecundity - uses mean of distribution if par not "Fnom" in runPopSims
# "rho" - density dependence parameter - uses mean of distribution if par not "rho" in runPopSims
# "br" - baseline reproduction - uses mean of distribution if par not "br" in runPopSims
# "por" - post spill reproductive success rate - uses mean of distribution if par not "por" in runPopSims
# "spos" - post oil spill survival  - uses mean of distribution if par not "spos" in runPopSims
# "ascS" -  age sex class survival -  uses mean of realizations of survival curves if par not "ascS" in runPopSims
# "PM" - proportion marked  -  uses mean of realizations of P(marked|age) if par not "PM" in runPopSims
# ---------------
# While pe is always 1 for Ttru, it is nonetheless evaluated, as a control: no variability is expected
# "pe" - proportion exposed - uses mean of realizations if par not "pe" in runPopSims
# ---------------

#list of parameters
parS <- c("ASM","N0","pe","per","Fmax","Fnom","rho","br","por","spos","ascS","PM")

#define the number of iterations to run
ns <- 1000
# define the number of years each iteration is run for
ny <- 75
# For each parameter we want to evaluate sensitivity for, see list above  

for (j in parS) {
  runPopSims(Sp = "Ttru", type = "Sens", nsims = ns, nyears = ny, par =j)
}
