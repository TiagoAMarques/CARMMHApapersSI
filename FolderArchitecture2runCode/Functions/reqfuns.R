# Author: TAM

# This file includes several required functions

# 1. getbeta    - computes the density dependent fecundity beta paramter
# 
# 2. ft         - computes the fecundity as a function of population size and other 
#                 parameters) depends on "getbeta()"
#                 
# 3. getRedFac  - computes the reduction factors; works for either fecundity or survival
# 
# 4. getM       - computes the transition matrix requires: all that affects 
#                 survival and reproduction
#                 
# 5. plotSims   - a plotting funtion, that plots the results coming outof running 
#                 the code inside e.g. Schwackeetal2017.R or CodeComparison.R
#
# 6. getSpData  - a function that allows you to read species details from a file
#                 definition list
#
# 7. plotIP     - plot inital population distribution, against the population exposed
#                 or the proportion exposed
#
# 8. getStableF0- get the F0 that leads to a constant population
#
# 9. getInjury  - using the object produced by a simulation run, compute the injury metrics 
#
# 10. getfilename - a bespoke function to make code shorter when defining file names used inside runPopSims
#
# 11. getBetaDistPars - returns beta distribution alpha and beta parameters given input p and var(p)

# 1. Function getbeta -------------------------------------
getbeta <- function(Fmax, rho, Fnom, Nnom) {
  
  # get the density dependent fecundity beta parameter
  # implements equation 8 in Schwacke et al 2017
  beta <- 1 / Nnom * ((Fmax - Fnom) / Fnom) ^ (1 / rho)
  return(beta)
}


#---- 2. Function ft ------------------------------------------
ft <- function(Nt, Fmax, rho, Nnom, Fnom) {
  
  # get fecundity given population size (and the other relevant parameters)
  # implements equation 7 in Schwacke et al 2017
  beta <- getbeta(Fmax, rho, Fnom, Nnom)
  res  <- Fmax / (1 + (beta * Nt)^rho)
  return(res)
}


# 3. Function getRedFac -----------------------------------
getRedFac <- function(y1, ny2r2n = 15, ny = 150) {
  # this function returns a yearly reduction factor for either fecundity or survival
  # Inputs:
  #       y1     - value in the first year post spill
  #       ny2r2n - number of years to return to normal (i.e., baseline) 
  #       ny     - total number of years to get reduction factor for
  redfac <- c(y1, y1 + 1:ny2r2n * (1 - y1) / ny2r2n, rep(1, ny - ny2r2n))
  return(redfac)
}


# 4. Function getM ----------------------------------------

# Creating the transition matrix
# inputs required: all that affects fecundity and survival
getM <- function(na, dimm, femalesur, malesur, srf, srf.yr1 = srf, 
                 ddfr, frf, frf.yr1 = frf, a1stR, alasR) {
  
  # given all the required inputs, this function structures the data to 
  # get a suitable transition matrix
  # 
  # Inputs:
  # na   - number of age classes
  # dimm - number of classes
  # 
  # Survival:
  # femalesur - female survival, per age class
  # malesur   - male survival, per age class
  # srf       - survival reduction factor
  # srf.yr1   - survival reduction factor in yr1 - for animals that don't recover
  # 
  # Fecundity:
  # ddfr    - density dependent fecundity rate
  # frf     - fecundity reduction factor
  # frf.yr1 - fecundity reduction factor in yr1 - for animals that don't recover
  # a1stR   - the age of first reproduction 8 means that females with 0 to 7
  #           years old do not reproduce
  # alasR   - the age of last reproduction
  #           48 means females 49 years old and older do not reproduce
  
  # FECUNDITY -------------------------
  # Transition probability
  TransM <- matrix(0, nrow = dimm, ncol = dimm)
  
  # values for fecundity --------------

  ddfrv <- rep(ddfr, na) # vector of unexposed females fecundities
  ddfrv[1:a1stR] <- 0 # females younger than a1stR have fecundity 0
  # exposed that recover - impact by fecundity reduction factor
  ddfrver <- ddfrv * frf
  # exposed that do not recover - impact by yr 1 fecundity reduction factor
  ddfrvenr <- ddfrv * frf.yr1
  
  # 1st row - produces females at time 0
  #unexposed animals
  ffecunfem   <- ddfrv / 2      # females
  ffecunma    <- rep(0, na)     # males fecundity - always 0
  # exposed that recover
  ffecexpfemar <- ddfrver / 2   # females
  ffecexpmar   <- rep(0, na)    # males fecundity - always 0
  # exposed that do not recover
  ffecexpfemanr <- ddfrvenr / 2 # females
  ffecexpmanr   <- rep(0, na)   # males fecundity - always 0
  
  TransM[1, ] <- c(ffecunfem, ffecunma, ffecexpfemar, ffecexpmar, ffecexpfemanr, ffecexpmanr)
  
  # na + 1 row - produces males at age 0 - just the same as above, but useful to have
  # here if variability in sex ratio is included later
  mfecunfem    <- ddfrv / 2  # unexposed females
  mfecunma     <- rep(0, na) # males fecundity - always 0
  mfecexpfemar  <- ddfrver / 2 # exposed
  mfecexpmar    <- rep(0, na) # males fecundity - always 0
  mfecexpfemanr  <- ddfrvenr / 2 # exposed
  mfecexpmanr    <- rep(0, na) # males fecundity - always 0
  
  TransM[na + 1, ] <- c(mfecunfem, mfecunma, mfecexpfemar, mfecexpmar, mfecexpfemanr, mfecexpmanr)
  
  # values for survival ---------------
  # constrained so all values less than or equal to 1
  # values > 1 can be produced otherwise in low salinity simulations
  #unexposed animals
  sfus <- femalesur       # females 
  smus <- malesur         # males 
  sfus[sfus > 1] <- 1; smus[smus > 1] <- 1
  # exposed that recover
  sfers <- femalesur * srf # females 
  smers <- malesur * srf   # males 
  sfers[sfers > 1] <- 1; smers[smers > 1] <- 1
  # exposed that do not recover
  sfenrs <- femalesur * srf.yr1 # females 
  smenrs <- malesur * srf.yr1 # males
  sfenrs[sfenrs > 1] <- 1; smenrs[smenrs > 1] <- 1
  
  # construct transition matrix
  for(m in 1:(na - 1)) {
    # females unexposed
    TransM[m + 1, m] <- sfus[m]
    # males unexposed
    TransM[na + 1 + m, na + m] <- smus[m]
    # females exposed that recover
    TransM[na + na + 1 + m, na + na + m] <- sfers[m]
    # males exposed that recover
    TransM[na + na + na + 1 + m, na + na + na + m] <- smers[m]
    # females exposed that recover
    TransM[na + na + na + na + 1 + m, na + na + na + na + m] <- sfenrs[m]
    # males exposed that recover
    TransM[na + na + na + na +na + 1 + m, na + na + na + na + na + m] <- smenrs[m]
  }
  return(TransM)
}


# 5. Function plotSims ------------------------------------
plotSims <- function(sr) {
  # ploting simulation results just a wrapper for code inside plots4sims.R
  # will break if anything changes, can easily be made more general
  # evolution of population size over time
  # 150 years - oil spill
  nyears <- dim(sr)[2]
  nsims <- dim(sr)[3]
  ylims <- c(min(c(colSums(sr[, , 1:nsims, 2]), colSums(sr[, , 1:nsims, 1])))*0.9, max(c(colSums(sr[, , 1:nsims, 2]), colSums(sr[, , 1:nsims, 1])))*1.1)
  par(mfrow = c(2, 2), mar = c(3, 4, 1, 0.5))
  # just first 50 years - NO oil spill
  plot(colSums(sr[, , 1, 2]), type = "n", ylim = ylims, xlim = c(0, 50), ylab = "Predicted population size", xlab = "", las = 1)
  # grid(lty=1)
  poptraj = matrix(NA, ncol = nyears, nrow = nsims)
  for(i in 1:nsims) {
    lines(colSums(sr[, , i, 2]), type = "l", col = rgb(0, 0, 0, 0.15), lwd = 0.7)
    poptraj[i, ] = colSums(sr[, , i, 2], na.rm = FALSE)
  }
  lines(colMeans(poptraj, na.rm = TRUE), type = "l", lwd = 3, col = "#1b2ac8")
  
  
  # just first 50 years - oil spill
  par(mar = c(3, 2, 1, 2.5))
  
  plot(colSums(sr[, , 1, 1]), type = "n", ylim = ylims, xlim = c(0, 50), yaxt = "n", ylab = "", xlab = "")
  axis(2, seq(0, 6000, 1000), labels = F)
  poptraj <- matrix(NA, ncol = nyears, nrow = nsims)
  for(i in 1:nsims) {
    lines(colSums(sr[, , i, 1]), type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
    poptraj[i, ] = colSums(sr[, , i, 1])
  }
  lines(colMeans(poptraj, na.rm = TRUE), type = "l", lwd = 3, col = "#C81B2A")
  
  
  # 150 years - NO oil spill
  par(mar = c(3, 4, 0.5, 0.5))
  
  plot(colSums(sr[, , 1, 2]), type = "n", ylim = ylims, xlab = "Years post spill", ylab = "Predicted population size", las = 1)
  title(xlab = "Years post-spill", line = 2)
  poptraj = matrix(NA, ncol = nyears, nrow = nsims)
  for(i in 1:nsims) {
    lines(colSums(sr[, , i, 2]), type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
    poptraj[i, ] = colSums(sr[, , i, 2])
  }
  lines(colMeans(poptraj, na.rm = TRUE), type = "l", lwd = 3, col = "#1b2ac8")
  
  
  par(mar = c(3, 2, 0.5, 2.5))
  
  plot(colSums(sr[, , 1, 1]), type = "n", yaxt = "n", ylim = ylims, xlab = "Years post spill", ylab = "")
  axis(2, seq(0, 6000, 1000), labels = F)
  title(xlab = "Years post-spill", line = 2)
  poptraj = matrix(NA, ncol = nyears, nrow = nsims)
  for(i in 1:nsims) {
    lines(colSums(sr[, , i, 1]), type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
    poptraj[i, ] = colSums(sr[, , i, 1])
  }
  #abline(h = 1000, col = "#2ac81b", lty = 2, lwd = 2)
  lines(colMeans(poptraj, na.rm = TRUE), type = "l", lwd = 3, col = "#C81B2A")
}


# 6. Function getSpData -------------------------------------
#this function reads in the information in a file
#to create a species object which can then be used 
#for a given simulation

#the file must have a specific fixed format
#and in this case was created by TAM
getSpData <- function(sp, file = "InputFiles/SpeciesDefinitionFile.xlsx"){
  # the 1st argument for this function is the species code
  # sp: a 4 letter character acronym
  # e.g. Pmac for sperm whale
  #if needed check file for acronym list 
  #(acronym is the first column on the file)
  # The second argument is the file with the species information
  require(readxl)
  SDF <- read_excel(file, na = "NA")
  Spinfo = SDF[SDF$Species==sp, ]
  Sp <- list(
    #Folder for reading and storing results
    folder = Spinfo$folder,
    #Total population size
    N = Spinfo$N,
    sdN = Spinfo$sdN,
    #maximum fecundity paramters
    minFmax = Spinfo$minFmax,
    modeFmax = Spinfo$modeFmax,
    maxFmax = Spinfo$maxFmax,
    #minimum fecundity parameters
    minFnom = Spinfo$minFnom,
    modeFnom = Spinfo$modeFnom,
    maxFnom = Spinfo$maxFnom,
    #Rho parameters
    rhoshift = Spinfo$rhoshift,
    rhoshape = Spinfo$rhoshape,
    rhoscale = Spinfo$rhoscale,
    #minimum inter-birth interval  
    ib.min = Spinfo$ibmin,
    #inter-birth interval resulting in a stable population
    ib.nom = Spinfo$ibnom,
    #f.max <- 1 / ib.min
    #f.nom <- 1 / ib.nom
    #proportion of the population exposed to oil
    pexposed = Spinfo$pe,  
    #Density dependence rho parameter
    #model to influence interbrith intervals
    rho = Spinfo$rho,
    #stage specific survival rates
    #survs = as.numeric(Spinfo[,(3+15+1):(3+15+1+15)]),
    #expected years in stage
    #T.init = as.numeric(Spinfo[,3:(3+15)]),
    #Gestation duration
    gd = Spinfo$gd,
    #proportion of exposed animals that recover
    prl = Spinfo$prl,
    pru = Spinfo$pru,
    pra = Spinfo$pra,
    prb = Spinfo$prb,
    #Survival reduction factor
    srl = Spinfo$srl,
    sru = Spinfo$sru,
    sra = Spinfo$sra,
    srb = Spinfo$srb,
    #added mortality due to oil
    mort.effect = Spinfo$me,
    #years it remains constant
    mc = Spinfo$mc,
    #years it takes to get back to baseline after constant
    yrbs = Spinfo$yrbs,
    #baseline reproduction
    meanbrs = Spinfo$meanbrs,
    sdbrs = Spinfo$sdbrs,
    meanpors = Spinfo$meanpors,
    sdpors = Spinfo$sdpors,
    #age at first reproduction
    sha1str = Spinfo$sha1str,
    sca1str = Spinfo$sca1str,
    #age at last reproduction
    alastRep = Spinfo$alastRep,
    #lowered reproduction due to oil
    repro.effect = Spinfo$re,
    #years it remains constant
    rc = Spinfo$rc,
    #years it takes to get back to baseline after constant
    rc2b = Spinfo$rc2b,
    # baseline survival, BND only
    absr = Spinfo$absr,
    bbsr = Spinfo$bbsr,
    #quantities from MTTIQ document
    LCY = Spinfo$LCY,
    MPD = Spinfo$MPD,
    YTR = Spinfo$YTR
  )
  return(Sp)
}

# 7. Function plotIP -------------------------------------
plotIP <- function(IP, EP, type = 1){
  #inputs
  #IP - a vector with the initial population size
  #EP - a vector with the exposed animals
  #type: sets the type of plot
  #       1: pop size versus exposed pop size
  #       2: pop size versus proportion exposed
  require(ggplot2)
  require(gridExtra)
  if(type==1){  
    hist_top <- ggplot()+geom_histogram(aes(IP))+xlab("Total stock")
    empty <- ggplot()+geom_point(aes(1, 1), colour = "white")+
      theme(axis.ticks = element_blank(), 
            panel.background = element_blank(), 
            axis.text.x = element_blank(), axis.text.y = element_blank(),           
            axis.title.x = element_blank(), axis.title.y = element_blank())
    
    scatter <- ggplot()+geom_point(aes(IP, EP))
    hist_right <- ggplot()+geom_histogram(aes(EP))+coord_flip()+xlab("Exposed animals")
    grid.arrange(hist_top, empty, scatter, hist_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  }
  if(type==2){  
    pExp <- EP/IP
    hist_top <- ggplot()+geom_histogram(aes(IP))+xlab("Total stock")
    empty <- ggplot()+geom_point(aes(1, 1), colour = "white")+
      theme(axis.ticks = element_blank(), 
            panel.background = element_blank(), 
            axis.text.x = element_blank(), axis.text.y = element_blank(),           
            axis.title.x = element_blank(), axis.title.y = element_blank())
    
    scatter <- ggplot()+geom_point(aes(IP, pExp))
    hist_right <- ggplot()+geom_histogram(aes(pExp))+coord_flip()+xlab("Proportion Exposed animals")
    grid.arrange(hist_top, empty, scatter, hist_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  }  
}

# 8. Function getStableF0 -------------------------------------
getStableF0 <- function(F0, ns1, nc1, na1, femalesur1, malesur1, a1stRsim1){
  #a new function that takes the parameters in explicitly
  M0 <- getM(na        = na1, 
             dimm      = ns1 * nc1 * na1,
             femalesur = femalesur1,
             malesur   = malesur1,
             srf       = 1,
             ddfr      = F0,
             frf       = 1,
             a1stR     = a1stRsim1,
             alasR     = na1)
  target<-eigen.analysis(M0)$lambda1 - 1
  return(target)
}

# 9. Function getInjury  -------------------------------------

getInjury <- function(SI, p4YTR = 0.95, plot = TRUE, show.plot = TRUE, median = TRUE){
  # using the object produced by a simulation run, compute the injury metrics 
  # 1. Lost Cetacean Years, the difference between the baseline and injured 
  #       population sizes, summed over the entire modeled time period
  # 2. Years to recovery (YTR), the number of years required before the 
  #       injured population trajectory reaches 95% of the baseline population trajectory; 
  # 3. maximum proportional decrease (MPD), the difference between the 2 population 
  #       trajectories when the injured trajectory is at its lowest point, divided by the baseline
  #Inputs: 
  #       SI:         the simulation object
  #       p4YTR:      1-p4YTR is the proportion of the recovery required to assume recovered
  #       plot:       defaults to true, plot is produced
  #       show.plot:  defaults to true, plot is shown in R. 
  #                       If FALSE, plot is printed as file in corresponding Sp folder
  #       median:     defaults to true, if FALSE means are presented  
  #  
  #Outputs:
  #       A data.frame with   
  #       LCY
  #       YTR
  #       MPD  
  #       and (by default) an optional injury results plot is produced 
  nsims <- dim(SI)[3]
  nyears<- dim(SI)[2]
  LCYs <- YTRs <- MPDs<- numeric(nsims)  
  for(i in 1:nsims){
    #for each simulation
    #get the differences in N, per year
    diff.Ns <- colSums(SI[, , i, 2])-colSums(SI[, , i, 1])
    #and get the initial population size
    N0 <- sum(SI[, 1, i, 1])
    lim4YTR<-(1-p4YTR)*colSums(SI[, , i, 2])
    # Lost Cetacean Years
    LCYs[i] <- sum(diff.Ns)
    # this line is executed if something went wrong in the simulation
    # and there is no clear reason for it to happen but I saw it once 
    # during sensitivity analysis
    # when a population went extinct and the code would otherwise break
    if(is.nan(LCYs[i])) next 
    # Years to recovery
    p <- nyears
    while (diff.Ns[p]<lim4YTR[p]){
      p <- p-1
      #if p = 0 it means that even in year 0 the population did not go below the % difference
      # of (1-p4YTR) times the baseline population
      if(p==0) break
    }
    YTRs[i] <- p
    #maximum proportional decrease
    MPDs[i] <- max(diff.Ns/colSums(SI[, , i, 2]))*100
  }
  if(plot)
  {
    # Note to self: Not a great way of defining the path to save the image in
    # as it reads path from the workspace!
    if(show.plot==FALSE) png(filename = paste0("InOutBySp/", SpInfo$folder, "/", Sp, "Injury", nsims, ".png"), width = 3*480,  height = 2*480)
    par(mfrow = c(1, 3), mar = c(4, 3, 4, 1))
    mLCY <- round(mean(LCYs), 1)
    mdLCY <- round(median(LCYs), 1)
    sdLCY <- round(sd(LCYs), 1)
    hist(LCYs, main = paste0(ifelse(median, "Median", "Mean"), "LCY=", ifelse(median, mdLCY, mLCY), "(", sdLCY, ")"))
    abline(v = mLCY, lty = 2)
    abline(v = quantile(LCYs, probs = c(0.025, 0.5, 0.975),na.rm=TRUE), lty = 2, col = c(4, 3, 4))
    mYTR <- round(mean(YTRs), 1)
    sdYTR <- round(sd(YTRs), 1)
    mdYTR <- round(median(YTRs), 1)
    hist(YTRs, main = paste0(ifelse(median, "Median", "Mean"), "YTR=", ifelse(median, mdYTR, mYTR), "(", sdYTR, ")"))
    abline(v = mYTR, lty = 2)
    abline(v = quantile(YTRs, probs = c(0.025, 0.5, 0.975),na.rm=TRUE), lty = 2, col = c(4, 3, 4))
    mMPD <- round(mean(MPDs), 1)
    sdMPD <- round(sd(MPDs), 1)
    mdMPD <- round(median(MPDs), 1)
    hist(MPDs, main = paste0(ifelse(median, "Median", "Mean"), "MPD=", ifelse(median, mdMPD, mMPD), "(", sdMPD, ")"))
    abline(v = mMPD, lty = 2)
    abline(v = quantile(MPDs, probs = c(0.025, 0.5, 0.975),na.rm=TRUE), lty = 2, col = c(4, 3, 4))
    if(show.plot==FALSE)  dev.off()
  }
  return(data.frame(LCY = LCYs, YTR = YTRs, MPD = MPDs))
}

# 10. Function getfilename  -------------------------------------
getfilename <- function(a1, a2, a3, a4) {
  return(paste0(a1, a2, a3, a4))
}

# 11. 
getBetaDistPars <- function (p, var.p){
  #Returns the alpha and beta parameters from a beta distribution given input mean (p) and variance (var.p).
  #Inputs:
  # p - mean of beta distribution
  # var.p - variance of beta distribution
  #Outputs: - named list containing
  # alpha - alpha parameter
  # beta - beta parameter
  alpha <- ((1 - p) / var.p - 1 / p) * p^2
  beta <- alpha * (1 / p - 1)
  if (! ((alpha >= 0) & (beta >= 0))) stop ("Values for p and var.p cannot produce a valid beta distribution")
  return(list(alpha = alpha, beta = beta))  
}