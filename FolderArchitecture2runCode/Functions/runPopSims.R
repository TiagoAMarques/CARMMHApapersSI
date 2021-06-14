runPopSims <- function(Sp, nsims, nyears, type = "Sim", parS = NULL, 
                       stratum = NULL, seed = NULL, verbose = TRUE) {
  # This function runs a set of simulations of the CARMMHA population dynamics model
  # for a given taxonomic unit Sp, a number of iterations nsims and for a number of years nyears
  # if under sensitivity mode (type!="Sim") all parameters at the mean nominal values
  # except for parameter parS that can be different things - see below for details
  # Requires
  #  R library progress
  # Inputs
  #       Sp:      the species code, see details below
  #       nsims:   the number of iterations to consider
  #       nyears:  the number of years to run each simulation for
  #       type:    one of 3 options
  #                 1. "Sim", the default which means it should simply run the simulations
  #                 2. "Sens" runs a sensitivity analysis (parameter parS varies, all other are kept fixed at mean)
  #                 3. "Elas" runs and elasticity analysis, used only for parameters that are not scalars, 
  #                     namely "ascS" or "PM", for that purpose renamed as "ascSE" and "PME", respectively  
  #       parS:    the parameter for which the simulations are being run for if under type==Sim, i.e. under sensitivity mode
  #       stratum: either NULL (default) or an integer from 1-4 (only valid for Ttru).  If an integer, runs on the 
  #                 proportion of the BSE BND population corresponding to that stratum number (see below). 
  #       seed:    if an integer value, then used to set the random number seed; if NULL (the default) no seed is set
  #       verbose: if TRUE (default) provides output to the console showing which iteration it is on
  #------------------------------------------------------------------------------------------------------------------
  # Details on "Sp" to define species to work with
  #--------------------------------------------------------------------------
  # Sp is a 4 letter character acronym (exceptions are Schwackeetal, LowSal)
  # See 'SpeciesDefinitionFile.xlsx'
  # e.g. 'Pmac' for sperm whale
  # if needed check file for acronym list correspondence to species
  # (acronym is the first column on the file 'SpeciesDefinitionFile.xlsx')
  #--------------------------------------------------------------------------
  # Details on stratum analyses
  #--------------------------------------------------------------------------
  # The stratum argument is only used for the dedicated analysis to evaluate the effect of increased freshwater
  # being diverted into Barataria Bay. That was an analysis conducted by the CARMMHA Team to respond to a request from the MMC
  # that might lead to a separate paper at some point.
  # Stratum is only valid for Sp = 'Ttru' or 'LowSal'.  In this case, the starting population
  # size is multiplied by the estimated proportion of the population in the corresponding stratum.
  # Stratum names are defined below.
  #--------------------------------------------------------------------------
  # Outputs
  #       There are several by-products of running the function (but no object returned at the end) 
  #       An R workspace is saved in a species specific folder
  #       The folder name is defined in "SpeciesDefinitionFile.xlsx"
  #       The object name reflects its contents
  #       The workspace contains the following objects
  #       simres: raw results, with population size by age, sex and class per year
  #               for each iteration of the simulation - see below details for object simres
  #       injury: object created by post-processing simres using function getInjury
  #               includes LCY, MPD and YTR for each iteration
  #       RunOnThe: time stamp for when the simulations were run
  #       TimeSpent: how long the simulations took to run
  #       pars2save: dataframe with parameter values used per iteration
  # Additionally, 3 sets of plots are produced and saved
  #         1. plot of population trajectories; 2. plot of injury metrics; 3. plot of fecundities at year 0
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # just a check to make sure argument type is well defined
  if (!(type %in% c("Sim", "Sens", "Elas"))) {
    stop("Wrong type used; only 'Sim', 'Sens' or 'Elas' allowed; please check argument type")
  }
  # check that if not under type="Sim" the parS argument is defined and not the NULL default 
  if (type != "Sim" & is.null(parS)){
    stop("Under sensitivity analysis, type='Sens' or type='Elas', parS cannot be NULL")
  }

  #--------------------------------------------------------------------------
  # Set an object to save parameter values used in each iteration
  pars2save <- data.frame(sim = 1:nsims, a1r = NA, N0 = NA, pe = NA, per = NA, Fmax = NA,
  Fnom = NA, rho = NA, br = NA, por = NA, spos = NA, ascS = NA, SR = NA, BrTt = NA, 
  PorTt = NA, y2R = NA, srTtru = NA, fecRed = NA, PM = NA)
  
  #--------------------------------------------------------------------------
  # to time the procedure
  start_time <- Sys.time()
  # set a random seed, to make results reproducible
  if (!is.null(seed)) set.seed(seed)
  
  #--------------------------------------------------------------------------
  # required libraries and functions
  #--------------------------------------------------------------------------
  # load required libraries
  library(mc2d) # for beta-pert distribution
  library(popbio) # for computing mean age and expected time to live
  if(verbose) require(progress) # progress bar library
  #--------------------------------------------------------------------------
  # load up required functions
  source("Functions/reqfuns.R") # most functions are here
  source("Functions/SilerFuns.R") # Siler model functions are here
  #--------------------------------------------------------------------------

  #--------------------------------------------------------------------------
  # Set counters and constants
  #--------------------------------------------------------------------------
  na <- 61 # number of age classes for BND
  ns <- 2 # number of sexes
  nc <- 3 # number of conditions (see details above)
  dimm <- ns * nc * na # number of classes (2 sexes * 3 conditions * number of age classes)
  multsE<-c(0.995,1,1.005) #this is only used in Elasticity mode for two parameters but the cost of setting here is negligible

  #Define all values of Sp that refer to BB BNDs.
  BB_BND_Sp <- c("Ttru", "Schwackeetal2017", "LowSal")
  
  #Define strata and check stratum being called appropriately
  strata <- c("Island", "Southeast", "Central", "West")
  if (is.null(stratum)){
    #Stratum must be specified for LowSal runs
    if(Sp == "LowSal")
      stop("Stratum must be specified for LowSal scenarios")
  } else {
    if(!(Sp %in% c("Ttru", "LowSal")))
      stop("Stratum analysis specified but species is not bottlenose dolphin")
    if(!(stratum %in% 1:4))
      stop("stratum argument must be NULL or an integer in the range 1-4")
  }
  #--------------------------------------------------------------------------

  #--------------------------------------------------------------------------
  # Get all the relevant species information required to run the simulation
  # taken from file "SpeciesDefinitionFile.xlsx"
  # This file is hardwired in the code to reside in the subfolder "InputFiles"
  if(Sp == "LowSal") {
    SpInfo <- getSpData("Ttru")
  } else {
    SpInfo <- getSpData(Sp)
  }
  #--------------------------------------------------------------------------

  #--------------------------------------------------------------------------
  # BND data is required even if this is not a BND simulation
  #--------------------------------------------------------------------------
  # Scaling factor with respect to BND
  #--------------------------------------------------------------------------
  # that happens for the gestation duration scaling but also
  # for values for which we only have for BND (like the Fmax distribution)
  Ttru <- getSpData("Ttru")
  # scaling for current species with respect to BND, to scale survivals
  # Scaling is by gestation duration (gd)
  scaling <- Ttru$gd / SpInfo$gd
  # matrix to hold results
  # might need to be bigger or smaller depending on scaling
  # store number of original age classes (the BB BND)
  naBND <- na
  # scale the number of age classes
  na <- ceiling(na / scaling)
  #reformat the dimensions of object given scaling
  dimm <- ns * nc * na
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # BB BND post oil Survival from Glennie et al 2021 SCR study
  #--------------------------------------------------------------------------
  # Survival
  SpostOilBB <- read.csv(file = "InOutBySp/Bottlenose_dolphin_BB/PostOilSurv.csv", header = FALSE, col.names = "S")
  SpostOilBB <- as.numeric(SpostOilBB$S)
  #forcing random drwas - there was an order in the draws
  SpostOilBB <- SpostOilBB[sample.int(n=length(SpostOilBB),replace = FALSE)]
  if (type != "Sim") {
    # if running a sensitivity analysis and parS is not "spos" just use the mean
    if(parS != "spos"){
      meanspos <- mean(SpostOilBB)
      SpostOilBB <- rep(meanspos,length(SpostOilBB))
    }
  }
  #--------------------------------------------------------------------------

  #--------------------------------------------------------------------------
  # Survival
  #--------------------------------------------------------------------------
  # reading in Siler parameters, for male and female BNDs
  # These files were provided by LS via email on the Fri 5/22/2020 4:10 PM
  pf <- read.csv("InputFiles/var_female 15May2020.csv")
  pm <- read.csv("InputFiles/var_male 15May2020.csv")
  # remove first column, which is just ID (not needed)
  pf <- pf[, -1]
  pm <- pm[, -1]
  # rename columns
  names(pf) <- names(pm) <- c("a1", "a2", "a3", "b1", "b3", "rmean")
  # check size of the two sets of samples is the same
  if(nrow(pf) != nrow(pm)) {
    stop ("Number of samples from Siler model for males and females do not match.")
  }
  # get a sample from the Siler model for each sim - which is a posterior distribution
  # number of Siler model draws
  nS <- nrow(pf)
  # sample from the Siler model to use for each of the sims
  # this is necessary to ensure that everything works even if
  # the number of simulations is different from the number of Siler posterior samples
  itS <- sample.int(nS, size = nsims, replace = TRUE)
  #--------------------------------------------------------------------------
  # Baseline survival
  #--------------------------------------------------------------------------
  # Obtaining P(marked|age) required for survival pre-oil spill
  # See SI file "SurvivalReduction" for a detailed description of this procedure
  # including the production of some of the objects used
  #--------------------------------------------------------------------------
  # get probability of marked as a function of age
  predictions <- read.csv(file = "InputFiles/pmarked.csv")
  predictions <- as.matrix(predictions)
  if(nrow(predictions) != nS) {
    # later code assumes the number of predictions is the same as the number of survival samples, 
    stop("Number of samples of probability marked does not equal the number of samples from Siler model.")
  }
  # for sensitivity only
  mpredictions <- colMeans(predictions)
  
  
  # get sex and age specific survivals
  #--------------------------------------------------------------------------
  # ages for Ttru
  agesBND <- 0:(naBND-1)
  # ages if scalling occurred
  ages <- 0:(na-1)

  # objects to hold female and male realizations - Ttru
  pxFs <- matrix(NA, nrow = length(agesBND), ncol = nS)
  pxMs <- matrix(NA, nrow = length(agesBND), ncol = nS)
  # objects to hold female and male realizations - scaled for species
  SpxFs <- matrix(NA, nrow = length(ages), ncol = nS)
  SpxMs <- matrix(NA, nrow = length(ages), ncol = nS)
  # for each observation of the posterior
  for (i in 1:nS) {
    # Survival Ttru
    pxFs[, i] <- px(agesBND, pf[i, 1], pf[i, 2], pf[i, 3], pf[i, 4], pf[i, 5])
    pxMs[, i] <- px(agesBND, pm[i, 1], pm[i, 2], pm[i, 3], pm[i, 4], pm[i, 5])
    # Scaled survival
    SpxFs[, i] <- px(ages * scaling, pf[i, 1], pf[i, 2], pf[i, 3], pf[i, 4], pf[i, 5], scaling = scaling)
    SpxMs[, i] <- px(ages * scaling, pm[i, 1], pm[i, 2], pm[i, 3], pm[i, 4], pm[i, 5], scaling = scaling)
  }
  # arrange as data frame
  # Ttru
  pxFs <- data.frame(pxFs)
  pxMs <- data.frame(pxMs)
  # make single object
  # note dim(pxMF)
  # [1]  122 4000
  pxMF <- rbind(pxFs, pxMs)
  # other species
  SpxFs <- data.frame(SpxFs) 
  SpxMs <- data.frame(SpxMs)
  #--------------------------------------------------------------------------

  # get survival reduction - note if BND this is changed after
  srsims <- with(SpInfo, rbeta(nsims, sra, srb) * (sru - srl) + srl)
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "SR"){
      # if doing sensitivity analysis just use the mean
      meanSR <- with(SpInfo, (sra/(sra+srb)) * (sru - srl) + srl)
      srsims <- rep(meanSR,nsims)
    }
  }

  if (Sp == "LowSal") {
    # Read in the survival changes from NOAA
    LowSalS <- read.csv("InOutBySp/LowSal/BB_AllRegions_Survival.csv", header = TRUE)
    # Survival reduction per stratum is survival under APA (preferred alternative)
    #  divided by survival under NAA (no action alternative)
    #Note - assumes stratum is set - as has to be for LowSal analysis
    LowSalSR <- LowSalS[LowSalS$scenario == "APA", strata[stratum]] / 
      LowSalS[LowSalS$scenario == "NAA", strata[stratum]]
    #Sample nsims scenarios for each year, and turn into a matrix nsims * nyears
    #But first save the random number seed, and restore afterwards - to
    # make realizations run with LowSal compatible with those
    # for Ttru
    if (exists(".Random.seed", .GlobalEnv)) oldseed <- .GlobalEnv$.Random.seed
    else oldseed <- NULL
    LowSalSR <- matrix(LowSalSR[sample.int(length(LowSalSR), size = nsims * nyears, replace = TRUE)], nsims, nyears)
    if (!is.null(oldseed)) .GlobalEnv$.Random.seed <- oldseed
    else rm(".Random.seed", envir = .GlobalEnv)
  }

  #--------------------------------------------------------------------------
  # Initial population size and proportion exposed
  #  also accounting for stratum if required (for BB bottlenose dolphin)
  #--------------------------------------------------------------------------
  # These files were provided by LT and the originals are hosted at
  # C:\Users\tam2\Dropbox\Trabalho\Funded\CARMMHA\integrativemodelling\OffshoreAbundance2
  #--------------------------------------------------------------------------
  # Initial population size
  Nstart <- read.csv(paste0("InOutBySp/", SpInfo$folder, "/N_boot.csv"), header = TRUE)
  Nstart <- Nstart[, 2]
  Nstart <- as.vector(Nstart)
  # Population exposed
  Nexp <- read.csv(paste0("InOutBySp/", SpInfo$folder, "/N_boot_in_oil.csv"), header = TRUE)
  Nexp <- Nexp[, 2]
  Nexp <- as.vector(Nexp)

  #If a stratum analysis, get the proportion of the population in that stratum
  # and multiply Nstart and Nexp appropriately
  if(!is.null(stratum)){
    Dstr <- read.csv(paste0("InOutBySp/", SpInfo$folder, "/D_str.csv"), header = TRUE)
    pstr <- Dstr[, strata[stratum]]
    if(length(pstr) != length(Nstart)) 
      stop("D_str.csv and N_boot.csv have differing number of rows")
    Nstart <- Nstart * pstr
    Nexp <- Nexp * pstr
  }

  # removing extreme outliers that arise on occasion from GAM simulations
  mult.sd <- 5
  lim.outlier <- mean(Nstart) + mult.sd * sd(Nstart)
  index.NOT.outlier <- Nstart < lim.outlier
  # selecting just non outliers
  Nstart <- Nstart[index.NOT.outlier]
  Nexp <- Nexp[index.NOT.outlier]
  
  # Get proportion of population exposed
  pexp<-Nexp/Nstart
  # get the mean values for 
  # Nstart
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "N0"){
      # if doing sensitivity analysis just use the mean
      meanN0 <- mean(Nstart)
      Nstart <- rep(meanN0,nsims)
    }
  }
  # proportion exposed
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "pe"){
      # if doing sensitivity analysis just use the mean
      meanpexp<-mean(pexp)
      pexp <- rep(meanpexp,length(pexp))
    }
  }

  
  #--------------------------------------------------------------------------
  # The proportion of the population exposed that recovers
  pexprecsims <- with(SpInfo, rbeta(nsims, pra, prb) * (pru - prl) + prl)
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "per"){
      # if doing sensitivity analysis just use the mean
      meanper <- with(SpInfo, (pra/(pra+prb)) * (pru - prl) + prl)
      pexprecsims <- rep(meanper,nsims)
    }
  }
  pars2save$per <- pexprecsims
  #--------------------------------------------------------------------------

  #--------------------------------------------------------------------------
  # Fecundity
  # get DD fecundity parameters
  # Maximum fecundity
  Fmaxsims <- with(Ttru, rpert(nsims, min = minFmax, mode = modeFmax, max = maxFmax, shape = 4))
  # Scaling needs to happen here for maximum fecundity (i.e. minimum IBI) if not BND
  Fmaxsims <- Fmaxsims * scaling
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "Fmax"){
      # if doing sensitivity analysis just use the mean
      meanFmax <- with(Ttru,(minFmax + 4 * modeFmax + maxFmax)/(4 + 2))
      Fmaxsims <- rep(meanFmax,nsims)
    }
  }
  pars2save$Fmax <- Fmaxsims
  # F_nominal
  FnomsimsTtru <- with(Ttru, rpert(nsims, min = minFnom, mode = modeFnom, max = maxFnom, shape = 4))
  Fnomsims <- FnomsimsTtru * scaling
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "Fnom"){
      # if doing sensitivity analysis just use the mean
      meanFnom <- with(Ttru,(minFnom + 4 * modeFnom + maxFnom)/(4 + 2))
      Fnomsims <- rep(meanFnom,nsims) * scaling
      FnomsimsTtru <- rep(meanFnom,nsims)
    }
  }
  pars2save$Fnom <- Fnomsims
  # rho
  rhosims <- with(SpInfo, rhoshift + rgamma(nsims, shape = rhoshape, scale = rhoscale))
  
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "rho"){
      # if doing sensitivity analysis just use the mean
      meanrho <- with(SpInfo, rhoshift + rhoshape*rhoscale)
      rhosims <- rep(meanrho,nsims)
    }
  }
  pars2save$rho <- rhosims

  #--------------------------------------------------------------------------
  # age at first reproduction
  # this is for BB BND, then needs to be scaled for other taxa
  # if a regular Sim or a1r (age at first reproduction) is the parameter 
  # being evaluated in a sensitivity analysis
  a1stRsimsTtru <- with(Ttru, rgamma(nsims, shape = sha1str, scale = sca1str))
  a1stRsims <- a1stRsimsTtru
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "a1r"){
      # if doing sensitivity analysis for other parameter not a1r just repeat the mean a1r nsims times
      meana1R <- with(Ttru, sha1str * sca1str)
      a1stRsims <- rep(meana1R, nsims)
    }
  }
  # the modeling requires an integer and to adjust for possible scaling for non Ttru
  # the -0.5 is justified by a need to adjust the age of maturity to IBI
  # see SI for details
  a1stRsims <- round((a1stRsims - 0.5) / scaling)
  pars2save$a1r <- a1stRsims

  #------------------------------------------------------------------------------
  # Baseline reproductive success rate
  # get parameters of binomial
  beta.pars <- getBetaDistPars(SpInfo$meanbrs, SpInfo$sdbrs^2)
  # get Baseline reproductive success rate
  pRepbasesims <- rbeta(nsims, beta.pars$alpha, beta.pars$beta)
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "br"){
      # if doing sensitivity analysis just use the mean
      meanpRepbase <- beta.pars$alpha/(beta.pars$alpha+beta.pars$beta)
      pRepbasesims <- rep(meanpRepbase,nsims)
    }
  }
  pars2save$br <- pRepbasesims
  #------------------------------------------------------------------------------
  # Reproduction reduction
  # Post spill reproductive success rate
  # get parameters of binomial
  beta.pars <- getBetaDistPars(SpInfo$meanpors, SpInfo$sdpors^2)
  # get post spill reproductive success rate
  pRepPostsims <- rbeta(nsims, beta.pars$alpha, beta.pars$beta)
  
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "por"){
      # if doing sensitivity analysis just use the mean
      meanpor <- beta.pars$alpha/(beta.pars$alpha+beta.pars$beta)
      pRepPostsims  <- rep(meanpor,nsims)
    }
  }
  pars2save$por <- pRepPostsims
  
  #--------------------------------------------------------------------------
  # For BB BND
  # fecundity reduction for Ttru
  beta.pars <- getBetaDistPars(Ttru$meanbrs, Ttru$sdbrs^2)
  # get Baseline reproductive success rate for Ttru
  TtrupRepbasesims <- rbeta(nsims, beta.pars$alpha, beta.pars$beta)
  
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "BrTt"){
      # if doing sensitivity analysis just use the mean
      meanTtrupRepbase <- beta.pars$alpha/(beta.pars$alpha+beta.pars$beta)
      TtrupRepbasesims <- rep(meanTtrupRepbase,nsims)
    }
  }
  pars2save$BrTt <- TtrupRepbasesims
  
  # get post spill reproductive success rate for Ttru
  beta.pars <- getBetaDistPars(Ttru$meanpors, Ttru$sdpors^2)
  TtrupRepPostsims <- rbeta(nsims, beta.pars$alpha, beta.pars$beta)
  if (type != "Sim") {
    # if running a sensitivity analysis
    if(parS != "PorTt"){
      # if doing sensitivity analysis just use the mean
      meanTtrupRepPost <- beta.pars$alpha/(beta.pars$alpha+beta.pars$beta)
      TtrupRepPostsims <- rep(meanTtrupRepPost,nsims)
    }
  }
  pars2save$PorTt <- TtrupRepPostsims
  
  #--------------------------------------------------------------------------
  # Define a suitable object to hold the simulation data (i.e. the class sizes)
  #--------------------------------------------------------------------------
  # The array has 4 dimensions:
  #
  # dimension 1 - sex (female and male) and exposure (exposed and recovering, exposed and not recovering and non exposed)
  #               1 to 6*na:  1 is female exposed, recovering
  #                           2 is male   exposed, recovering
  #                           3 is female exposed, not recovering
  #                           4 is male   exposed, not recovering
  #                           5 is female non exposed, and
  #                           6 is male   non exposed
  #
  # dimension 2 - year after spill
  #               1 to nyears, year 1 is just after oil spill (so oil spill would have been year 0)
  #
  # dimension 3 - simulation number
  #               1 to number of sims = nsims
  #
  # dimension 4 - impact vs non impact
  #               1 is oil spill occurred, 2 is no oil spill occurred
  #--------------------------------------------------------------------------
  # simres[class, year, sim, reality] - names of the dimension
  # simres[dimm, nyears, nsims, 2] - size of the dimension
  simres <- array(data = NA, dim = c(dimm, nyears, nsims, 2), 
    dimnames = list(class = paste0("c", 1:dimm), 
                    year = paste0("y", 1:nyears), 
                    sim = paste0("s", 1:nsims), 
                    reality = c("spill", "nospill")))

  #--------------------------------------------------------------------------
  # Simulation implementation
  #--------------------------------------------------------------------------
  if (verbose) {
    cat("Species ",Sp," [", nsims, " simulations under mode ", type, 
        ifelse(is.null(parS), "", paste0(", parameter ", parS)), "]\n", sep = "")
    progbar <- progress::progress_bar$new(format = " Computing [:bar] :percent eta: :eta ",
            total = nsims)
  }


  # for each simulation iteration
  for (i in 1:nsims) {
    #Increment the progress bar
    if(verbose) progbar$tick()
    
    # Here go things that are constant across all years for a given iteration
    #------------------------------------------------------------
    # get initial conditions
    #------------------------------------------------------------

    # define iteration for the Siler model survival parameters
    iS <- itS[i]

    #------------------------------------------------------------
    # get initial population size
    # index to link pop size and survival post oil spill, if BND
    iterN <- sample.int(length(Nstart), 1)
    #--------------------------------------------------------------------------
    # Initial population size: N0
    # if a regular Sim or N0 is the parameter being evaluated in a sensitivity analysis
    N0sim <- Nstart[iterN]
    pars2save$N0[i] <- N0sim
    #------------------------------------------------------------------------------
    # proportion of animals in the population exposed to oil: pe
    # if a regular Sim or pe is the parameter being evaluated in a sensitivity analysis
    # First get the number of exposed animals
    pexpsim <- pexp[iterN]
    pars2save$pe[i] <- pexpsim

    #------------------------------------------------------------------------------
    # proportion of the population exposed that recovers: per
    # if a regular Sim or per is the parameter being evaluated in a sensitivity analysis
    pexprecsim <- pexprecsims[i]
    
    #------------------------------------------------------------------------------
    # Survival
    # if a regular Sim or ascS is the parameter being evaluated in a sensitivity analysis
    femalesur <- SpxFs[, iS]
    malesur <- SpxMs[, iS]
    if (type != "Sim") {
      # if running a sensitivity analysis
      if(parS != "ascS"){
        # if doing sensitivity analysis just use the mean
        femalesur <- rowMeans(SpxFs)
        malesur <- rowMeans(SpxMs)
      }
    }
    
    if (type != "Sim") {
      # if running a sensitivity analysis
      if(parS == "ascSE"){
        #only 3 sims are being run by definition
        femalesur <- multsE[i]*rowMeans(SpxFs)
        malesur <- multsE[i]*rowMeans(SpxMs)
      }
    }
    
    
    #saves meanS, only useful to plot elasticity for ascS
    #mean of female survival
    pars2save$ascS[i] <- mean(femalesur)
    
    #------------------------------------------------------------------------------
    # Survival reduction - for the current Sp
    # note if BND srsim is changed later
    srsim <- srsims[i]
    
    #------------------------------------------------------------
    # Fecundity
    # get DD fecundity parameters
    # if a regular Sim or Fmax is the parameter being evaluated in a sensitivity analysis
    Fmaxsim <- Fmaxsims[i]
  
    #------------------------------------------------------------
    # Fecundity
    # get DD fecundity parameters
    # if a regular Sim or Fnom is the parameter being evaluated in a sensitivity analysis
    Fnomsim <- Fnomsims[i]
    
    #------------------------------------------------------------
    # Fecundity
    # get DD fecundity parameters
    # if a regular Sim or rho is the parameter being evaluated in a sensitivity analysis
    rhosim <- rhosims[i]
    
    # obtain sims i age at first reproduction
    a1stRsim <- a1stRsims[i]
    
    # br Baseline reproductive success rate
    # if a regular Sim or br is the parameter being evaluated in a sensitivity analysis
    pRepbasesim <- pRepbasesims[i]

    #------------------------------------------------------------
    # Get transition matrix
    # note dimm = number of sexes * number of conditions * number of age classes
    # Nothe this object is only used way below but it's here to stay close to M0BND 
    M0 <- getM(na = na, dimm = ns * nc * na, femalesur = femalesur, malesur = malesur, 
               srf = 1, ddfr = Fnomsim, frf = 1, a1stR = a1stRsim, alasR = na)
    #------------------------------------------------------------------------------
    # Do the same for BND
    # note dimm = number of sexes * number of conditions * number of age classes
    M0BND <- getM(na = naBND, dimm = ns * nc * naBND, femalesur = pxFs[, iS], malesur = pxMs[, iS], 
                  srf = 1, ddfr = FnomsimsTtru[i], frf = 1, a1stR = a1stRsimsTtru[i], alasR = naBND)
    #------------------------------------------------------------   
    # get nominal distribution in age and sex class
    ev0 <- eigen(M0BND)
    distNominal <- Re(ev0$vectors[, 1]) / sum(Re(ev0$vectors[, 1]))
    
    # get population averaged survival pre-oil spill
    # see file "survivalReduction" for details
    # note we use the index iS even though there's no link between the pmarked and the siler model, 
    # but there is both 10000 samples in predictions and pxMF
    # and so that saves having to set a new index
    qijs <- c(distNominal[1:61] * predictions[iS, ], distNominal[62:122] * predictions[iS, ]) / sum(c(distNominal[1:61] * predictions[iS, ], distNominal[62:122] * predictions[iS, ]))
    # population baseline, i.e. pre oil spill, average survival consistent with SCR study
    meanS <- sum(pxMF[, iS] * qijs)
    #saves mean P(marked|age), only useful to plot elasticity for PM
    pars2save$PM[i] <- mean(predictions[iS, ])
    
    if (type != "Sim") {
      # if running a sensitivity analysis
      if(parS != "PM"){
        # note dimm = number of sexes * number of conditions * number of age classes
        M0BND <- getM(na = naBND, dimm = ns * nc * naBND, femalesur = rowMeans(pxFs), malesur = rowMeans(pxMs), 
                      srf = 1, ddfr = FnomsimsTtru[i], frf = 1, a1stR = a1stRsimsTtru[i], alasR = naBND)
        #------------------------------------------------------------   
        # get nominal distribution in age and sex class
        ev0 <- eigen(M0BND)
        distNominal <- Re(ev0$vectors[, 1]) / sum(Re(ev0$vectors[, 1]))
        # if doing sensitivity analysis just use the mean P(marked|age)
        qijs <- c(distNominal[1:61] * mpredictions, distNominal[62:122] * mpredictions) / sum(c(distNominal[1:61] * mpredictions, distNominal[62:122] * mpredictions))
        # population baseline, i.e. pre oil spill, average survival consistent with SCR study
        meanS <- sum(rowMeans(pxMF) * qijs)
        pars2save$PM[i] <- mean(mpredictions)
      }
    }
    
    if (type != "Sim") {
      # if running a sensitivity analysis
      if(parS == "PME"){
        #only 3 sims are being run by definition
        mpredictionsPME <- mpredictions * multsE[i]
        qijs <- c(distNominal[1:61] * mpredictions, distNominal[62:122] * mpredictionsPME) / sum(c(distNominal[1:61] * mpredictionsPME, distNominal[62:122] * mpredictionsPME))
        # population baseline, i.e. pre oil spill, average survival consistent with SCR study
        meanS <- sum(rowMeans(pxMF) * qijs)
        pars2save$PM[i] <- mean(mpredictionsPME)
      }
    }
    
    # spos survival post oil-spill for BB BND with Glennie et al. SCR analysis
    Spostspill <- SpostOilBB[iterN]
    pars2save$spos[i] <- Spostspill

    # Ttru survival reduction factor - to be used to scale fecundity reduction later
    srsimTtru <- Spostspill / meanS
    pars2save$srTtru[i] <- srsimTtru
    # constrain survival reduction - this could be changed to sample values of Spostspill that are lower than meanS
    if (srsimTtru > 1) srsimTtru <- 1
    
    # if we are in a Ttru sim, need to replace survival reduction with Ttru value
    if (Sp %in% BB_BND_Sp) {
      srsim <- srsimTtru
    }
    pars2save$SR[i] <- srsim

    # make population stable if not Ttru in Barataria Bay
    if (!(Sp %in% BB_BND_Sp)) {
      # this procedure is inspired by the algorithm described in an email
      # sent by LT to TAM on the Tue 6/2/2020 1:21 AM
      # solve for the fecundity that would lead to a constant population
      Fnomsim <- uniroot(getStableF0, interval = c(0.01, 1), ns1 = ns, nc1 = nc, na1 = na, femalesur1 = femalesur, 
                    malesur1 = malesur, a1stRsim1 = a1stRsim)$root
      # get transition matrix
      M0 <- getM(na = na, dimm = ns * nc * na, femalesur = femalesur, malesur = malesur, 
                 srf = 1, ddfr = Fnomsim, frf = 1, a1stR = a1stRsim, alasR = na)
    }

    #------------------------------------------------------------
    # Accounting for the impact of oil
    #------------------------------------------------------------

    #------------------------------------------------------------
    # number of years to return to baseline
    # (applies to both survival and fecundity)
    # _____________________________
    # Both sexes
    sa <- eigen.analysis(M0[1:(2 * na), 1:(2 * na)])$stable.stage
    # FEMALES
    # index
    FEMi <- 1:na
    # proportion
    propFEMALE <- sum(sa[FEMi])
    fmFEM <- fundamental.matrix(M0[FEMi, FEMi])
    # mean age
    maFEM <- sum(sa[FEMi] / sum(sa[FEMi]) * 0:(na - 1))
    # mean time to live
    mt2lFEM <- sum(sa[FEMi] / sum(sa[FEMi]) * fmFEM$meaneta)
    # MALES
    # index
    MALi <- (na + 1):(2 * na)
    # proportion
    propMALE <- sum(sa[MALi])
    fmMAL <- fundamental.matrix(M0[MALi, MALi])
    # mean age
    maMAL <- sum(sa[MALi] / sum(sa[MALi]) * 0:(na - 1))
    #  mean time to live
    mt2lMAL <- sum(sa[MALi] / sum(sa[MALi]) * fmMAL$meaneta)
    # average taking account of males and females
    mt2l <- mt2lFEM * propFEMALE + mt2lMAL * propMALE
    # _____________________________
    # years that an animal that recovers takes to recover to baseline survival levels
    yrbssim <- round(mt2l)
    pars2save$y2R[i] <- yrbssim
    #------------------------------------------------------------

    # yearly reduction factor for survival
    rfssim <- getRedFac(y1 = srsim, ny2r2n = yrbssim, ny = nyears)

    # por: post spill reproductive success rate
    # if a regular Sim or por is the parameter being evaluated in a sensitivity analysis
    pRepPostsim <- pRepPostsims[i]
    
    # reproduction reduction factor
    # if BB Ttru, no scaling required
    fecRed <- 1 - pRepPostsim / pRepbasesim
    # but if other species/stock not BND, some scaling is required
    # by multiplying the species survival reduction
    # by the ratio of the fecundity to survival reductions in Ttru
    # not to be confused with the scaling of survival for different species
    if (!(Sp %in% BB_BND_Sp)) {
      # BrTt Baseline reproductive success rate for Ttru
      TtrupRepbasesim <- TtrupRepbasesims[i]
      # PorTt post spill reproductive success rate for Ttru
      TtrupRepPostsim <- TtrupRepPostsims[i]
      #fecundity reduction      
      fecRedTtru <- 1 - TtrupRepPostsim / TtrupRepbasesim
      fecRed <- srsim * fecRedTtru / srsimTtru
    }
    pars2save$fecRed[i] <- fecRed

    # yearly reduction factor for fecundity
    rrfsim <- getRedFac(y1 = fecRed, ny2r2n = yrbssim, ny = nyears)
    
    #-----------------------------------------------------------
    # Initial population numbers per class
    # Nominal age distribution is normalized eigenvector associated with dominant eigenvalue
    ev0 <- eigen(M0)
    distNominal <- Re(ev0$vectors[, 1]) / sum(Re(ev0$vectors[, 1]))
    # initial proportions per age class
    inipopsim <- N0sim * distNominal
    # oil spill happened - classes where some (or all, e.g. BB) are exposed
    simres[, 1, i, 1] <- c(inipopsim[1:(2 * na)] * (1 - pexpsim), inipopsim[1:(2 * na)] * pexpsim * (1 - pexprecsim), inipopsim[1:(2 * na)] * pexpsim * pexprecsim)
    # oil spill did not happen - classes where everyone is not exposed
    simres[, 1, i, 2] <- inipopsim

    # Here go things that change each year ------------------
    for (j in 1:(nyears - 1)) {

      # first, get the current N needed for the density dependent fecundity calculations
      currNoil <- sum(simres[, j, i, 1]) # classes from oil spill scenario
      currNnooil <- sum(simres[, j, i, 2]) # classes from no-oil spill scenario

      # get DD fecundity rate -------------------------------
      # DD in an oil spill scenario
      ftsimoil <- ft(Nt = currNoil, Fmax = Fmaxsim, rho = rhosim, Nnom = N0sim, Fnom = Fnomsim)
      
      # get appropriate low salinity survival multiplier
      LowSalSR_Yr <- 1
      if (Sp == "LowSal") {
        if(j >= (2027 - 2010 + 1)) LowSalSR_Yr <- LowSalSR[i, j]
      }
      # Get the actual transition matrix --------------------
      # oil spill scenario
      Moil <- getM(na = na, dimm = ns * nc * na, femalesur = femalesur * LowSalSR_Yr, 
                   malesur = malesur * LowSalSR_Yr, srf = rfssim[j], srf.yr1 = rfssim[1],
                   ddfr = ftsimoil, frf = rrfsim[j], frf.yr1 = rrfsim[1], 
                   a1stR = a1stRsim, alasR = na)
      # If we are under constant pop, i.e. under any stock not BND Barataria Bay
      Mnooil <- M0
      # no-oil spill scenario (survival and fecundity reduction factors = 1)
      if (Sp %in% BB_BND_Sp) {
        # since we are under no constant pop, we need to do this step
        # DD for no-oil spill scenario
        ftsimnooil <- ft(Nt = currNnooil, Fmax = Fmaxsim, rho = rhosim, Nnom = N0sim, Fnom = Fnomsim)
        Mnooil <- getM(na = na, dimm = ns * nc * na, femalesur = femalesur, malesur = malesur, srf = 1, 
                       ddfr = ftsimnooil, frf = 1, a1stR = a1stRsim, alasR = na)
      }
      # Propagate population forward ------------------------
      # oil spill scenario at j+1 years
      simres[, j + 1, i, 1] <- Moil %*% simres[, j, i, 1]
      # no-oil spill scenario at j+1 years
      simres[, j + 1, i, 2] <- Mnooil %*% simres[, j, i, 2]
    }
    gc()
  }
  if(verbose) progbar$terminate()

  # time the procedure
  end_time <- Sys.time()
  TimeSpent <- end_time - start_time
  RunOnThe <- Sys.time()

  # calculate and plot injury measures
  injury <- getInjury(simres)

  #------------------------------------------------------------------
  # save results to corresponding species folder
  # file path
  if(Sp == "LowSal") {
    fp <- paste0("InOutBySp/LowSal/", Sp)
  } else {
    fp <- paste0("InOutBySp/", SpInfo$folder, "/", Sp)
  }
  # file type
  ft <- "simres"
  # file specifics - number of sims, Sim or Sens depending on type and parameter
  # note if parS is NULL i.e. when type=="Sim" then name is just as paste0(nsims,type)
  fs <- paste0(nsims, type, parS)
  # file extension
  fe <- ".RData"
  save(simres, pars2save, injury, RunOnThe, TimeSpent, file = getfilename(fp, ft, fs, fe))
  #Return nothing to calling function
  invisible()
}
