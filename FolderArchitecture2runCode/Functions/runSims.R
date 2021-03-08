runSims <- function(Sp,nsims,nyears){
  # This function runs a simulation, for 
  # a given species, a number of iterations and for nyears
  # Inputs
  #       Sp - the species code, see details below
  #       nsims - the number of iterations to consider
  #       nyears - the number of years to run each simulation for
  #--------------------------------------------------------------------------
  # Details on "Sp" to define species to work with
  #--------------------------------------------------------------------------
  # Sp must be a (almost always) 4 letter character acronym
  # e.g. Pmac for sperm whale
  # if needed check file for acronym list correspondence to species
  # (acronym is the first column on the file "SpeciesDefinitionFile4offshore.xlsx")
  # Sp <- "Ttru"
  # Sp <- "Pmac"
  # Sp <- "Kosp"
  # Sp <- "Scly"
  # Sp <- "Bbry"
  # Sp <- "Pcra"
  # Sp <- "Pele"
  # Sp <- "Ttro"
  # Sp <- "Gmac"
  # Sp <- "Fatt"
  # Sp <- "Gris"
  # Sp <- "Sbre"
  # Sp <- "Satt"
  # Sp <- "Ttrs"
  # Sp <- "Slon"
  # Sp <- "Scoe"
  # To compare as close as possible with Schwackeetal2017
  # This required a number of approximations
  # Sp <- "Schwackeetal2017"
  #--------------------------------------------------------------------------
  # Outputs
  #       There is no single object produced at the end, 
  #       but there are several by-products of running the function 
  #       An R workspace is saved in a species specific folder
  #       The folder name is defined in "SpeciesDefinitionFile4offshore.xlsx"
  #       The workspace containins the following objects
  #       simres: raw results, with population size by age, sex and class per year 
  #               for each iteration of the simulation - see below details for object simres
  #       injury: object created by preprocesing simres using function getInjury that
  #               includes LCY, MPD and YTR for each iteration
  #       RunOnThe: time stamp for when the simulations were run
  #       TimeSpent: how long the simulations took to run
  #       fecs: the initial fecundities used in each iteration 
  # Additionally, 3 sets of plots are produced and saved
  #         1. plot of population trajectories
  #         2. plot of injury metrics
  #         3. plot of fecundities at year 0
  #--------------------------------------------------------------------------
  #--------------------------------------------------------------------------
  #--------------------------------------------------------------------------
  # to time the procedure 
  start_time <- Sys.time()
  # set a random seed, just because
  set.seed(1234)
  #--------------------------------------------------------------------------
  # required libraries and functions
  #--------------------------------------------------------------------------
  # load required libraries
  library(mc2d)     # for beta-pert distribution
  library(popbio)   #for computing mean age and expected time to live
  #--------------------------------------------------------------------------
  # load up required functions
  source("Functions/reqfuns.R")         # most functions are here
  source("Functions/SilerFuns.R")       # Siler model functions are here
  #--------------------------------------------------------------------------

  #--------------------------------------------------------------------------
  # Get all the relevant species information required to run the simulation
  # taken from file "SpeciesDefinitionFile4offshore.xlsx"
  # This file is hardwired in the code to reside in the subfolder "InputFiles"
  SpInfo <- getSpData(Sp)
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # Set counters and constants
  #--------------------------------------------------------------------------
  i       <- 1             # counter for the simulation iteration
  j       <- 1             # counter for the age class
  na      <- 61            # number of age classes for BND
  ns      <- 2             # number of sexes
  nc      <- 3             # number of conditions (see details above)
  dimm    <- ns * nc * na  # number of classes (2 sexes * 3 conditions * number of age classes)
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # BND data is required even if this is not a BND simulation
  #--------------------------------------------------------------------------
  # Get BB BND post oil Survival and population size from Glennie et al 2021 SCR study
  #--------------------------------------------------------------------------
  # Survival 
  SpostOilBB <- read.csv(file="InOutBySp/Bottlenose_dolphin_BB/PostOilSurv.csv",header=FALSE,col.names = "S")
  SpostOilBB <- as.numeric(SpostOilBB$S)
  # Population size
  N0BB <-read.csv(file="InOutBySp/Bottlenose_dolphin_BB/N_boot.csv",header=TRUE)
  N0BB <- as.numeric(N0BB$x)
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  #Scaling factor with respect to BND
  #--------------------------------------------------------------------------
  # that happens for the gestation duration scaling but also
  # for values for which we only have for BND (like the Fmax distribution)
  Ttru <- getSpData("Ttru")
  # gestation duration for Ttru
  gdTtru <- Ttru$gd
  # scaling for current species with respect to BND, to scale survivals
  scaling<-gdTtru/SpInfo$gd
  # matrix to hold results
  # might need to be bigger or smaller depending on scaling
  naBND  <- na
  na     <- ceiling(na/scaling)     
  dimm   <- ns * nc * na
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # Survival
  #--------------------------------------------------------------------------
  # reading in Siler parameters, for males and females
  # These files were provided by LS via email on the Fri 5/22/2020 4:10 PM
  pf <- read.csv("InputFiles/var_female 15May2020.csv")
  pm <- read.csv("InputFiles/var_male 15May2020.csv")
  # remove useless ID first column 
  pf<-pf[, -1]
  pm<-pm[, -1]
  # rename columns
  names(pf) <- names(pm) <- c("a1", "a2", "a3", "b1", "b3", "rmean")
  # get iteration of the siler model for each sim - which is a posterior distrbution
  #number of Siler model draws
  nS <- nrow(pf)
  #iteration of the Siler model to use for each of the sims
  itS <- sample(1:nS,size = nsims,replace=TRUE)
  #--------------------------------------------------------------------------
  #Baseline survival
  #--------------------------------------------------------------------------
  # Obtaining P(marked|age) required for survival pre-oil spill
  # See SI file "SurvivalReduction" for a detailed description of this procedure
  # including the production of some of the objects used
  #--------------------------------------------------------------------------
  #get probability of marked as a function of age
  predictions<-read.csv(file="InputFiles/pmarked.csv")
  predictions<-as.matrix(predictions)
  #how many predictions
  pred2use<-sample(1:nrow(predictions),size = nsims,replace=TRUE)
  #get sex and age specific survivals
  #--------------------------------------------------------------------------
  ages <- 0:60
  #replaced nrow(pf) with nsims
  ncols <- nS
  nrows <- length(ages)
  #get a dataframe to use ggplot2
  #objects to hold female and male realizations
  pxFs <- matrix(NA,nrow=nrows,ncol=ncols)
  pxMs <- matrix(NA,nrow=nrows,ncol=ncols)
  #for each observation of the posterior
  for(i in itS){
    #get the funtion
    pxFs[,i] <- px(ages,pf[i,1],pf[i,2],pf[i,3],pf[i,4],pf[i,5])
    pxMs[,i] <- px(ages,pm[i,1],pm[i,2],pm[i,3],pm[i,4],pm[i,5])
  }
  #arrange as data frame
  pxFs2 <- data.frame(pxFs,ages=ages)
  pxMs2 <- data.frame(pxMs,ages=ages)
  # add means and relevant quantiles
  pxFs2$mean=rowMeans(pxFs)
  pxMs2$mean=rowMeans(pxMs)
  #add the sex
  pxFs2$sex="Female"
  pxMs2$sex="Male"
  #make single object
  #note dim(pxMF)
  #[1]  122 4001
  #4001=4000 iterations plus sex
  pxMF <- rbind(pxFs2,pxMs2)
  #--------------------------------------------------------------------------
  ## get proportion of animals in each sex and age class
  # not here yet!
  #get survival reduction - note if BND this is changed after
  srsims <- with(SpInfo, rbeta(nsims, sra, srb)*(sru-srl)+srl)
  
  #--------------------------------------------------------------------------
  #Initial population size and proportion exposed
  #--------------------------------------------------------------------------
  # These files were provided by LT and the originals are hosted at
  # C:\Users\tam2\Dropbox\Trabalho\Funded\CARMMHA\integrativemodelling\OffshoreAbundance2
  #--------------------------------------------------------------------------
  # Initial population size
  Nstart <- read.csv(paste0("InOutBySp/",SpInfo$folder,"/N_boot.csv"), header = TRUE)
  Nstart <- Nstart[,2]
  Nstart<-as.vector(Nstart)
  # Population exposed
  Nexp <- read.csv(paste0("InOutBySp/",SpInfo$folder,"/N_boot_in_oil.csv"), header = TRUE)
  Nexp <- Nexp[,2]
  Nexp<-as.vector(Nexp)
  # Note: naturally the proportion exposed is Nexp/Nstart
  # removing Rogue iterations given bad GAM extrapolations
  # see document "IniPopSizes.html" for rationale
  mult.sd <- 5
  lim.rogue <- mean(Nstart) + mult.sd * sd(Nstart)
  index.NOT.rogue <- Nstart < lim.rogue
  #selecting just non rogue values
  Nstart <- Nstart[index.NOT.rogue]
  Nexp <- Nexp[index.NOT.rogue]
  #--------------------------------------------------------------------------
  #The proportion of the population exposed that recovers
  pexprecsims  <- with(SpInfo, rbeta(nsims, pra, prb)*(pru-prl)+prl)
  #--------------------------------------------------------------------------
 
  #--------------------------------------------------------------------------
  # Fecundity
  # get DD fecundity parameters
  Fmaxsims <- with(Ttru, rpert(nsims, min = minFmax, mode = modeFmax, max = maxFmax, shape = 4))
  Fnomsims <- with(Ttru, rpert(nsims, min = minFnom, mode = modeFnom, max = maxFnom, shape = 4))
  rhosims  <- with(SpInfo, rhoshift+rgamma(nsims, shape = rhoshape, scale = rhoscale))
  # SCALING needs to happen here for Maximum fecundity (i.e. minimum IBI) if not BND
  Fmaxsims <- Fmaxsims*scaling
  #--------------------------------------------------------------------------
  # obtain random deviate of age at first reproduction
  # this is for Tursiops, then scaled if other species
  # sample from right distribution
  a1stRsims <- with(Ttru, rgamma(nsims, shape=sha1str, scale=sca1str))
  # the modeling requires an integer
  # and to adjust for possible scaling for non Ttru
  # the -0.5 is justified by a need to adjust the age of maturity to IBI - see emails
  # between LT, LS and TAM (thread ends Sun 6/21/2020 8:28 PM)
  a1stRsims <- round((a1stRsims-0.5)/scaling)
  # and age at last reproduction - currently no senescent assumed, so commented out
  # alasRsim <- Ttru$alastRep
  # alasRsim <- round(alasRsim/scaling)
  #------------------------------------------------------------------------------
  # Baseline reproductive success rate
  # get parameters of binomial
  p <- SpInfo$meanbrs
  var.p <- SpInfo$sdbrs
  # transform to corresponding beta parameters
  alpha <- ((1 - p) / var.p - 1 / p) * p ^ 2 
  beta <- alpha * (1 / p - 1) 
  # get Baseline reproductive success rate
  pRepbasesims <- rbeta(nsims, alpha, beta) 
  #-----------------------------------------------------------
  # Get reproduction reduction ----------------------------
  # Post spill reproductive success rate
  # get parameters of binomial
  p <- SpInfo$meanpors
  var.p <- SpInfo$sdpors
  # transform to corresponding beta parameters
  alpha <- ((1 - p) / var.p - 1 / p) * p ^ 2 
  beta <- alpha * (1 / p - 1) 
  # get post spill reproductive success rate
  pRepPostsims <- rbeta(nsims, alpha, beta)
  #--------------------------------------------------------------------------
  #For BB BND 
  # fecundity reduction for Ttru
  p <- Ttru$meanbrs
  var.p <- Ttru$sdbrs
  # transform to corresponding beta parameters
  alpha <- ((1 - p) / var.p - 1 / p) * p ^ 2 
  beta <- alpha * (1 / p - 1) 
  # get Baseline reproductive success rate for Ttru
  TtrupRepbasesims <- rbeta(nsims, alpha, beta) 
  # get parameters of binomial
  p <- Ttru$meanpors
  var.p <- Ttru$sdpors
  # transform to corresponding beta parameters
  alpha <- ((1 - p) / var.p - 1 / p) * p ^ 2 
  beta <- alpha * (1 / p - 1) 
  # get post spill reproductive success rate for Ttru
  TtrupRepPostsims <- rbeta(nsims, alpha, beta) 
  #--------------------------------------------------------------------------
  
  
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
  #               1 to 150, year 1 is just after oil spill (so oil spill would have been year 0)
  # 
  # dimension 3 - simulation number
  #               1 to number of sims = nsims
  # 
  # dimension 4 - impact vs non impact
  #               1 is oil spill occurred, 2 is no oil spill occurred
  #--------------------------------------------------------------------------
  # set.seed(12345)
  # simres[class, year, sim, reality] - names of the dimension
  # simres[dimm, nyears, nsims, 2] - size of the dimension  
  simres <- array(data=NA, dim=c(dimm, nyears, nsims, 2),dimnames=list(class=paste0("c",1:dimm),year=paste0("y",1:nyears),sim=paste0("s",1:nsims),reality=c("spill","nospill")))
  #--------------------------------------------------------------------------
  #also saving fecundities as a check these are sensible
  fecs <- numeric(nsims)
  
  #--------------------------------------------------------------------------
  # Simulation implementation
  #--------------------------------------------------------------------------
  # for each simulation iteration 
  # uncomment next line to run SIMS, else this is in "quick debug" mode!
  for(i in 1:nsims) {
    # Here go things that are constant across all years for a given iteration
    #------------------------------------------------------------
    #get initial conditions
    #------------------------------------------------------------
    
    #define iteration for the Siler parameters
    iS<-itS[i]
    
    #------------------------------------------------------------
    # get initial population size
    # Note: for smooth functioning "Nstart" and "SpostOilBB" must have same length 
    # (currently 1000 values in each vector)
    # index to link pop size and survival post oil spill, if BND
    iterSCR <- sample(1:length(Nstart),1)
    N0sim <- Nstart[iterSCR]
    # get the population exposed
    Nexpsim <- Nexp[iterSCR]
    # get the proportion of the population exposed
    pexpsim <-Nexpsim/N0sim
    #------------------------------------------------------------
    
    #------------------------------------------------------------
    # get the proportion of the population exposed that recovers
    pexprecsim  <- pexprecsims[i]
    #------------------------------------------------------------
    
    #------------------------------------------------------------
    # Survival  
    # if sampling instead of using each iteration of the Siler posterior in turn
    # niterSiler <- nrow(pm)
    # iS <- sample(1:niterSiler, 1)
    # to make these above work replace i by iS in the 3 lines below 
    ages <- 0:(na-1)
    femalesur <- px(ages*scaling, pf[iS, 1], pf[iS, 2], pf[iS, 3], pf[iS, 4], pf[iS, 5], scaling=scaling)
    malesur   <- px(ages*scaling, pm[iS, 1], pm[iS, 2], pm[iS, 3], pm[iS, 4], pm[iS, 5], scaling=scaling)
    #------------------------------------------------------------
    
    #Survival reduction - for Sp
    srsim <- srsims[i]
    
    #------------------------------------------------------------
    # Fecundity
    # get DD fecundity parameters
    Fmaxsim <- Fmaxsims[i]
    Fnomsim <- Fnomsims[i]
    rhosim  <- rhosims[i]
    # I need this to make sure don't hit problems with Fnom > Fmax which leads to NaNs!
    # but is this OK ? 
    if(Fmaxsim < Fnomsim) Fnomsim <- Fmaxsim*0.99
    # note below we use the same N0sim twice 
    # this is somewhat strange as it means that F0 becomes Fnomsim - but here anyway
    # Get baseline population fecundity - note this will only be used as such for BB BND
    # for all other taxonomic units we will change the IBI to reach a stable population
    F0 <- ft(Nt = N0sim, Fmax = Fmaxsim, rho  = rhosim, Nnom = N0sim, Fnom = Fnomsim)
    # obtain random deviate of age at first reproduction
    a1stRsim <- a1stRsims[i]
    #get Baseline reproductive success rate
    pRepbasesim <- pRepbasesims[i]
    #------------------------------------------------------------
    
    #------------------------------------------------------------
    # Get transition matrix
    # note dimm = number of sexes * number of conditions * number of age classes
    M0 <- getM(na        = na, dimm = ns*nc*na, femalesur = femalesur, malesur = malesur, srf = 1, N = N0sim, ddfr = F0, frf = 1, a1stR = a1stRsim, alasR     = na)
    #------------------------------------------------------------
    # Do the same for BND
    # note dimm = number of sexes * number of conditions * number of age classes
    M0BND <- getM(na        = naBND, dimm = ns*nc*naBND, femalesur = pxFs[,iS], malesur = pxMs[,iS], srf = 1, N = N0sim, ddfr = F0, frf = 1, a1stR = a1stRsim, alasR     = naBND)
    #------------------------------------------------------------
    #get nominal distribution in age and sex class
    ev0 <- eigen(M0BND)
    distNominal <- Re(ev0$vectors[, 1]) / sum(Re(ev0$vectors[, 1]))
    
    #get population averaged survival pre-oil spill
    # see file "survivalReduction" for details 
    # note we use the index iS even though there's no link between the pmarked and the siler model, but there's both 4000 iterations in each
    # and so that saves having to set a new index
    qijs<-c(distNominal[1:61]*predictions[iS,],distNominal[62:122]*predictions[iS,])/sum(c(distNominal[1:61]*predictions[iS,],distNominal[62:122]*predictions[iS,]))
    #population baseline, i.e. pre oil spill, average survival consistent with SCR study
    meanS<-sum(pxMF[,iS]*qijs)
    #-----------------------------------------------------------
    
    # survival post oil-spill for BB BND with Richard Glennie's SCR analysis
    Spostspill <- SpostOilBB[iterSCR]
    # survival reduction factor - to be used to scale fecundity reduction later
    srsimTtru <- Spostspill/meanS
    # constrain survival reduction - this could be changed to sample values of Spostspill that are lower than meanS
    if(srsimTtru>1) srsimTtru <- 1

    
    # make population stable if not Ttru in Barataria Bay
    if(Sp!="Ttru" & Sp!="Schwackeetal2017"){
      # this procedure is inspired by the algorithm described in an email
      # sent by LT to TAM on the Tue 6/2/2020 1:21 AM
      
      # first, solve for the fecundity that would lead to a constant population
      F0<-uniroot(getStableF0,interval=c(0.01, 1), ns1=ns, nc1=nc, na1=na, femalesur1=femalesur, malesur1=malesur, N0sim1=N0sim, a1stRsim1=a1stRsim)$root
      
      # beta4rho <- getbeta(Fmax = Fmaxsim,rho = rhosim,Fnom = F0,Nnom = N0sim)
      # then, using that Fecundity as the nominal fecundity, 
      # get an updated fecundity that given the existing Fmax and N0 would
      # lead to a stable population
      Fnomsim <- F0
      F0 <- ft(Nt = N0sim, Fmax = Fmaxsim, rho = rhosim, Nnom = N0sim, Fnom = Fnomsim)
      # in some cases the required F0 might be higher than Fmax
      # make then F0 to be Fmax
      # probably not the best solution? needs thinking!
      if(is.nan(F0)) F0 <- Fmaxsim
      #save initial fecundity
      fecs[i]<-F0
      M0 <- getM(na = na, dimm = ns*nc*na, femalesur = femalesur, malesur = malesur, srf = 1, N = N0sim, ddfr = F0, frf = 1, a1stR = a1stRsim, alasR = na)
    }
 
    #------------------------------------------------------------
    # Accounting for the impact of oil
    #------------------------------------------------------------
    
    #------------------------------------------------------------
    # number of years to return to baseline 
    # (applies to both survival and fecundity)
    #_____________________________
    #Both sexes
    sa <- eigen.analysis(M0[1:(2*na), 1:(2*na)])$stable.stage
    #FEMALES
    # index
    FEMi <- 1:na
    #proportion
    propFEMALE <- sum(sa[FEMi])
    fmFEM <- fundamental.matrix(M0[FEMi, FEMi])
    #mean age
    maFEM <- sum(sa[FEMi]/sum(sa[FEMi]) * 0:(na-1))
    # mean time to live 
    mt2lFEM <- sum(sa[FEMi]/sum(sa[FEMi]) * fmFEM$meaneta)
    # MALES
    # index
    MALi <- (na+1):(2*na)
    #proportion 
    propMALE <- sum(sa[MALi])
    fmMAL <- fundamental.matrix(M0[MALi, MALi])
    #mean age
    maMAL <- sum(sa[MALi]/sum(sa[MALi]) * 0:(na-1))
    #  mean time to live 
    mt2lMAL <- sum(sa[MALi]/sum(sa[MALi]) * fmMAL$meaneta)
    #average taking account of males and females
    mt2l <- mt2lFEM*propFEMALE+mt2lMAL*propMALE
    #_____________________________
    # this would be if hardwired
    # but this now depends on the mean time an animal has to live
    # yrbssim <- SpInfo$yrbs
    #_____________________________
    yrbssim <- round(mt2l)
    #------------------------------------------------------------
    
    # yearly reduction factor for survival
    rfssim <- getRedFac(y1 = srsim, ny2r2n = yrbssim, ny = nyears)  
    
    # get post spill reproductive success rate
    pRepPostsim <- pRepPostsims[i]
    # reproduction reduction factor
    # if BB Ttru, no scaling required
    fecRed <- 1-pRepPostsim/pRepbasesim
    # but if other species/stock not BND, some scaling is required
    # by multiplying the species survival reduction
    # by the ratio of the fecundity to survival reductions in Ttru
    # not to be confused with the scaling of survival for different species
    if (Sp!="Ttru" & Sp!="Schwackeetal2017"){
      # get Baseline reproductive success rate for Ttru
      TtrupRepbasesim <- TtrupRepbasesims[i]
      # get post spill reproductive success rate for Ttru
      TtrupRepPostsim <- TtrupRepPostsims[i]
      fecRedTtru <- 1-TtrupRepPostsim/TtrupRepbasesim
      fecRed <- srsim * fecRedTtru/srsimTtru
    }
    
    rrfsim <- getRedFac(y1 = fecRed, ny2r2n = yrbssim, ny = nyears)
    #-----------------------------------------------------------
    
    #-----------------------------------------------------------
    # Initial population numbers per class
    # Nominal age distribution is normalized eigen vector associated with dominant eigen value
    ev0 <- eigen(M0)
    distNominal <- Re(ev0$vectors[, 1]) / sum(Re(ev0$vectors[, 1]))
    # initial proportions per age class
    inipopsim <- N0sim * distNominal
    # oil spill happened - classes where some (or all, e.g. BB) are exposed
    simres[, 1, i, 1] <- c(inipopsim[1:(2*na)]*(1-pexpsim), inipopsim[1:(2*na)]*pexpsim*(1-pexprecsim), inipopsim[1:(2*na)]*pexpsim*pexprecsim)
    # oil spill did not happen - classes where everyone is not exposed
    simres[, 1, i, 2] <- inipopsim
    
    
    # Here go things that change each year ------------------
    for(j in 1:(nyears - 1)) {
      
      # first, get the current N for the density dependence
      currNoil   <- sum(simres[, j, i, 1]) # classes from oil spill scenario
      currNnooil <- sum(simres[, j, i, 2]) # classes from no-oil spill scenario
      
      # get DD fecundity rate -------------------------------
      # DD in an oil spill scenario
      ftsimoil   <- ft(Nt = currNoil, Fmax = Fmaxsim, rho = rhosim, Nnom = N0sim, Fnom = Fnomsim)
      # Get the actual transition matrix --------------------
      # oil spill scenario
      Moil   <- getM(na = na, dimm = ns* nc * na, femalesur = femalesur, malesur = malesur, srf = rfssim[j], N = currNoil, ddfr = ftsimoil, frf = rrfsim[j], a1stR = a1stRsim, alasR= na)
      # If we are under constant pop, i.e. under any stock not BND Barataria Bay
      Mnooil <- M0
      # no-oil spill scenario (survival reduction factor = 1)
      if(Sp=="Ttru" | Sp=="Schwackeetal2017"){
        # since we are under no constant pop, we need to do this step
        # DD for no-oil spill scenario
        ftsimnooil <-   ft(Nt = currNnooil, Fmax = Fmaxsim, rho = rhosim, Nnom = N0sim, Fnom = Fnomsim)
        Mnooil <- getM(na = na, dimm = ns*nc*na, femalesur = femalesur, malesur = malesur, srf = 1, N = currNnooil, ddfr = ftsimnooil, frf = 1, a1stR = a1stRsim, alasR= na)
      }
      # Propagate population forward ------------------------
      # oil spill scenario at j+1 years
      simres[, j + 1, i, 1] <- Moil %*% simres[, j, i, 1]
      # no-oil spill scenario at j+1 years
      simres[, j + 1, i, 2] <- Mnooil %*% simres[, j, i, 2]
    }
    if( i == 1) cat("Counting simulations (from a total of ", nsims,"):\n", sep="")
    if( i %% 10 == 0) cat(i, "\n", sep="")
    gc()
  }
  
  # save results to corresponding species folder
  save(simres, file=paste0("InOutBySp/", SpInfo$folder, "/", Sp, "simres", nsims, ".RData"))
  
  # time the procedure
  end_time <- Sys.time()
  TimeSpent <- end_time - start_time
  RunOnThe <- Sys.time()
  #------------------------------------------------------------------
  # post-processing
  #------------------------------------------------------------------
  
  # plot initial populations
  plotIP(Nstart, Nexp, type=1)
  plotIP(Nstart, Nexp, type=2)
  
  # plot simulation summary
  plotSims(simres)
  
  # calculate and plot the injury measures
  injury <- getInjury(simres)
  #------------------------------------------------------------------
  
  
  # SAVING RESULTS
  #--------------------------------------------------------------------------
  # make plot of initial population size and populatio/proportion expose
  png(filename = paste0("InOutBySp/", SpInfo$folder, "/", Sp, "TSvsEA", nsims, ".png"), width = 3*480,height = 2*480)
  # plot initial populations
  plotIP(Nstart, Nexp, type=1)
  dev.off()
  png(filename = paste0("InOutBySp/", SpInfo$folder, "/", Sp, "TSvsPEA", nsims, ".png"), width = 3*480,height = 2*480)
  # plot initial populations
  plotIP(Nstart, Nexp, type=2)
  dev.off()
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # make plot of population trajectories
  png(filename = paste0("InOutBySp/", SpInfo$folder, "/", Sp, "PopTrajs", nsims, ".png"), width = 3*480,height = 2*480)
  plotSims(simres)
  dev.off()
  #--------------------------------------------------------------------------
  #--------------------------------------------------------------------------
  # make plot of injury
  png(filename = paste0("InOutBySp/", SpInfo$folder, "/", Sp, "Injury", nsims, ".png"), width = 3*480,height = 2*480)
  injury <- getInjury(simres)
  dev.off()
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # make plot of fecundities at year 0
  png(filename = paste0("InOutBySp/", SpInfo$folder, "/", Sp, "Fecundity", nsims, ".png"), width = 3*480,height = 2*480)
  hist(fecs)
  dev.off()
  #--------------------------------------------------------------------------

  # save results to corresponding species folder
  save(simres, injury, RunOnThe, fecs, TimeSpent, file=paste0("InOutBySp/", SpInfo$folder, "/", Sp, "simres", nsims, ".RData"))  
    
}