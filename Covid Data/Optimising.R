##Loading Data:
load("Covid.Rda")
##Loading Libraries:
library(EpiDataAug)
library(mcmcse)
#random seed
set.seed(93829)
#maximum runs for keep track
MaxRUNS <- 11^2 + 7^2 + 5^2 + 15
#hiding progress bars
nimbleOptions(
  MCMCprogressBar = FALSE
)
##Defining Sampler so that it can be called:
sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = stepSampler_setup,
  run = stepSampler_run,
  methods = list(
    reset = function() {}
  )
)
##Constant Hyper Parameters:
hyperParameters <- list(
  `Initial Values` = list(
    Beta = 1,
    Gamma = 0.01,
    Runs = 2000
  ),
  Priors = list(
    Beta = list( #Uninformative, low mean, high variance
      Shape = 0.0001,
      Rate = 0.0001
    ),
    Gamma = list(
      Shape = 0.0001*0.01,
      Rate = 0.0001
    )
  ),
  RWM = list(
    propCov = matrix(c(1,0,0,0.01), nrow = 2)
  ),
  `N+-Delta` = list(
    R = 1
  )
)
samples <- 5000
burnin <- 15000
thin <- 1
samplesPer <- 100
##Optimising TMax and DeltaMax:
#defining values to check:
TMaxs <- seq(1,101, by = 10)
DeltaMaxs <- seq(1,101, by = 10)
#Setting up results data frame
Results <- data.frame(TMax = rep(0, length(TMaxs)*length(DeltaMaxs)), 
                      DeltaMax = rep(0, length(TMaxs)*length(DeltaMaxs)),
                      MSJD = rep(0, length(TMaxs)*length(DeltaMaxs))
                      )
#row where will be store the next values generated
ResultsRow <- 1
#setting up model
model <- SIR(newR = CovidData$Simulated$newR,
             N = CovidData$Simulated$Pop,
             Frequency = CovidData$Simulated$Frequency)
#finding initial values
model <- initialValues(model, hyperParameters)
#storing initial value so we don't have to repeat
initialNewI <- model@Model$newI

#setting up MCMC Config, we can reuse this
mainMCMC <- configureMCMC(model@Model, nodes = NULL)
mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
mainMCMC$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
#Generating values
#Note: we don't use the methods defined in the library to save time recompile
#the base epidemic model
for(TMax in TMaxs){
  for(DeltaMax in DeltaMaxs){
    #building MCMC
    mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
    mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
    MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
    Results[ResultsRow, 1:2] <- c(TMax, DeltaMax)
    #burnining in
    model@Model$Beta <- hyperParameters$`Initial Values`$Beta
    model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    model@Model$newI <- initialNewI
    model@Model$calculate()
    MCMC$run(niter = burnin, nburnin = burnin)
    #running
    for(k in 1:samplesPer){
      MCMC$mvSamples$resize(rows  = 0)
      MCMC$run(niter = samples*thin, thin = thin, reset = FALSE)
      Samples <- as.matrix(MCMC$mvSamples)
      #storing
      Results[ResultsRow,3] <- Results[ResultsRow,3] + 
        MeanSquareJumpDistance(Samples[,-c(1,2)])/samplesPer
    }
    ResultsRow <- ResultsRow + 1
    print(MaxRUNS - ResultsRow)
  }
}
nimble:::clearCompiled(model@Model)
#this clears all compiled code so we don't hit dll limit
model <- SIR(newR = CovidData$Simulated$newR,
             N = CovidData$Simulated$Pop,
             Frequency = CovidData$Simulated$Frequency)
mainMCMC <- configureMCMC(model@Model, nodes = NULL)
mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
mainMCMC$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
#recompiles model code
#Focusing on Maxima
set.seed(882)
MaxEp <- Results[which.max(Results$MSJD),1:2]
MaxEp <- c(MaxEp$TMax, MaxEp$DeltaMax)
seqProducer <- function(centrevalue){
  return(c(rev(seq(centrevalue, max(1, centrevalue - 9), by = -3)),
           seq(centrevalue + 3, centrevalue + 9, by = 3)))
}
TMaxs <- seqProducer(MaxEp[1])
DeltaMaxs <- seqProducer(MaxEp[2])
#rerunning
for(TMax in TMaxs){
  for(DeltaMax in DeltaMaxs){
    if(!identical(c(TMax, DeltaMax), MaxEp)){
      #building MCMC
      mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
      mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
      MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
      Results[ResultsRow,] <- c(TMax, DeltaMax, 0)
      #burnining in
      model@Model$Beta <- hyperParameters$`Initial Values`$Beta
      model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
      model@Model$newI <- initialNewI
      model@Model$calculate()
      MCMC$run(niter = burnin, nburnin = burnin)
      #running
      for(k in 1:samplesPer){
        MCMC$mvSamples$resize(rows  = 0)
        MCMC$run(niter = samples*thin, thin = thin, reset = FALSE)
        Samples <- as.matrix(MCMC$mvSamples)
        #storing
        Results[ResultsRow,3] <- Results[ResultsRow,3] + 
          MeanSquareJumpDistance(Samples[,-c(1,2)])/samplesPer
      }
      ResultsRow <- ResultsRow + 1
      print(MaxRUNS - ResultsRow)
    }
  }
}
nimble:::clearCompiled(model@Model)
#this clears all compiled code so we don't hit dll limit
model <- SIR(newR = CovidData$Simulated$newR,
             N = CovidData$Simulated$Pop,
             Frequency = CovidData$Simulated$Frequency)
mainMCMC <- configureMCMC(model@Model, nodes = NULL)
mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
mainMCMC$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
#recompiles model code
#Focusing on Maxima
MaxEp <- Results[which.max(Results$MSJD),1:2]
MaxEp <- c(MaxEp$TMax, MaxEp$DeltaMax)
seqProducer <- function(centrevalue){
  return(seq(max(1, centrevalue - 2), centrevalue + 2))
}
TMaxs <- seqProducer(MaxEp[1])
DeltaMaxs <- seqProducer(MaxEp[2])
#rerunning
for(TMax in TMaxs){
  for(DeltaMax in DeltaMaxs){
    if(!identical(c(TMax, DeltaMax), MaxEp)){
      #building MCMC
      mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
      mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
      MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
      Results[ResultsRow, 1:3] <- c(TMax, DeltaMax, 0)
      #burnining in
      model@Model$Beta <- hyperParameters$`Initial Values`$Beta
      model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
      model@Model$newI <- initialNewI
      model@Model$calculate()
      MCMC$run(niter = burnin, nburnin = burnin)
      #running
      for(k in 1:samplesPer){
        MCMC$mvSamples$resize(rows  = 0)
        MCMC$run(niter = samples*thin, thin = thin, reset = FALSE)
        Samples <- as.matrix(MCMC$mvSamples)
        #storing
        Results[ResultsRow,3] <- Results[ResultsRow,3] + 
          MeanSquareJumpDistance(Samples[,-c(1,2)])/samplesPer
      }
      ResultsRow <- ResultsRow + 1
      print(MaxRUNS - ResultsRow)
    }
  }
}
nimble:::clearCompiled(model@Model)
#this clears all compiled code so we don't hit dll limit
model <- SIR(newR = CovidData$Simulated$newR,
             N = CovidData$Simulated$Pop,
             Frequency = CovidData$Simulated$Frequency)
mainMCMC <- configureMCMC(model@Model, nodes = NULL)
mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
mainMCMC$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
#recompiles model code
##Optimising R:
MaxEp <- Results[which.max(Results$MSJD),1:2]
MaxEp <- c(MaxEp$TMax, MaxEp$DeltaMax)
Rs <- seq(1, 15, by = 1)
Results2 <- data.frame(R = Rs, Neff = rep(0, length(Rs)), Error = rep(0, length(Rs)))
mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- MaxEp[1]
mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- MaxEp[2]
for(j in 1:length(Rs)){
  mainMCMC$samplerConfs[[2]]@.xData$control$R <- Rs[j]
  MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
  #burnining in
  model@Model$Beta <- hyperParameters$`Initial Values`$Beta
  model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  model@Model$newI <- initialNewI
  model@Model$calculate()
  MCMC$run(niter = burnin, nburnin = burnin)
  #running
  MCMC$run(niter = samples*thin*samplesPer, thin = thin)
  Neffs <- rep(NA, samplesPer)
  for(k in 1:samplesPer){
    #setting starting values values
    MCMC$mvSamples$resize(rows  = 0)
    MCMC$run(niter = samples*thin, thin = thin, reset = FALSE)
    Samples <- as.matrix(MCMC$mvSamples)
    #storing
    Neffs[k] <- multiESS(Samples[,c("Beta","Gamma")])
  }
  Results2[j,2:3] <- c(mean(Neffs),
                       (qt(0.975,df=samplesPer - 1)/sqrt(samplesPer)) * sd(Neffs))
  print(15 - j)
}
##Storing for later analysis:
Results <- list(TDeltaMax = Results, R = Results2)
save(Results, file = "OptimisationResults.Rda")