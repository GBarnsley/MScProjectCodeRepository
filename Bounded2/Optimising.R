#Please see naive or Covid version for more detailed annotations
load("SimulatedEpidemics.Rda")
library(EpiDataAugV2)
library(mcmcse)
set.seed(9301)
MaxRUNS <- 11^2 + 2*7^2 + 2*5^2 + 2*15
#just for keeping track
nimbleOptions(
    MCMCprogressBar = FALSE
)
###Setting epidemic length to be equal, so that we can use the same SIR Model on both
maxLength = max(length(SimulatedEpidemics[[1]]$newI), length(SimulatedEpidemics[[2]]$newI))
for(i in 1:length(SimulatedEpidemics)){
  newI <- rep(0, maxLength)
  newR <- rep(0, maxLength)
  for(j in 1:length(SimulatedEpidemics[[i]]$newI)){
    newI[j] <- SimulatedEpidemics[[i]]$newI[j]
    newR[j] <- SimulatedEpidemics[[i]]$newR[j]
  }
  SimulatedEpidemics[[i]]$newI <- newI
  SimulatedEpidemics[[i]]$newR <- newR
}
###Defining Samplers
sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = stepSampler_setup,
  run = stepSampler_run_bounded,
  methods = list(
    reset = function() {}
  )
)
###Constant Hyper Parameters
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
samples <- 1000
burnin <- 10000
thin <- 1
samplesPer <- 10
##Optimising TMax and DeltaMax
TMaxs <- seq(1,101, by = 10)
DeltaMaxs <- seq(1,101, by = 10)
#Setting up results data frame
Results <- data.frame(TMax = rep(0, length(TMaxs)*length(DeltaMaxs)), 
                      DeltaMax = rep(0, length(TMaxs)*length(DeltaMaxs)),
                      AcceptanceEpi1 = rep(0, length(TMaxs)*length(DeltaMaxs)),
                      MSJDEpi1 = rep(0, length(TMaxs)*length(DeltaMaxs)),
                      AcceptanceEpi2 = rep(0, length(TMaxs)*length(DeltaMaxs)),
                      MSJDEpi2 = rep(0, length(TMaxs)*length(DeltaMaxs))
)
ResultsRow <- 1
#setting up model
model <- SIR(newR = SimulatedEpidemics[[1]]$newR,
             N = SimulatedEpidemics[[1]]$Pop,
             Frequency = SimulatedEpidemics[[1]]$Frequency)
#finding initial values
initialNewIs <- list()
for(i in 1:length(SimulatedEpidemics)){
  model@Model$newR <- SimulatedEpidemics[[i]]$newR
  model <- initialValues(model, hyperParameters)
  initialNewIs[[i]] <- model@Model$newI
}
#setting up MCMC Config
mainMCMC <- configureMCMC(model@Model, nodes = NULL)
mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
mainMCMC$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
#Generating values
for(TMax in TMaxs){
  for(DeltaMax in DeltaMaxs){
    #building MCMC
    mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
    mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
    MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
    Results[ResultsRow, 1:2] <- c(TMax, DeltaMax)
    for(i in 1:length(SimulatedEpidemics)){
      model@Model$newR <- SimulatedEpidemics[[i]]$newR
      for(k in 1:samplesPer){
        #setting starting values values
        model@Model$Beta <- hyperParameters$`Initial Values`$Beta
        model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
        model@Model$newI <- initialNewIs[[i]]
        model@Model$calculate()
        #running MCMC
        MCMC$run(niter = samples*thin + burnin, nburnin = burnin, thin = thin)
        Samples <- as.matrix(MCMC$mvSamples)
        #storing
        Results[ResultsRow,(3:4) + 2*(i - 1)] <- Results[ResultsRow,(3:4) + 2*(i - 1)] + 
          c(AcceptanceProbability(Samples[,-c(1,2)]), 
            MeanSquareJumpDistance(Samples[,-c(1,2)]))/samplesPer
      }
    }
    ResultsRow <- ResultsRow + 1
    print(MaxRUNS - ResultsRow)
  }
}
nimble:::clearCompiled(model@Model)
#this clears all compiled code so we don't hit dll limit
model <- SIR(newR = SimulatedEpidemics[[1]]$newR,
             N = SimulatedEpidemics[[1]]$Pop,
             Frequency = SimulatedEpidemics[[1]]$Frequency)
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
MaxEp1 <- Results[which.max(Results$MSJDEpi1),1:2]
MaxEp1 <- c(MaxEp1$TMax, MaxEp1$DeltaMax)
MaxEp2 <- Results[which.max(Results$MSJDEpi2),1:2]
MaxEp2 <- c(MaxEp2$TMax, MaxEp2$DeltaMax)
seqProducer <- function(centrevalue){
  return(c(rev(seq(centrevalue, max(1, centrevalue - 9), by = -3)),
           seq(centrevalue + 3, centrevalue + 9, by = 3)))
}
TMaxs <- seqProducer(MaxEp1[1])
DeltaMaxs <- seqProducer(MaxEp1[2])
#rerunning
for(TMax in TMaxs){
  for(DeltaMax in DeltaMaxs){
    if(!identical(c(TMax, DeltaMax), MaxEp1)){
      #building MCMC
      mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
      mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
      MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
      Results[ResultsRow,] <- c(TMax, DeltaMax, 0, 0, 0, 0)
      for(i in 1:length(SimulatedEpidemics)){
        model@Model$newR <- SimulatedEpidemics[[i]]$newR
        for(k in 1:samplesPer){
          #setting starting values values
          model@Model$Beta <- hyperParameters$`Initial Values`$Beta
          model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
          model@Model$newI <- initialNewIs[[i]]
          model@Model$calculate()
          #running MCMC
          MCMC$run(niter = samples*thin + burnin, nburnin = burnin, thin = thin)
          Samples <- as.matrix(MCMC$mvSamples)
          #storing
          Results[ResultsRow,(3:4) + 2*(i - 1)] <-Results[ResultsRow,(3:4) + 2*(i - 1)] +
            c(AcceptanceProbability(Samples[,-c(1,2)]), 
              MeanSquareJumpDistance(Samples[,-c(1,2)]))/samplesPer
        }
      }
      ResultsRow <- ResultsRow + 1
      print(MaxRUNS - ResultsRow)
    }
  }
}
#recompiles model code
if(!identical(MaxEp1, MaxEp2)){
  nimble:::clearCompiled(model@Model)
  #this clears all compiled code so we don't hit dll limit
  model <- SIR(newR = SimulatedEpidemics[[1]]$newR,
               N = SimulatedEpidemics[[1]]$Pop,
               Frequency = SimulatedEpidemics[[1]]$Frequency)
  mainMCMC <- configureMCMC(model@Model, nodes = NULL)
  mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                      type = sampler_RW_block,
                      control = hyperParameters[["RWM"]])
  mainMCMC$addSampler(target = "newI",
                      type = sampler,
                      control = hyperParameters[["N+-Delta"]])
  mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
  
  TMaxs <- seqProducer(MaxEp2[1])
  DeltaMaxs <- seqProducer(MaxEp2[2])
  #rerunning
  for(TMax in TMaxs){
    for(DeltaMax in DeltaMaxs){
      if(!identical(c(TMax, DeltaMax), MaxEp1)){
        #building MCMC
        mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
        mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
        MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
        Results[ResultsRow,] <- c(TMax, DeltaMax, 0, 0, 0, 0)
        for(i in 1:length(SimulatedEpidemics)){
          model@Model$newR <- SimulatedEpidemics[[i]]$newR
          for(k in 1:samplesPer){
            #setting starting values values
            model@Model$Beta <- hyperParameters$`Initial Values`$Beta
            model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
            model@Model$newI <- initialNewIs[[i]]
            model@Model$calculate()
            #running MCMC
            MCMC$run(niter = samples*thin + burnin, nburnin = burnin, thin = thin)
            Samples <- as.matrix(MCMC$mvSamples)
            #storing
            Results[ResultsRow,(3:4) + 2*(i - 1)] <-Results[ResultsRow,(3:4) + 2*(i - 1)] +
              c(AcceptanceProbability(Samples[,-c(1,2)]), 
                MeanSquareJumpDistance(Samples[,-c(1,2)]))/samplesPer
          }
        }
        ResultsRow <- ResultsRow + 1
        print(MaxRUNS - ResultsRow)
      }
    }
  }
}
#Focusing on Maxima
MaxEp1 <- Results[which.max(Results$MSJDEpi1),1:2]
MaxEp1 <- c(MaxEp1$TMax, MaxEp1$DeltaMax)
MaxEp2 <- Results[which.max(Results$MSJDEpi2),1:2]
MaxEp2 <- c(MaxEp2$TMax, MaxEp2$DeltaMax)
seqProducer <- function(centrevalue){
  return(seq(max(1, centrevalue - 2), centrevalue + 2))
}
TMaxs <- seqProducer(MaxEp1[1])
DeltaMaxs <- seqProducer(MaxEp1[2])
#rerunning
for(TMax in TMaxs){
  for(DeltaMax in DeltaMaxs){
    if(!identical(c(TMax, DeltaMax), MaxEp1)){
      #building MCMC
      mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
      mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
      MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
      Results[ResultsRow,] <- c(TMax, DeltaMax, 0, 0, 0, 0)
      for(i in 1:length(SimulatedEpidemics)){
        model@Model$newR <- SimulatedEpidemics[[i]]$newR
        for(k in 1:samplesPer){
          #setting starting values values
          model@Model$Beta <- hyperParameters$`Initial Values`$Beta
          model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
          model@Model$newI <- initialNewIs[[i]]
          model@Model$calculate()
          #running MCMC
          MCMC$run(niter = samples*thin + burnin, nburnin = burnin, thin = thin)
          Samples <- as.matrix(MCMC$mvSamples)
          #storing
          Results[ResultsRow,(3:4) + 2*(i - 1)] <-Results[ResultsRow,(3:4) + 2*(i - 1)] +
            c(AcceptanceProbability(Samples[,-c(1,2)]), 
              MeanSquareJumpDistance(Samples[,-c(1,2)]))/samplesPer
        }
      }
      ResultsRow <- ResultsRow + 1
      print(MaxRUNS - ResultsRow)
    }
  }
}
#recompiles model code
if(!identical(MaxEp1, MaxEp2)){
  nimble:::clearCompiled(model@Model)
  #this clears all compiled code so we don't hit dll limit
  model <- SIR(newR = SimulatedEpidemics[[1]]$newR,
               N = SimulatedEpidemics[[1]]$Pop,
               Frequency = SimulatedEpidemics[[1]]$Frequency)
  mainMCMC <- configureMCMC(model@Model, nodes = NULL)
  mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                      type = sampler_RW_block,
                      control = hyperParameters[["RWM"]])
  mainMCMC$addSampler(target = "newI",
                      type = sampler,
                      control = hyperParameters[["N+-Delta"]])
  mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
  
  TMaxs <- seqProducer(MaxEp2[1])
  DeltaMaxs <- seqProducer(MaxEp2[2])
  #rerunning
  for(TMax in TMaxs){
    for(DeltaMax in DeltaMaxs){
      if(!identical(c(TMax, DeltaMax), MaxEp1)){
        #building MCMC
        mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- TMax
        mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- DeltaMax
        MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
        Results[ResultsRow,] <- c(TMax, DeltaMax, 0, 0, 0, 0)
        for(i in 1:length(SimulatedEpidemics)){
          model@Model$newR <- SimulatedEpidemics[[i]]$newR
          for(k in 1:samplesPer){
            #setting starting values values
            model@Model$Beta <- hyperParameters$`Initial Values`$Beta
            model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
            model@Model$newI <- initialNewIs[[i]]
            model@Model$calculate()
            #running MCMC
            MCMC$run(niter = samples*thin + burnin, nburnin = burnin, thin = thin)
            Samples <- as.matrix(MCMC$mvSamples)
            #storing
            Results[ResultsRow,(3:4) + 2*(i - 1)] <-Results[ResultsRow,(3:4) + 2*(i - 1)] +
              c(AcceptanceProbability(Samples[,-c(1,2)]), 
                MeanSquareJumpDistance(Samples[,-c(1,2)]))/samplesPer
          }
        }
        ResultsRow <- ResultsRow + 1
        print(MaxRUNS - ResultsRow)
      }
    }
  }
}
nimble:::clearCompiled(model@Model)
#this clears all compiled code so we don't hit dll limit
model <- SIR(newR = SimulatedEpidemics[[1]]$newR,
             N = SimulatedEpidemics[[1]]$Pop,
             Frequency = SimulatedEpidemics[[1]]$Frequency)
mainMCMC <- configureMCMC(model@Model, nodes = NULL)
mainMCMC$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
mainMCMC$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
mainMCMC$addMonitors(c('Beta', 'Gamma', 'newI'))
#recompiles model code

##Optimising R
MaxEp1 <- Results[which.max(Results$MSJDEpi1),1:2]
MaxEp2 <- Results[which.max(Results$MSJDEpi2),1:2]
MaxEpi <- list(c(MaxEp1$TMax, MaxEp1$DeltaMax), c(MaxEp2$TMax, MaxEp2$DeltaMax))
Rs <- seq(1, 15, by = 1)
Results2 <- data.frame(R = Rs, Epi1 = rep(0, length(Rs)), Epi2 = rep(0, length(Rs)))
for(i in 1:length(MaxEpi)){
  mainMCMC$samplerConfs[[2]]@.xData$control$TMax <- MaxEpi[[i]][1]
  mainMCMC$samplerConfs[[2]]@.xData$control$DeltaMax <- MaxEpi[[i]][2]
  model@Model$newR <- SimulatedEpidemics[[i]]$newR
  for(j in 1:length(Rs)){
    mainMCMC$samplerConfs[[2]]@.xData$control$R <- Rs[j]
    MCMC <- compileNimble(buildMCMC(mainMCMC), project = model@Model, resetFunctions = TRUE)
    for(k in 1:samplesPer){
      model@Model$Beta <- hyperParameters$`Initial Values`$Beta
      model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
      model@Model$newI <- initialNewIs[[i]]
      model@Model$calculate()
      #running MCMC
      start <- Sys.time()
      MCMC$run(niter = samples*thin + burnin, nburnin = burnin, thin = thin)
      end <- Sys.time()
      #Calculating computationally adjusted CNeff
      Samples <- as.matrix(MCMC$mvSamples)
      CNeff <- multiESS(Samples[,c("Beta","Gamma")])/as.double(end-start, units = "secs")
      #storing
      Results2[j,i+1] <- Results2[j,i+1] + CNeff/samplesPer
    }
    print(30 - j - 15*(i-1))
  }
}
#Storing for analysis
Results <- list(TDeltaMax = Results, R = Results2)
save(Results, file = "Results.Rda")
