##Loading Libraries
library(ggplot2)
library(mcmcse)
set.seed(193021)
##Loading Data:
setwd("../Naive")
load("SimulatedEpidemics.Rda")
#Setting up list of optimal parameters
parameters <- list()
parameters[[1]] <- list() #for each epidemic
parameters[[2]] <- list()
parameters[[1]][[1]] <- list(
  R = 1,
  TMax = 4,
  DeltaMax = 19
)
parameters[[2]][[1]] <- list(
  R = 1,
  TMax = 70,
  DeltaMax = 5
)
parameters[[1]][[2]] <- list(
  R = 1,
  TMax = 5,
  DeltaMax = 18,
  Bounded = TRUE
)
parameters[[2]][[2]] <- list(
  R = 1,
  TMax = 80,
  DeltaMax = 4,
  Bounded = TRUE
)
parameters[[1]][[3]] <- list(
  R = 6,
  TMax = 4,
  DeltaMax = 19,
  Bounded = TRUE
)
parameters[[2]][[3]] <- list(
  R = 12,
  TMax = 100,
  DeltaMax = 8,
  Bounded = TRUE
)
samples <- 1000
burnin <- 15000
thin <- 1
runs <- 250
TotalRuns <- runs*6
#hyper parameters
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
  )
)
##initialising so we don't have to repeat:
library(EpiDataAug)
nimbleOptions(
  MCMCprogressBar = FALSE
)
SIRS <- list()
initialNewIs <- list()
for(i in 1:length(SimulatedEpidemics)){
  SIRS[[i]] <- SIR(newR = SimulatedEpidemics[[i]]$newR,
                 N =SimulatedEpidemics[[i]]$Pop,
                 Frequency = SimulatedEpidemics[[i]]$Frequency)
  initialNewIs[[i]] <- initialValues(SIRS[[i]],
                                     hyperParameters)@Model$newI
}
##Results Data Frame
Results <- data.frame(Algorithm = rep(c("Naive", "Bounded 1", "Bounded 2"),2),
                      Epidemic = c(rep("Short", 3), rep("Long", 3)),
                      MSJD = rep(0, 6),
                      UpperBound = rep(0, 6),
                      LowerBound = rep(0, 6),
                      Neff = rep(0, 6),
                      UpperBound = rep(0, 6),
                      LowerBound = rep(0, 6),
                      CNeff = rep(0, 6),
                      UpperBound = rep(0, 6),
                      LowerBound = rep(0, 6)
                      )
##Naive Algorithm:
RowIndexes <- c(1,4) #positions to store results
ParameterIndex <- 1
for(i in 1:length(RowIndexes)){
  #for each epidemic
  hyperParameters$`N+-Delta` <- parameters[[i]][[ParameterIndex]]
  SIRS[[i]]@MCMC <- buildMCMCInternal(SIRS[[i]], hyperParameters)
  MSJDs <- rep(0, runs)
  NEFFs <- rep(0, runs)
  CNEFFs <- rep(0, runs)
  for(j in 1:runs){
    SIRS[[i]]@Model$newI <- initialNewIs[[i]]
    SIRS[[i]]@Model$Beta <- hyperParameters$`Initial Values`$Beta
    SIRS[[i]]@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    SIRS[[i]]@Model$calculate()
    start <- Sys.time()
    SIRS[[i]]@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    end <- Sys.time()
    Samples <- as.matrix(SIRS[[i]]@MCMC$mvSamples)
    MSJDs[j] <- MeanSquareJumpDistance(Samples[,-c(1,2)])
    NEFFs[j] <- multiESS(Samples[,c("Beta","Gamma")])
    CNEFFs[j] <- NEFFs[j]/as.double(end-start, units = "secs")
    TotalRuns <- TotalRuns - 1
    print(TotalRuns)
  }
  error <- qt(0.975,df=runs-1)/sqrt(runs)
  Results[RowIndexes[i],3:11] <- c(
    mean(MSJDs),
    mean(MSJDs) + error*sd(MSJDs),
    mean(MSJDs) - error*sd(MSJDs),
    mean(NEFFs),
    mean(NEFFs) + error*sd(NEFFs),
    mean(NEFFs) - error*sd(NEFFs),
    mean(CNEFFs),
    mean(CNEFFs) + error*sd(CNEFFs),
    mean(CNEFFs) - error*sd(CNEFFs)
  )
}
##Bounded Version 1:
library(EpiDataAugV1)
nimbleOptions(
  MCMCprogressBar = FALSE
)
RowIndexes <- RowIndexes + 1
ParameterIndex <- ParameterIndex + 1
for(i in 1:length(RowIndexes)){
  hyperParameters$`N+-Delta` <- parameters[[i]][[ParameterIndex]]
  SIRS[[i]]@MCMC <- buildMCMCInternal(SIRS[[i]], hyperParameters)
  MSJDs <- rep(0, runs)
  NEFFs <- rep(0, runs)
  CNEFFs <- rep(0, runs)
  for(j in 1:runs){
    SIRS[[i]]@Model$newI <- initialNewIs[[i]]
    SIRS[[i]]@Model$Beta <- hyperParameters$`Initial Values`$Beta
    SIRS[[i]]@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    SIRS[[i]]@Model$calculate()
    start <- Sys.time()
    SIRS[[i]]@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    end <- Sys.time()
    Samples <- as.matrix(SIRS[[i]]@MCMC$mvSamples)
    MSJDs[j] <- MeanSquareJumpDistance(Samples[,-c(1,2)])
    NEFFs[j] <- multiESS(Samples[,c("Beta","Gamma")])
    CNEFFs[j] <- NEFFs[j]/as.double(end-start, units = "secs")
    TotalRuns <- TotalRuns - 1
    print(TotalRuns)
  }
  error <- qt(0.975,df=runs-1)/sqrt(runs)
  Results[RowIndexes[i],3:11] <- c(
    mean(MSJDs),
    mean(MSJDs) + error*sd(MSJDs),
    mean(MSJDs) - error*sd(MSJDs),
    mean(NEFFs),
    mean(NEFFs) + error*sd(NEFFs),
    mean(NEFFs) - error*sd(NEFFs),
    mean(CNEFFs),
    mean(CNEFFs) + error*sd(CNEFFs),
    mean(CNEFFs) - error*sd(CNEFFs)
  )
}
##Bounded Version 2:

#warning: takes a long time since R = 12 for one of the set ups!!

library(EpiDataAugV2)
nimbleOptions(
  MCMCprogressBar = FALSE
)
RowIndexes <- RowIndexes + 1
ParameterIndex <- ParameterIndex + 1
for(i in 1:length(RowIndexes)){
  hyperParameters$`N+-Delta` <- parameters[[i]][[ParameterIndex]]
  SIRS[[i]]@MCMC <- buildMCMCInternal(SIRS[[i]], hyperParameters)
  MSJDs <- rep(0, runs)
  NEFFs <- rep(0, runs)
  CNEFFs <- rep(0, runs)
  for(j in 1:runs){
    SIRS[[i]]@Model$newI <- initialNewIs[[i]]
    SIRS[[i]]@Model$Beta <- hyperParameters$`Initial Values`$Beta
    SIRS[[i]]@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    SIRS[[i]]@Model$calculate()
    start <- Sys.time()
    SIRS[[i]]@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    end <- Sys.time()
    Samples <- as.matrix(SIRS[[i]]@MCMC$mvSamples)
    MSJDs[j] <- MeanSquareJumpDistance(Samples[,-c(1,2)])
    NEFFs[j] <- multiESS(Samples[,c("Beta","Gamma")])
    CNEFFs[j] <- NEFFs[j]/as.double(end-start, units = "secs")
    TotalRuns <- TotalRuns - 1
    print(TotalRuns)
  }
  error <- qt(0.975,df=runs-1)/sqrt(runs)
  Results[RowIndexes[i],3:11] <- c(
    mean(MSJDs),
    mean(MSJDs) + error*sd(MSJDs),
    mean(MSJDs) - error*sd(MSJDs),
    mean(NEFFs),
    mean(NEFFs) + error*sd(NEFFs),
    mean(NEFFs) - error*sd(NEFFs),
    mean(CNEFFs),
    mean(CNEFFs) + error*sd(CNEFFs),
    mean(CNEFFs) - error*sd(CNEFFs)
  )
}
##Bounded again but R=1
Results[7,1:2] <- c("Bounded2, R = 1", "Short") 
Results[8,1:2] <- c("Bounded2, R = 1", "Long") 
library(EpiDataAugV2)
nimbleOptions(
  MCMCprogressBar = FALSE
)
RowIndexes <- c(7,8)
for(i in 1:length(RowIndexes)){
  hyperParameters$`N+-Delta` <- parameters[[i]][[ParameterIndex]]
  hyperParameters$`N+-Delta`$R <- 1
  SIRS[[i]]@MCMC <- buildMCMCInternal(SIRS[[i]], hyperParameters)
  MSJDs <- rep(0, runs)
  NEFFs <- rep(0, runs)
  CNEFFs <- rep(0, runs)
  for(j in 1:runs){
    SIRS[[i]]@Model$newI <- initialNewIs[[i]]
    SIRS[[i]]@Model$Beta <- hyperParameters$`Initial Values`$Beta
    SIRS[[i]]@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    SIRS[[i]]@Model$calculate()
    start <- Sys.time()
    SIRS[[i]]@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    end <- Sys.time()
    Samples <- as.matrix(SIRS[[i]]@MCMC$mvSamples)
    MSJDs[j] <- MeanSquareJumpDistance(Samples[,-c(1,2)])
    NEFFs[j] <- multiESS(Samples[,c("Beta","Gamma")])
    CNEFFs[j] <- NEFFs[j]/as.double(end-start, units = "secs")
    TotalRuns <- TotalRuns - 1
    print(TotalRuns)
  }
  error <- qt(0.975,df=runs-1)/sqrt(runs)
  Results[RowIndexes[i],3:11] <- c(
    mean(MSJDs),
    mean(MSJDs) + error*sd(MSJDs),
    mean(MSJDs) - error*sd(MSJDs),
    mean(NEFFs),
    mean(NEFFs) + error*sd(NEFFs),
    mean(NEFFs) - error*sd(NEFFs),
    mean(CNEFFs),
    mean(CNEFFs) + error*sd(CNEFFs),
    mean(CNEFFs) - error*sd(CNEFFs)
  )
}
##Saving for later comparison:
setwd("../Comparison")
save(Results, file = "ResultsOptimising.Rda")

##better? Neff v R for naive:
library(EpiDataAug)
set.seed(2109)
Rs <- 15
CNEFFss <- matrix(0, nrow = Rs, ncol = 4)
runs <- 100
Total <- 15
burnin <- 10000
hyperParameters$`N+-Delta` <- parameters[[1]][[1]]
CNEFFs <- rep(0, runs)
for(R in 1:Rs){
  hyperParameters$`N+-Delta`$R <- R
  SIRS[[1]]@MCMC <- buildMCMCInternal(SIRS[[1]], hyperParameters)
  for(j in 1:runs){
    SIRS[[1]]@Model$newI <- initialNewIs[[1]]
    SIRS[[1]]@Model$Beta <- hyperParameters$`Initial Values`$Beta
    SIRS[[1]]@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    SIRS[[1]]@Model$calculate()
    SIRS[[1]]@MCMC$run(niter = burnin)
    start <- Sys.time()
    SIRS[[1]]@MCMC$run(niter = samples, reset = FALSE)
    end <- Sys.time()
    Samples <- as.matrix(SIRS[[1]]@MCMC$mvSamples)[burnin+(1:samples),]
    CNEFFs[j] <- multiESS(Samples[,c("Beta","Gamma")])/as.double(end-start, units = "secs")
  }
  error <- qt(0.975,df=runs-1)/sqrt(runs)
  CNEFFss[R,] <- c(
    R,
    mean(CNEFFs),
    mean(CNEFFs) + error*sd(CNEFFs),
    mean(CNEFFs) - error*sd(CNEFFs)
  )
  Total <- Total - 1
  print(Total)
}
width <- 1/8
ggplot() + 
  geom_point(aes(
    x = CNEFFss[,1],
    y = CNEFFss[,2]
  )) + 
  geom_line(aes(
    x = CNEFFss[,1],
    y = CNEFFss[,2]
  )) +
  geom_linerange(aes(
    x = CNEFFss[,1],
    ymax = CNEFFss[,3],
    ymin = CNEFFss[,4]
  ),
  linetype = "dashed") +
  geom_segment(aes(
    x = CNEFFss[,1] - width,
    y = CNEFFss[,4],
    xend = CNEFFss[,1] + width,
    yend = CNEFFss[,4]
  )) +
  geom_segment(aes(
    x = CNEFFss[,1] - width,
    y = CNEFFss[,3],
    xend = CNEFFss[,1] + width,
    yend = CNEFFss[,3]
  )) +
  labs(x= "R", y = "CNeff")

##Long epi comparing naive to b2:
library(EpiDataAugV2)
set.seed(210)
CNEFFs <- matrix(0, nrow = runs, ncol = 2)
runs <- 100
Total <- 15
burnin <- 15000
samples <- 1000
for(i in 1:2){
  hyperParameters$`N+-Delta` <- parameters[[2]][[c(1,3)[i]]]
  SIRS[[2]]@MCMC <- buildMCMCInternal(SIRS[[2]], hyperParameters)
  for(j in 1:runs){
    SIRS[[2]]@Model$newI <- initialNewIs[[2]]
    SIRS[[2]]@Model$Beta <- hyperParameters$`Initial Values`$Beta
    SIRS[[2]]@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    SIRS[[2]]@Model$calculate()
    SIRS[[2]]@MCMC$run(niter = burnin)
    start <- Sys.time()
    SIRS[[2]]@MCMC$run(niter = samples, reset = FALSE)
    end <- Sys.time()
    Samples <- as.matrix(SIRS[[2]]@MCMC$mvSamples)[burnin+(1:samples),]
    CNEFFs[j,i] <- multiESS(Samples[,c("Beta","Gamma")])/as.double(end-start, units = "secs")
  }
}
t.test(x = CNEFFs[,1], y = CNEFFs[,2])
qt(0.975,df=runs-1)/sqrt(runs)*sd(CNEFFs[,1])
qt(0.975,df=runs-1)/sqrt(runs)*sd(CNEFFs[,2])
