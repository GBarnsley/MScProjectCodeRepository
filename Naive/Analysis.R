##Loading Data:
load("Results.Rda")
##Loading Libraries:
library(ggplot2)
##Epidemic 1:
limits <- c(0,max(Results$TDeltaMax$MSJDEpi1))
#max MSJD for colour scales
#Tile Plots
ggplot() + 
  geom_tile(data = Results$TDeltaMax[1:121,], aes(x = TMax, y = DeltaMax, fill = MSJDEpi1),
            colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi1[1:121]),],
            aes(xmin = TMax - 5,xmax = TMax + 5,
                ymin = DeltaMax - 5, ymax = DeltaMax + 5,
                fill = MSJDEpi1),
            colour = "red") +
  labs(fill = "MSJD") +
  scale_fill_continuous(type = "viridis",
                        limits = limits) +
  theme(legend.position = "none")
ggplot() + 
  geom_rect(data = Results$TDeltaMax[c(14,122:169),], aes(xmin = TMax-1.5,xmax = TMax+1.5,
                                                          ymin = DeltaMax-1.5, ymax=DeltaMax+1.5, fill = MSJDEpi1),
            colour  = "black") +
  geom_rect(data = Results$TDeltaMax[219:243,], aes(xmin = TMax - 0.5,
                                                    xmax = TMax + 0.5,
                                                    ymin = DeltaMax - 0.5,
                                                    ymax = DeltaMax + 0.5,
                                                    fill = MSJDEpi1), colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi1),],
            aes(xmin = TMax - .5,xmax = TMax + 0.5,
                ymin = DeltaMax - .5, ymax = DeltaMax + .5,
                fill = MSJDEpi1),
            colour = "red") +
  labs(fill = "MSJD", x = "TMax", y = "DeltaMax") +
  scale_fill_continuous(type = "viridis",
                        limits = limits)
#acceptance rate against msjd:
ggplot() + 
  geom_point(data = Results$TDeltaMax,
             aes(x = MSJDEpi1,
                 y = AcceptanceEpi1),
             alpha = 0.7) + 
  labs(x = "MSJD", y = "Acceptance Probability")
#runs:
ggplot(data = Results$R,
       aes(x = R,
           y = Epi1)) + 
  geom_point() + 
  geom_line() +
  labs(x = "R", y = "CNeff")

##Epidemic 2:
limits <- c(0,max(Results$TDeltaMax$MSJDEpi2))
#Tile Plot
ggplot() + 
  geom_tile(data = Results$TDeltaMax[1:121,], aes(x = TMax, y = DeltaMax, fill = MSJDEpi2),
            colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi2[1:121]),],
            aes(xmin = TMax - 5,xmax = TMax + 5,
                ymin = DeltaMax - 5, ymax = DeltaMax + 5,
                fill = MSJDEpi2),
            colour = "red") +
  labs(fill = "MSJD") +
  scale_fill_continuous(type = "viridis",
                        limits = limits) +
  theme(legend.position = "none")
ggplot() + 
  geom_rect(data = Results$TDeltaMax[c(170:218),],
            aes(xmin = TMax-1.5,
                xmax = TMax+1.5,
                ymin = DeltaMax - 1.5,
                ymax=DeltaMax+1.5,
                fill = MSJDEpi2),
            colour  = "black") +
  geom_rect(data = Results$TDeltaMax[244:268,],
            aes(xmin = TMax - 0.5,
                xmax = TMax + 0.5,
                ymin = DeltaMax - 0.5,
                ymax = DeltaMax + 0.5,
                fill = MSJDEpi2),
            colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi2),],
            aes(xmin = TMax - .5,xmax = TMax + 0.5,
                ymin = DeltaMax - .5, ymax = DeltaMax + .5,
                fill = MSJDEpi2),
            colour = "red") +
  labs(fill = "MSJD", x = "TMax", y = "DeltaMax") +
  scale_fill_continuous(type = "viridis",
                        limits = limits)
#accept against msjd
ggplot() + 
  geom_point(data = Results$TDeltaMax,
             aes(x = MSJDEpi2,
                 y = AcceptanceEpi2),
             alpha = 0.7) + 
  labs(x = "MSJD", y = "Acceptance Probability")
#runs
ggplot(data = Results$R,
       aes(x = R,
           y = Epi2)) + 
  geom_point() + 
  geom_line() +
  labs(x = "R", y = "CNeff")
##Comparing R values

#Warning this will take a long time!

#Loading libraries:
library(EpiDataAug)
library(mcmcse)
nimbleOptions(
  MCMCprogressBar = FALSE
)
load("SimulatedEpidemics.Rda")
set.seed(193)
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
    TMax = Results$TDeltaMax$TMax[which.max(Results$TDeltaMax$MSJDEpi2)],
    DeltaMax = Results$TDeltaMax$DeltaMax[which.max(Results$TDeltaMax$MSJDEpi2)]
  )
)
samples <- 1000
burnin <- 10000
thin <- 1
model <- SIR(Frequency = SimulatedEpidemics[[2]]$Frequency,
             N = SimulatedEpidemics[[2]]$Pop,
             newR = SimulatedEpidemics[[2]]$newR)
model <- initialValues(model, hyperParameters)
initialNewI <- model@Model$newI
Rs <- c(1,4)
runs <- 500
cneff <- matrix(0,ncol = length(Rs), nrow = runs)
for(j in 1:length(Rs)){
  hyperParameters$`N+-Delta`$R <- Rs[j]
  model@MCMC <- buildMCMCInternal(model, hyperParameters)
  for(i in 1:runs){
    model@Model$newI <- initialNewI
    model@Model$Beta <- hyperParameters$`Initial Values`$Beta
    model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    start <- Sys.time()
    model@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    end <- Sys.time()
    Samples <- as.matrix(model@MCMC$mvSamples)
    
    cneff[i,j] <- multiESS(Samples[,c("Beta","Gamma")])/as.double(end-start, units = "secs")
  }
}
colMeans(cneff)
t.test(x = cneff[,1], y = cneff[,2])
#variances
sqrt(colSums((cneff - matrix(c(rep(mean(cneff[,1]),runs),
                        rep(mean(cneff[,2]),runs)),
                        ncol = length(Rs)))^2)/(runs - 1))
#t and neff seperately
set.seed(1039)
runs <- 100
hyperParameters$`N+-Delta`$R <- 1
model@MCMC <- buildMCMCInternal(model, hyperParameters)
res <- matrix(0, ncol = 2, nrow = runs)
for(i in 1:runs){
    model@Model$newI <- initialNewI
    model@Model$Beta <- hyperParameters$`Initial Values`$Beta
    model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
    #model@Model$calculate()
    start <- Sys.time()
    model@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    end <- Sys.time()
    Samples <- as.matrix(model@MCMC$mvSamples)
    res[i,]  <- c(multiESS(Samples[,c("Beta","Gamma")]), as.double(end-start, units = "secs"))
}
#stds
sqrt(colSums((res - matrix(c(rep(mean(res[,1]),runs),
                               rep(mean(res[,2]),runs)),
                             ncol = 2))^2)/(runs - 1))
#Checking neff
set.seed(1050)
runs <- 100
samples <- 3000
hyperParameters$`N+-Delta`$R <- 1
model@MCMC <- buildMCMCInternal(model, hyperParameters)
res <- c()
for(i in 1:runs){
  model@Model$newI <- initialNewI
  model@Model$Beta <- hyperParameters$`Initial Values`$Beta
  model@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  #model@Model$calculate()
  start <- Sys.time()
  model@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
  end <- Sys.time()
  Samples <- as.matrix(model@MCMC$mvSamples)
  res[i]  <- multiESS(Samples[,c("Beta","Gamma")])
}
#stds
sqrt(var(res))
