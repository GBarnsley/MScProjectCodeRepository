##Loading Libraries
library(ggplot2)
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
plotNewIestimates <- function(Samples, oneInEvery, realValue, UpperLim = 50){
  plottingSamples <- seq(1, nrow(Samples), by = oneInEvery)
  timeMax <- ncol(Samples[,-c(1,2)])
  rows <- length(plottingSamples)*timeMax
  newIs <- data.frame(sample = rep(0,rows),
                      time = rep(0,rows),
                      newI = rep(0,rows))
  rowPos <- 1
  for(i in plottingSamples){
    for(j in 1:timeMax){
      newIs[rowPos,1:3] <- c(i,j,Samples[i,j+2][[1]])
      rowPos <- rowPos + 1
    }
  }
  print(
  ggplot() + 
    geom_line(data = newIs,
              aes(x = time, y = newI, group = sample),
              alpha = 0.03,
              size = 1) +
    geom_line(aes(
      x = 1:length(realValue),
      y = realValue),
      colour = "red",
      size = 1) +
    xlim(1,UpperLim) +
    labs(x = "Time", y = "I*")
  )
}
#Naive
library(EpiDataAug)
#short epi
hyperParameters$`N+-Delta` <- parameters[[1]][[1]]
model <- SIR(newR = SimulatedEpidemics[[1]]$newR, N = SimulatedEpidemics[[1]]$Pop)
model <- metropolisHastings(model, hyperParameters, samples, burnin, thin = 1)
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[1]]$newI, UpperLim = 15)
#long epi
hyperParameters$`N+-Delta` <- parameters[[2]][[1]]
model <- SIR(newR = SimulatedEpidemics[[2]]$newR, N = SimulatedEpidemics[[2]]$Pop)
model <- metropolisHastings(model, hyperParameters, samples, burnin, thin = 1)
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[2]]$newI, UpperLim = 55)
#B1
library(EpiDataAugV1)
#short epi
hyperParameters$`N+-Delta` <- parameters[[1]][[2]]
model <- SIR(newR = SimulatedEpidemics[[1]]$newR, N = SimulatedEpidemics[[1]]$Pop)
model <- metropolisHastings(model, hyperParameters, samples, burnin, thin = 1)
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[1]]$newI, UpperLim = 15)
#long epi
hyperParameters$`N+-Delta` <- parameters[[2]][[2]]
model <- SIR(newR = SimulatedEpidemics[[2]]$newR, N = SimulatedEpidemics[[2]]$Pop)
model <- metropolisHastings(model, hyperParameters, samples, burnin, thin = 1)
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[2]]$newI, UpperLim = 55)
#B2
library(EpiDataAugV2)
#short epi
hyperParameters$`N+-Delta` <- parameters[[1]][[3]]
model <- SIR(newR = SimulatedEpidemics[[1]]$newR, N = SimulatedEpidemics[[1]]$Pop)
model <- metropolisHastings(model, hyperParameters, samples, burnin, thin = 1)
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[1]]$newI, UpperLim = 15)
#long epi
hyperParameters$`N+-Delta` <- parameters[[2]][[3]]
model <- SIR(newR = SimulatedEpidemics[[2]]$newR, N = SimulatedEpidemics[[2]]$Pop)
model <- metropolisHastings(model, hyperParameters, samples, burnin, thin = 1)
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[2]]$newI, UpperLim = 55)
