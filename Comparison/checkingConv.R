##Libraries
library(ggplot2)
set.seed(193021)
##Loading Data
setwd("C:/Users/gregb/Documents/Data Science/Project/R/Bounded2")
load("Results.Rda")
B2 <- Results
##Checking Burnin on epi 2
load("SimulatedEpidemics.Rda")
library(EpiDataAugV2)
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
    R = 1,
    TMax = B2$TDeltaMax$TMax[which.max(B2$TDeltaMax$MSJDEpi2)],
    DeltaMax = B2$TDeltaMax$DeltaMax[which.max(B2$TDeltaMax$MSJDEpi2)],
    Bounded = TRUE
  )
)
samples <- 10000 * 3
burnin <- 0
thin <- 1
model <- metropolisHastings(SIR(Frequency = SimulatedEpidemics[[2]]$Frequency,
                                N = SimulatedEpidemics[[2]]$Pop,
                                newR = SimulatedEpidemics[[2]]$newR),
                            hyperParameters,
                            samples,
                            burnin,
                            thin)
ggplot() +
  geom_line(
    aes(y = model@Samples[,"Beta"], x = 1:samples),
    alpha = 0.8
  ) +
  geom_line(
    aes(y = rep(1/2,2), x = c(1, samples)),
    colour = "red"
  ) +
  geom_vline(
    aes(xintercept = 10000),
    linetype = "dashed"
  ) +
  geom_vline(
    aes(xintercept = 11000),
    linetype = "dashed"
  ) +
  labs(x = "Iteration", y = "Beta")
ggplot() +
  geom_line(
    aes(y = model@Samples[,"Gamma"], x = 1:samples),
    alpha = 0.8
  ) +
  geom_line(
    aes(y = rep(1/5,2), x = c(1, samples)),
    colour = "red"
  ) +
  geom_vline(
    aes(xintercept = 10000),
    linetype = "dashed"
  ) +
  geom_vline(
    aes(xintercept = 11000),
    linetype = "dashed"
  ) +
  labs(x = "Iteration", y = "Gamma")
