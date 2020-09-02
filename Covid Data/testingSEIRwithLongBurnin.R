##Loading Data:
load("Covid.Rda")
##Loading Libraries:
library(EpiDataAug)
library(ggplot2)
library(mcmcse)
##Setting hyper parameters:
hyperParameters <- list(
  `Initial Values` = list(
    Beta = 1,
    Gamma = 0.01,
    k = 1/7,
    Runs = 10000
  ),
  Priors = list(
    R0 = list(
      Mean = 0.62,
      SD = 0.26
    ),
    Gamma = list(
      Shape = 1/400,
      Rate = 1/200
    ),
    k = list(
      Shape = 200/(51^2), #2/(7^2)
      Rate =  20/51 #2/7
    )
  ),
  RWM = list(
    propCov = matrix(c(1,0,0,
                       0,0.01,0,
                       0,0,1), nrow = 3)
  ),
  `N+-Delta` = list(
    R = 1,
    TMax = 109,
    DeltaMax = 20
  )
)
samples <- 5000
burnin <- 5000000
##Lancaster:
set.seed(10920)
#setting up model
model <- SEIR(newR =CovidData$Lancaster$newR,
              N = CovidData$Lancaster$Pop,
              Frequency = CovidData$Lancaster$Frequency)
#running algorithm
model <- metropolisHastings(model, hyperParameters, samples, burnin, 10)
#plot
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"k"])) + 
  labs(x = "Sample Index", y = "k")
#calculating R0 estimates
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
quantile(R0, c(0.025,0.975))
quantile(model@Samples[,"k"], c(0.025,0.975))