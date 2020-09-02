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
    Runs = 5000
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
    R = 6,
    TMax = 109,
    DeltaMax = 2
  )
)
samples <- 5000
burnin <- 750000
##Lancaster:
set.seed(10920)
#setting up model
model <- SEIR(newR = CovidData$Lancaster$newR,
              N = CovidData$Lancaster$Pop,
              Frequency = CovidData$Lancaster$Frequency)
#running algorithm
model <- metropolisHastings(model, hyperParameters, samples, burnin, 1)
#calcuting the effective sample size
Neff <- multiESS(model@Samples)
#seting thinning interval
thin <- round(samples/Neff)
#generating new samples
model@MCMC$run(samples*thin, thin = thin, reset = FALSE)
#extracting new samples
model@Samples <- as.matrix(model@MCMC$mvSamples)[(samples + 1):(2*samples),]
#making plots
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"k"])) + 
  labs(x = "Sample Index", y = "k")
#calculating R0 estimates
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
quantile(R0, c(0.025,0.975))
quantile(model@Samples[,"k"], c(0.025,0.975))

##Colchester:
#setting up model
model <- SEIR(newR = CovidData$Colchester$newR,
              N = CovidData$Colchester$Pop,
              Frequency = CovidData$Colchester$Frequency)
#running algorithm
model <- metropolisHastings(model, hyperParameters, samples, burnin, 1)
#calcuting the effective sample size
Neff <- multiESS(model@Samples)
#seting thinning interval
thin <- round(samples/Neff)
#generating new samples
model@MCMC$run(samples*thin, thin = thin, reset = FALSE)
#extracting new samples
model@Samples <- as.matrix(model@MCMC$mvSamples)[(samples + 1):(2*samples),]
#calculating R0 estimates
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
quantile(R0, c(0.025,0.975))
quantile(model@Samples[,"k"], c(0.025,0.975))


##Lancashire:
set.seed(1010)
#changing NpmDelta parameters for larger size
hyperParameters$`N+-Delta`$DeltaMax <- 20
hyperParameters$`N+-Delta`$R <- 5
model <- SEIR(newR = CovidData$Lancashire$newR,
              N = CovidData$Lancashire$Pop,
              Frequency = CovidData$Lancashire$Frequency)
#running algorithm
model <- metropolisHastings(model, hyperParameters, samples, burnin, 1)
#calcuting the effective sample size
Neff <- multiESS(model@Samples)
#seting thinning interval
thin <- round(samples/Neff)
#generating new samples
model@MCMC$run(samples*thin, thin = thin, reset = FALSE)
#extracting new samples
model@Samples <- as.matrix(model@MCMC$mvSamples)[(samples + 1):(2*samples),]
#calculating R0 estimates
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
quantile(R0, c(0.025,0.975))
quantile(model@Samples[,"k"], c(0.025,0.975))
