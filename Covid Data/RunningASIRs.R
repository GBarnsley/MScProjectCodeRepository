##Loading Data:
load("Covid.Rda")
##Loading Libraries:
library(EpiDataAug)
library(ggplot2)
library(mcmcse)
##Setting hyper parameters:
hyperParameters <- list(
  `Initial Values` = list(
    Betas = c(1,1),
    UGamma = 0.01,
    DGamma = 0.01,
    Runs = 500000#20000
  ),
  Priors = list(
    R0 = list(
      Means = c(6.94,0.62),#c(2.63,0.62),
      SDs = c(0.62,0.54)
    ),
    RecoveryRate = list(
      Shape = 1/320,
      Rate = 1/40
    ),
    DetectionRate = list(
      Shape = 1/400,
      Rate = 1/200
    )
  ),
  RWM = list(
  ),
  `N+-Delta` = list(
    R = 2,
    TMax = 109,
    DeltaMax = 6
  )
)
samples <- 5000
burnin <- 750000
totalOverObserved <- 2.37
##Lancaster:
set.seed(10920)
#setting up model
model <- ASIR(newR = CovidData$Lancaster$newR,
              N = CovidData$Lancaster$Pop,
              Frequency = CovidData$Lancaster$Frequency,
              ChangePoint = CovidData$Lancaster$ChangePoint,
              TotalInfections = round(sum(CovidData$Lancaster$newR)*totalOverObserved))
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
  geom_line(aes(x = 1:samples, y = model@Samples[,"Betas[1]"])) + 
  labs(x = "Sample Index", y = "Beta 0")
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Betas[2]"])) + 
  labs(x = "Sample Index", y = "Beta 1")
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"DGamma"])) + 
  labs(x = "Sample Index", y = "Gamma d")
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"UGamma"])) + 
  labs(x = "Sample Index", y = "Gamma u")
#calculating R0 estimates
R0Pre <- model@Samples[,"Betas[1]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
R0Post <- model@Samples[,"Betas[2]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
quantile(R0Pre, c(0.025,0.975))
quantile(R0Post, c(0.025,0.975))
ggplot() +
  geom_line(aes(x = 1:length(model@Model$I), y = model@Model$I), colour = "red", size = 1) + 
  labs(x = "Time", y = "I")

##Colchester:
#setting up model
model <- ASIR(newR = CovidData$Colchester$newR,
              N = CovidData$Colchester$Pop,
              Frequency = CovidData$Colchester$Frequency,
              ChangePoint = CovidData$Colchester$ChangePoint,
              TotalInfections = round(sum(CovidData$Colchester$newR)*totalOverObserved))
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
  geom_line(aes(x = 1:samples, y = model@Samples[,"Betas[1]"])) + 
  labs(x = "Sample Index", y = "Beta 0")
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Betas[2]"])) + 
  labs(x = "Sample Index", y = "Beta 1")
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"DGamma"])) + 
  labs(x = "Sample Index", y = "Gamma d")
ggplot() +
  geom_line(aes(x = 1:samples, y = model@Samples[,"UGamma"])) + 
  labs(x = "Sample Index", y = "Gamma u")
#calculating R0 estimates
R0Pre <- model@Samples[,"Betas[1]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
R0Post <- model@Samples[,"Betas[2]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
quantile(R0Pre, c(0.025,0.975))
quantile(R0Post, c(0.025,0.975))
ggplot() +
  geom_line(aes(x = 1:length(model@Model$I), y = model@Model$I), colour = "red", size = 1) + 
  labs(x = "Time", y = "I")
##Colchester with a later starting time
set.seed(12)
model <- ASIR(newR = CovidData$Colchester$newR[15:length(CovidData$Colchester$newR)],
              N = CovidData$Colchester$Pop,
              Frequency = CovidData$Colchester$Frequency,
              ChangePoint = CovidData$Colchester$ChangePoint-15,
              TotalInfections = round(sum(CovidData$Colchester$newR[15:length(CovidData$Colchester$newR)])*totalOverObserved))
model <- metropolisHastings(model, hyperParameters, samples, burnin, 1)
Neff <- multiESS(model@Samples)
thin <- round(samples/Neff)
model@MCMC$run(samples*thin, thin = thin, reset = FALSE)
model@Samples <- as.matrix(model@MCMC$mvSamples)[(samples + 1):(2*samples),]
R0Pre <- model@Samples[,"Betas[1]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
R0Post <- model@Samples[,"Betas[2]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
quantile(R0Pre, c(0.025,0.975))
quantile(R0Post, c(0.025,0.975))
ggplot() +
  geom_line(aes(x = 1:length(model@Model$I), y = model@Model$I), colour = "red", size = 1) + 
  labs(x = "Time", y = "I")


##Lancashire:
set.seed(1010)
#changing NpmDelta parameters for larger size
hyperParameters$`N+-Delta`$DeltaMax <- 20
hyperParameters$`N+-Delta`$R <- 5
model <- ASIR(newR = CovidData$Lancashire$newR,
              N = CovidData$Lancashire$Pop,
              Frequency = CovidData$Lancashire$Frequency,
              ChangePoint = CovidData$Lancashire$ChangePoint,
              TotalInfections = round(sum(CovidData$Lancashire$newR)*totalOverObserved))
model <- metropolisHastings(model, hyperParameters, samples, burnin, 1)
Neff <- multiESS(model@Samples)
thin <- round(samples/Neff)
model@MCMC$run(samples*thin, thin = thin, reset = FALSE)
model@Samples <- as.matrix(model@MCMC$mvSamples)[(samples + 1):(2*samples),]
ggplot()+
  geom_line(aes(x = 1:samples, y = model@Samples[,"Betas[1]"])) + 
  labs(x = "Sample Index", y = "Beta 0")
ggplot()+
  geom_line(aes(x = 1:samples, y = model@Samples[,"Betas[2]"])) + 
  labs(x = "Sample Index", y = "Beta 1")
ggplot()+
  geom_line(aes(x = 1:samples, y = model@Samples[,"DGamma"])) + 
  labs(x = "Sample Index", y = "Gamma d")
ggplot()+
  geom_line(aes(x = 1:samples, y = model@Samples[,"UGamma"])) + 
  labs(x = "Sample Index", y = "Gamma u")
R0Pre <- model@Samples[,"Betas[1]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
R0Post <- model@Samples[,"Betas[2]"]/(model@Samples[,"DGamma"] + model@Samples[,"UGamma"])
quantile(R0Pre, c(0.025,0.975))
quantile(R0Post, c(0.025,0.975))
#producing ACF/trace plots for R0
ggplot() + 
  geom_line(aes(x = 1:samples, y = R0Pre)) + 
  labs(x = "Sample Index", y = "Pre-Lock-Down R0")
ggplot() + 
  geom_line(aes(x = 1:samples, y = R0Post)) + 
  labs(x = "Sample Index", y = "Post-Lock-Down R0")
ggplot() + 
  geom_segment(aes(x = 0:100, xend = 0:100,
                   y = acf(R0Pre, plot = FALSE, lag.max = 100)$acf), yend = rep(0,101)) + 
  geom_hline(aes(yintercept = c(-0.1,0.1)), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(x = 0:100, xend = 0:100,
                   y = acf(R0Post, plot = FALSE, lag.max = 100)$acf), yend = rep(0,101)) + 
  geom_hline(aes(yintercept = c(-0.1,0.1)), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() +
  geom_line(aes(x = 1:length(model@Model$I), y = model@Model$I), colour = "red", size = 1) + 
  labs(x = "Time", y = "I")