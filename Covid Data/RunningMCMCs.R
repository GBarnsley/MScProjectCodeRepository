##Loading Data:
load("Covid.Rda")
##Loading Library:
library(EpiDataAug)
library(ggplot2)
library(mcmcse)
##SIR with R0 Prior:
SIRclass <- setClass(
  "SIR",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY"
  )
)
newISIRclass <- setClass(
  "iSIR",
  contains = "SIR"
)
metropolisHastings <- function(epiModel,
                               hyperParameters,
                               samples = 1000,
                               burnin = 500,
                               thin = 10){
  epiModel <- initialValues(epiModel, hyperParameters)
  epiModel@MCMC <- buildMCMCInternal(epiModel, hyperParameters)
  epiModel@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
  epiModel@Samples <- as.matrix(epiModel@MCMC$mvSamples)
  return(epiModel)
}
SIR <- function(S = NULL,
                I = NULL,
                R = NULL,
                newI = NULL,
                newR = NULL,
                N = NULL,
                Beta = NULL,
                Gamma = NULL,
                t.step = 1,
                Frequency = TRUE){
  #calculating initial values from given dataS
  if(is.null(N)){
    print("Error: N must be specified")
    return(NA)
  }
  if(is.null(S)&!is.null(newI)){
    S <- rep(N-1, length(newI) + 1)
    for(i in 2:length(S)){
      S[i] <- N - 1 - sum(newI[1:(i-1)])
    }
  }
  if(is.null(R)&!is.null(newR)){
    R <- rep(0, length(newR) + 1)
    for(i in 2:length(R)){
      R[i] <- sum(newR[1:(i-1)])
    }
  }
  if(is.null(I)&!is.null(S)&!is.null(R)){
    I <- N - S - R
  }
  if(is.null(newI)&!is.null(S)){
    newI <- -diff(S)
  }
  if(is.null(newR)&!is.null(R)){
    newR <- diff(R)
  }
  if(Frequency){
    Frequency <- 1
  }
  else{
    Frequency <- 0
  }
  tempCode <- nimbleCode({
    # Set priors
    Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
    Beta ~ dnorm(mean = R0Mean*Pop^(Frequency == 0)*Gamma, sd = R0SD) #new Prior!!
    # likelihood
    S[1] <- Pop - 1
    I[1] <- 1
    for(i in 1:TimePeriod){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Beta*t.step/(Pop^Frequency)))
      newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
      S[i+1] <- S[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newR[i]
    }
  })
  return(newISIRclass(
    Model = compileNimble(
      nimbleModel(
        code = tempCode,
        constants = list(TimePeriod = length(newR)),
        data = list(newR = newR,
                    t.step = t.step,
                    Pop = N,
                    R0Mean = 1,
                    R0SD = 1,
                    GammaShape = 1,
                    GammaRate = 1,
                    Frequency = Frequency),
        inits = list(Beta = 1,
                     Gamma = 1,
                     newI = rep(0, length(newR))
        ),
        calculate = FALSE
      )
    ),
    MCMC = NA,
    Samples = NA
  )
  )
}
#changing initialisation method to add new priors
initialValues.iSIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$R0Mean <- hyperParameters$Priors$R0$Mean
  epiModel@Model$R0SD <- hyperParameters$Priors$R0$SD
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  epiModel@Model$newI <- rep(0, length(epiModel@Model$newR))
  epiModel@Model$newI[1] <- sum(epiModel@Model$newR) - 1
  mcmc <- configureMCMC(epiModel@Model, nodes = NULL)
  mcmc$addSampler(target = "newI",
                  type = sampler,
                  control = list(
                    TMax = 20,
                    DeltaMax = 20,
                    R = hyperParameters$`Initial Values`$Runs
                  ))
  mcmc <- buildMCMC(
    mcmc
  )
  mcmc <- compileNimble(mcmc, project = epiModel@Model, resetFunctions = TRUE)
  mcmc$run(1)
  return(
    epiModel
  )
}
set.seed(01001)
acfLagMax <- 100
##Constant Parameters:
samples <- 1000
burnin <- 20000
hyperParameters <- list(
  `Initial Values` = list(
    Beta = 1,
    Gamma = 0.01,
    Runs = 2000
  ),
  Priors = list(
    R0 = list(
      Mean = 0.62,
      SD = 0.26
    ),
    Gamma = list(
      Shape = 1/400,
      Rate = 1/200
    )
  ),
  RWM = list(
    propCov = matrix(c(1,0,0,0.01), nrow = 2)
  ),
  `N+-Delta` = list(
    R = 6,
    TMax = 109,
    DeltaMax = 2
  )
)

##Lancaster:
#initial run
model <- SIR(newR = CovidData$Lancaster$newR,
             Frequency = CovidData$Lancaster$Frequency,
             N = CovidData$Lancaster$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
Neff <- multiESS(model@Samples[,c(1,2)])
thin <- round(samples/Neff)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = thin
)
#producing plots
ggplot() + 
  geom_line(aes(x = 1:(samples), y = model@Samples[,"Beta"])) +
  labs(x = "Sample Index", y = "Beta")
ggplot() + 
  geom_line(aes(x = 1:(samples), y = model@Samples[,"Gamma"])) +
  labs(x = "Sample Index", y = "Gamma")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(model@Samples[,"Beta"], lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(model@Samples[,"Gamma"], lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
#estimates for R0
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
ggplot() + 
  geom_line(aes(x = 1:(samples), y = R0)) +
  labs(x = "Sample Index", y = "R0")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(R0, lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_line(aes(x = 1:length(model@Model$I), y =  model@Model$I), colour = "red", size = 1) +
  labs(x = "Time", y = "I")
mean(R0)
quantile(R0, c(0.025,0.975))
#correlation
ggplot() +
  geom_point(aes(y = model@Samples[,"Beta"], x = model@Samples[,"Gamma"]),
             alpha = 0.5) + 
  geom_abline(aes(intercept = 0, slope = mean(R0)), colour = "red", size = 1,
              alpha = 0.5) +
  labs(x = "Gamma", y = "Beta")

##Colchester:
#initial run
model <- SIR(newR = CovidData$Colchester$newR,
             Frequency = CovidData$Colchester$Frequency,
             N = CovidData$Colchester$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
Neff <- multiESS(model@Samples[,c(1,2)])
thin <- round(samples/Neff)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = thin
)
ggplot() + 
  geom_line(aes(x = 1:(samples), y = model@Samples[,"Beta"])) +
  labs(x = "Sample Index", y = "Beta")
ggplot() + 
  geom_line(aes(x = 1:(samples), y = model@Samples[,"Gamma"])) +
  labs(x = "Sample Index", y = "Gamma")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(model@Samples[,"Beta"], lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(model@Samples[,"Gamma"], lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
ggplot() + 
  geom_line(aes(x = 1:(samples), y = R0)) +
  labs(x = "Sample Index", y = "R0")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(R0, lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_line(aes(x = 1:length(model@Model$I), y =  model@Model$I), colour = "red", size = 1) +
  labs(x = "Time", y = "I")
mean(R0)
quantile(R0, c(0.025,0.975))

##Lancashire:
#changin NpmDelta parameters
hyperParameters$`N+-Delta`$DeltaMax <- 25
hyperParameters$`N+-Delta`$R <- 10
burnin <- 75000 ##Takes a while to run!!
set.seed(1020)
#initial run
model <- SIR(newR = CovidData$Lancashire$newR,
             Frequency = CovidData$Lancashire$Frequency,
             N = CovidData$Lancashire$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
Neff <- multiESS(model@Samples[,c(1,2)])
thin <- round(samples/Neff)
model@MCMC$run(samples*thin, thin = thin, reset = FALSE)
model@Samples <- as.matrix(model@MCMC$mvSamples)[(samples + 1):(2*samples),]
ggplot() + 
  geom_line(aes(x = 1:(samples), y = model@Samples[,"Beta"])) +
  labs(x = "Sample Index", y = "Beta")
ggplot() + 
  geom_line(aes(x = 1:(samples), y = model@Samples[,"Gamma"])) +
  labs(x = "Sample Index", y = "Gamma")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(model@Samples[,"Beta"], lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(model@Samples[,"Gamma"], lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
ggplot() + 
  geom_line(aes(x = 1:(samples), y = R0)) +
  labs(x = "Sample Index", y = "R0")
ggplot() + 
  geom_segment(aes(x = 0:acfLagMax, xend = 0:acfLagMax, y = rep(0, acfLagMax + 1),
                   yend = acf(R0, lag.max = acfLagMax, plot = FALSE)$acf)) +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") +
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_line(aes(x = 1:length(model@Model$I), y =  model@Model$I), colour = "red", size = 1) +
  labs(x = "Time", y = "I")
mean(R0)
quantile(R0, c(0.025,0.975))


##Colchester starting only 3 days before first recovery
hyperParameters$`N+-Delta`$DeltaMax <- 2
hyperParameters$`N+-Delta`$R <- 6
burnin <- 20000
#initial run
model <- SIR(newR = CovidData$Colchester$newR[5:length(CovidData$Colchester$newR)],
             Frequency = CovidData$Colchester$Frequency,
             N = CovidData$Colchester$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
Neff <- multiESS(model@Samples[,c(1,2)])
thin <- round(samples/Neff)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = thin
)
ggplot() + 
  geom_line(aes(x = 1:length(model@Model$I), y =  model@Model$I), colour = "red", size = 1) +
  labs(x = "Time", y = "I")
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
mean(R0)
quantile(R0, c(0.025,0.975))

##optimising R for Lancashire
bestNeff <- 0
hyperParameters$`N+-Delta`$DeltaMax <- 25
burnin <- 75000
model <- SIR(newR = CovidData$Lancashire$newR,
             Frequency = CovidData$Lancashire$Frequency,
             N = CovidData$Lancashire$Pop)
for(R in 1:20){
  hyperParameters$`N+-Delta`$R <- R
  ##Takes a while to run!!
  set.seed(1020)
  #initial run
  model <- metropolisHastings(
    model,
    hyperParameters,
    samples = 10000,
    burnin = burnin,
    thin = 1
  )
  Neff <- multiESS(model@Samples[,c(1,2)])
  if(Neff > bestNeff){
    bestNeff <- Neff
    bestR <- R
    bestSamples <- model@Samples
  }
}
R0 <- model@Samples[,"Beta"]/model@Samples[,"Gamma"]
ggplot() + 
  geom_line(aes(x = 1:10000, y = R0)) +
  labs(x = "Sample Index", y = "R0")
mean(R0)
quantile(R0, c(0.025,0.975))
