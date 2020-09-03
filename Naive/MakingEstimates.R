##Loading Libraries:
library(ggplot2)
library(mcmcse)
library(EpiDataAug)
##Loading Data:
load("SimulatedEpidemics.Rda")
set.seed(01001)
acfLagMax <- 100
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
##Constant Parameters:
samples <- 1000
burnin <- 15000
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
##Short Epidemic:
hyperParameters$`N+-Delta` <- list(
  R = 1,
  TMax = 4,
  DeltaMax = 19
)
#initial run
model <- SIR(newR = SimulatedEpidemics[[1]]$newR,
             Frequency = SimulatedEpidemics[[1]]$Frequency,
             N = SimulatedEpidemics[[1]]$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
#traceplots
ggplot() + 
  geom_line(aes(x = c(1,samples), y = rep(SimulatedEpidemics[[1]]$Beta,2)),
            colour = "red") +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Beta"])) +
  labs(x = "Sample Index", y = "Beta") 
ggplot() + 
  geom_line(aes(x = c(1,samples), y = rep(SimulatedEpidemics[[1]]$Gamma,2)),
            colour = "red") +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Gamma"])) +
  labs(x = "Sample Index", y = "Gamma")
#acf plots
ggplot() + 
  geom_segment(aes(
    x = 0:acfLagMax,
    xend = 0:acfLagMax,
    y = acf(model@Samples[,"Beta"], lag.max = acfLagMax, plot = FALSE)$acf,
    yend = rep(0, acfLagMax + 1)
  )) + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") + 
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(
    x = 0:acfLagMax,
    xend = 0:acfLagMax,
    y = acf(model@Samples[,"Gamma"], lag.max = acfLagMax, plot = FALSE)$acf,
    yend = rep(0, acfLagMax + 1)
  )) + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") + 
  labs(x = "Lag", y = "ACF")
#finding effective sample size
neff <- multiESS(model@Samples[,c("Beta","Gamma")])
thin <- round(samples/neff)
#reruning
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = thin
)
ggplot() + 
  geom_line(aes(x = c(1,samples), y = rep(SimulatedEpidemics[[1]]$Beta,2)),
            colour = "red") +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Beta"])) +
  labs(x = "Sample Index", y = "Beta") 
ggplot() + 
  geom_line(aes(x = c(1,samples), y = rep(SimulatedEpidemics[[1]]$Gamma,2)),
            colour = "red") +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Gamma"])) +
  labs(x = "Sample Index", y = "Gamma")
#acf plots
ggplot() + 
  geom_segment(aes(
    x = 0:acfLagMax,
    xend = 0:acfLagMax,
    y = acf(model@Samples[,"Beta"], lag.max = acfLagMax, plot = FALSE)$acf,
    yend = rep(0, acfLagMax + 1)
  )) + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") + 
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(
    x = 0:acfLagMax,
    xend = 0:acfLagMax,
    y = acf(model@Samples[,"Gamma"], lag.max = acfLagMax, plot = FALSE)$acf,
    yend = rep(0, acfLagMax + 1)
  )) + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") + 
  labs(x = "Lag", y = "ACF")
#newI plot
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[1]]$newI, UpperLim = 15)

##Summaries:
colMeans(model@Samples[,c(1,2)])
#credible interval
quantile(model@Samples[,"Beta"], c(0.025,0.975))
quantile(model@Samples[,"Gamma"], c(0.025,0.975))

##Long Epidemic:
set.seed(219090)
hyperParameters$`N+-Delta` <- list(
  R = 1,
  TMax = 70,
  DeltaMax = 5
)
#initial run
model <- SIR(newR = SimulatedEpidemics[[2]]$newR,
             Frequency = SimulatedEpidemics[[2]]$Frequency,
             N = SimulatedEpidemics[[2]]$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
neff <- multiESS(model@Samples[,c(1,2)])
thin <- round(samples/neff)
#final
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = thin
)
ggplot() + 
  geom_line(aes(x = c(1,samples), y = rep(SimulatedEpidemics[[2]]$Beta,2)),
            colour = "red") +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Beta"])) +
  labs(x = "Sample Index", y = "Beta") 
ggplot() + 
  geom_line(aes(x = c(1,samples), y = rep(SimulatedEpidemics[[2]]$Gamma,2)),
            colour = "red") +
  geom_line(aes(x = 1:samples, y = model@Samples[,"Gamma"])) +
  labs(x = "Sample Index", y = "Gamma")
#acf plots
ggplot() + 
  geom_segment(aes(
    x = 0:acfLagMax,
    xend = 0:acfLagMax,
    y = acf(model@Samples[,"Beta"], lag.max = acfLagMax, plot = FALSE)$acf,
    yend = rep(0, acfLagMax + 1)
  )) + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") + 
  labs(x = "Lag", y = "ACF")
ggplot() + 
  geom_segment(aes(
    x = 0:acfLagMax,
    xend = 0:acfLagMax,
    y = acf(model@Samples[,"Gamma"], lag.max = acfLagMax, plot = FALSE)$acf,
    yend = rep(0, acfLagMax + 1)
  )) + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.1), linetype = "dashed") + 
  labs(x = "Lag", y = "ACF")
#newI plot
plotNewIestimates(model@Samples, 1, SimulatedEpidemics[[2]]$newI, UpperLim = 55)

##Summaries:
colMeans(model@Samples[,c(1,2)])
#credible interval
quantile(model@Samples[,"Beta"], c(0.025,0.975))
quantile(model@Samples[,"Gamma"], c(0.025,0.975))


##R
set.seed(219090)
hyperParameters$`N+-Delta` <- list(
  R = 5,
  TMax = 70,
  DeltaMax = 5
)
#initial run
model <- SIR(newR = SimulatedEpidemics[[2]]$newR,
             Frequency = SimulatedEpidemics[[2]]$Frequency,
             N = SimulatedEpidemics[[2]]$Pop)
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = 1
)
neff <- multiESS(model@Samples[,c(1,2)])
thin <- round(samples/neff)
#final
model <- metropolisHastings(
  model,
  hyperParameters,
  samples = samples,
  burnin = burnin,
  thin = thin
)
#credible interval
quantile(model@Samples[,"Beta"], c(0.025,0.975))
quantile(model@Samples[,"Gamma"], c(0.025,0.975))
