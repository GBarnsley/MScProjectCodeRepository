##Loading Libraries:
library(EpiDataAug)
library(ggplot2)
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
    R = 6,
    TMax = 40,
    DeltaMax = 10
  )
)
samples <- 1000
burnin <- 750000
plotNewIestimates <- function(Samples, oneInEvery, realValue, UpperLim = 50, Name = "I*"){
  plottingSamples <- seq(1, nrow(Samples), by = oneInEvery)
  timeMax <- ncol(Samples)
  rows <- length(plottingSamples)*timeMax
  newIs <- data.frame(sample = rep(0,rows),
                      time = rep(0,rows),
                      newI = rep(0,rows))
  rowPos <- 1
  for(i in plottingSamples){
    for(j in 1:timeMax){
      newIs[rowPos,1:3] <- c(i,j,Samples[i,j][[1]])
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
      labs(x = "Time", y = Name)
  )
}
##Simulating
set.seed(1001)
Beta <- 1/2
K <- 1/5
Gamma <- 1/3
pass <- F
while(!pass){
  S <- c(1000 -1) 
  E <- c(0)
  I <- c(1)
  newE <- c()
  newI <- c()
  newR <- c()
  time <- 1
  while(I[time] > 0){
    newE[time] <- rbinom(1,S[time],probGen(Beta*I[time]/1000))
    newI[time] <- rbinom(1,E[time],probGen(K))
    newR[time] <- rbinom(1,I[time],probGen(Gamma))
    time <- time + 1
    S[time] <- S[time - 1] - newE[time-1]
    E[time] <- E[time - 1] + newE[time-1] - newI[time-1]
    I[time] <- I[time - 1] + newI[time-1] - newR[time-1]
  }
  if(S[time] < 500){
    pass <- TRUE
  }
}

#running
newIs <- matrix(0,ncol = length(newR), nrow = samples)
newEs <- matrix(0,ncol = length(newR), nrow = samples)
ks <- rep(0, samples)
model <- SEIR(newR = newR, N = 1000)
model <- metropolisHastings(model, hyperParameters, burnin, burnin = 0, thin = 1)
for(i in 1:samples){
  model@MCMC$run(1, reset = FALSE)
  newIs[i,] <- model@Model$newI
  newEs[i,] <- model@Model$newE
  ks[i] <- model@Model$k
}
for(i in 1:samples){
  model@MCMC$run(1, reset = FALSE)
  newIs[i,] <- model@Model$newI
  newEs[i,] <- model@Model$newE
  ks[i] <- model@Model$k
}
ggplot() +
  geom_line(aes(x = 1:samples, y = ks)) +
  geom_hline(yintercept = K, colour = "red") +
  labs(x = "Sample", y = "K")
plotNewIestimates(newIs, 5, newI, 118)
plotNewIestimates(newEs, 5, newE, 118, Name = "E*")
  
