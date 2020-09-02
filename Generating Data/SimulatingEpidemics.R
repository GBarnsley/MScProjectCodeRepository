##Loading Library:
library(EpiDataAug)
##Setting constant parameters:
set.seed(1020203)
#setting population
Pop <- 1000
#Using frequency based transmission
Frequency <- TRUE
#maximum length of epidemics
t.max <- 250
##Generating Epidemics:
#Note, this will be slow since the functions are built with long set up times
#Short Epidemic
Beta <- 5/2
Gamma <- 1/2
#simulating epidemic
epi1 <- simulateSIRs(Beta = Beta, Gamma = Gamma, Pop = Pop, N = 1, t.step = 1, t.max = t.max, Frequency = Frequency)
#extracting values in a list
epi1 <- list(
  Beta = Beta,
  Gamma = Gamma,
  Pop = Pop,
  Frequency = Frequency,
  newI = epi1$newI[[1]],
  newR = epi1$newR[[1]]
)
#Long Epidemic
Beta <- 1/2
Gamma <- 1/5
#simulating epidemic
epi2 <- simulateSIRs(Beta = Beta, Gamma = Gamma, Pop = Pop, N = 1, t.step = 1, t.max = t.max, Frequency = Frequency)
#extracting values in a list
epi2 <- list(
  Beta = Beta,
  Gamma = Gamma,
  Pop = Pop,
  Frequency = Frequency,
  newI = epi2$newI[[1]],
  newR = epi2$newR[[1]]
)
##Merging into one list:
SimulatedEpidemics <- list()
SimulatedEpidemics[[1]] <- epi1
SimulatedEpidemics[[2]] <- epi2
##Saving to file:
setwd("..")
save(SimulatedEpidemics, file = "Naive/SimulatedEpidemics.Rda")
save(SimulatedEpidemics, file = "Bounded1/SimulatedEpidemics.Rda")
save(SimulatedEpidemics, file = "Bounded2/SimulatedEpidemics.Rda")
##Producing Plots:
library(ggplot2)
#wrapper function to sort data and then plot
plotSIR <- function(epi){
  temp <- data.frame(
    Time = rep(0,3), 
    Population = c(epi$Pop - 1 , 1, 0), 
    State = factor(c("S", "I", "R"), levels = c("S","I","R"))
  )
  for(i in 2:(length(epi$newI)+1)){
    temp[(i-1)*3 + 1,1:2] <- c(i-1,temp[(i-2)*3 + 1,2] - epi$newI[i-1])
    temp[(i-1)*3 + 1,3] <- "S"
    temp[(i-1)*3 + 2,1:2] <- c(i-1,temp[(i-2)*3 + 2,2] + epi$newI[i-1] - epi$newR[i-1])
    temp[(i-1)*3 + 2,3] <- "I"
    temp[(i-1)*3 + 3,1:2] <- c(i-1,temp[(i-2)*3 + 3,2] + epi$newR[i-1])
    temp[(i-1)*3 + 3,3] <- "R"
  }
  print(
    ggplot(temp) +
      geom_line(aes(x = Time, y = Population, colour = State), size = 1) + 
      scale_color_manual(values = c("Green", "Red", "Blue"))
    )
}
plotSIR(epi1)
plotSIR(epi2)