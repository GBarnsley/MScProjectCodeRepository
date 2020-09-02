##Loading Libraries:
library(tidyverse)
library(lubridate)
##Setting up data:
#Loading csvs
Lancashire <- read.csv("Covid Lancashire.csv")
Lancaster <- read.csv("Covid Lancaster.csv")
Colchester <- read.csv("Covid Colchester.csv")
#extracting data
Lancashire <- data.frame(Time = ymd(Lancashire$date),
                         newR = Lancashire$newCasesBySpecimenDate)
Lancaster<- data.frame(Time = ymd(Lancaster$date),
                         newR = Lancaster$newCasesBySpecimenDate)
Colchester <- data.frame(Time = ymd(Colchester$date),
                         newR = Colchester$newCasesBySpecimenDate)
#adjusting start of epidemic
lockdownStart <- ymd("2020-03-23")
startingTime <- 7 #week ahead of first case
start <- Lancashire$Time[length(Lancashire$Time)] - startingTime
Lancashire$Time <- as.numeric(Lancashire$Time - start)
LancashireChange <- as.numeric(lockdownStart - start)
start <- Lancaster$Time[length(Lancaster$Time)] - startingTime
Lancaster$Time <- as.numeric(Lancaster$Time - start)
LancasterChange <- as.numeric(lockdownStart - start)
start <- Colchester$Time[length(Colchester$Time)] - startingTime
Colchester$Time <- as.numeric(Colchester$Time - start)
ColchesterChange <- as.numeric(lockdownStart - start)
#filling in 0s
Lancashire[nrow(Lancashire) + (1:(startingTime-1)),] <- matrix(c((startingTime-1):1,
                                                             rep(0,startingTime-1)),
                                                           ncol = 2,
                                                           nrow = startingTime-1)
Lancaster[nrow(Lancaster) + (1:(startingTime-1)),] <- matrix(c((startingTime-1):1,
                                                                 rep(0,startingTime-1)),
                                                               ncol = 2,
                                                               nrow = startingTime-1)
Colchester[nrow(Colchester) + (1:(startingTime-1)),] <- matrix(c((startingTime-1):1,
                                                                 rep(0,startingTime-1)),
                                                               ncol = 2,
                                                               nrow = startingTime-1)

#flipping to be in time order
Lancashire$Time <- rev(Lancashire$Time)
Lancashire$newR <- rev(Lancashire$newR)
Lancaster$Time <- rev(Lancaster$Time)
Lancaster$newR <- rev(Lancaster$newR)
Colchester$Time <- rev(Colchester$Time)
Colchester$newR <- rev(Colchester$newR)
#plots
ggplot(Lancashire, aes(x = Time, y = newR)) + 
  geom_line(size = 1) +
  labs(y = "R*")
ggplot(Lancaster, aes(x = Time, y = newR)) + 
  geom_line(size = 1) +
  labs(y = "R*")
ggplot(Colchester, aes(x = Time, y = newR)) + 
  geom_line(size = 1) +
  labs(y = "R*")
##Putting into list format
CovidData <- list(
  Lancashire = list(
    newR = Lancashire$newR,
    Frequency = TRUE,
    Pop = 1.21*10^6,
    ChangePoint = LancashireChange
  ),
  Lancaster = list(
    newR = Lancaster$newR,
    Frequency = TRUE,
    Pop = 138375,
    ChangePoint = LancasterChange
  ),
  Colchester = list(
    newR = Colchester$newR,
    Frequency = TRUE,
    Pop = 194706,
    ChangePoint = ColchesterChange
  )
)
##Simulated Epidemic
set.seed(101)
Pop <- round((CovidData$Lancaster$Pop + CovidData$Colchester$Pop)/2)
maxT <- max(c(length(CovidData$Lancaster$newR), length(CovidData$Colchester$newR)))
Gamma <- 1/7
Beta <- 105/700
Frequency <- TRUE
prob <- function(x){
  return(1 - exp(-x))
}
epiLength <- 1
epiMax <- 1
epiSize <- 1
while(epiLength < 100 | epiLength > 300 | epiSize < 500){
  epidemic <- data.frame(S = Pop - 1,
                         newI = NA,
                         I = 1,
                         newR = NA,
                         R = 0)
  rowPos <- 1
  while(epidemic$I[rowPos] > 0){
    epidemic$newI[rowPos] <- rbinom(1, epidemic$S[rowPos], prob(Beta*epidemic$I[rowPos]/(Pop^Frequency)))
    epidemic$newR[rowPos] <- rbinom(1, epidemic$I[rowPos], prob(Gamma))
    rowPos <- rowPos + 1
    epidemic[rowPos,] <- c(epidemic$S[rowPos-1] - epidemic$newI[rowPos - 1],
                           NA,
                           epidemic$I[rowPos-1] + epidemic$newI[rowPos - 1] - epidemic$newR[rowPos - 1],
                           NA,
                           epidemic$R[rowPos-1] + epidemic$newR[rowPos - 1])
  }
  #conditions to try and ensure it follows the correct layout
  epiLength <- nrow(epidemic)
  epiMax <- max(epidemic$newR, na.rm = TRUE)
  epiSize <- sum(epidemic$newR, na.rm = TRUE)
}
ggplot(data = epidemic, aes(x = 1:nrow(epidemic), y = newR)) +
        geom_line(size = 1) +
  labs(x = "Time", y= "R*")

##Saving
CovidData$Simulated <- list(
  Beta = Beta,
  Gamma = Gamma,
  newI = na.omit(epidemic$newI),
  newR = na.omit(epidemic$newR),
  Frequency = Frequency,
  Pop = Pop
)
save(CovidData, file = "Covid.Rda")
