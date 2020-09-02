##Loading Data:
load("OptimisationResults.Rda")
##Loading Libraries:
library(ggplot2)
##Producing plots:
#finding maximum MSJD for use in colour scales
limits <- c(0,max(Results$TDeltaMax$MSJD))
#Tile Plot
ggplot() + 
  geom_tile(data = Results$TDeltaMax[1:121,], aes(x = TMax, y = DeltaMax, fill = MSJD),
            colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJD[1:121]),],
            aes(xmin = TMax - 5,xmax = TMax + 5,
                ymin = DeltaMax - 5, ymax = DeltaMax + 5,
                fill = MSJD),
            colour = "red") +
  labs(fill = "MSJD") +
  scale_fill_continuous(type = "viridis",
                        limits = limits) +
  theme(legend.position = "none")
ggplot() + 
  geom_rect(data = Results$TDeltaMax[c(111,122:148),], aes(xmin = TMax-1.5,xmax = TMax+1.5,
                                                          ymin = DeltaMax-1.5, ymax=DeltaMax+1.5, fill = MSJD),
            colour  = "black") +
  geom_rect(data = Results$TDeltaMax[149:163,], aes(xmin = TMax - 0.5,
                                                    xmax = TMax + 0.5,
                                                    ymin = DeltaMax - 0.5,
                                                    ymax = DeltaMax + 0.5,
                                                    fill = MSJD), colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJD),],
            aes(xmin = TMax - .5,xmax = TMax + 0.5,
                ymin = DeltaMax - .5, ymax = DeltaMax + .5,
                fill = MSJD),
            colour = "red") +
  labs(fill = "MSJD", x = "TMax", y = "DeltaMax") +
  scale_fill_continuous(type = "viridis",
                        limits = limits)
#Runs against effective sample size
ggplot(data = Results$R,
       aes(x = R,
           y = Neff)) + 
  geom_point() + 
  geom_line() +
  geom_segment(aes(y = Neff + Error, yend = Neff - Error, xend = R), linetype = "dashed") +
  geom_segment(aes(y = Neff + Error, yend = Neff + Error, x = R-0.1, xend = R+0.1)) +
  geom_segment(aes(y = Neff - Error, yend = Neff - Error, x = R-0.1, xend = R+0.1)) +
  labs(x = "R", y = "Effective Sample Size")
##Performing a test on values:
#wrapper function to deduce values and return p-value
t.tester <- function(i1,i2){
  x1 <- Results$R$Neff[i1]
  x2 <- Results$R$Neff[i2]
  s1 <- Results$R$Error[i1]/qt(0.975,df=100 - 1)
  s2 <- Results$R$Error[i2]/qt(0.975,df=100 - 1)
  t <- (x1 - x2)/sqrt(s1^2 + s2^2)
  df <- (s1^2 + s2^2)^2/(s1^4/99 + s2^4/99)
  return(pt(t, df))
}
#R=6 v R=10
t.tester(6, 10)
#R=6 v R=14
t.tester(6, 14)