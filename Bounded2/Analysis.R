#Please see naive or Covid version for more detailed annotations
load("Results.Rda")
library(ggplot2)
##Epidemic 1
limits <- c(0,max(Results$TDeltaMax$MSJDEpi1))
#Tile Plot
ggplot() + 
  geom_tile(data = Results$TDeltaMax[1:121,], aes(x = TMax, y = DeltaMax, fill = MSJDEpi1),
            colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi1[1:121]),],
            aes(xmin = TMax - 5,xmax = TMax + 5,
                ymin = DeltaMax - 5, ymax = DeltaMax + 5,
                fill = MSJDEpi1),
            colour = "red") +
  labs(fill = "MSJD") +
  scale_fill_continuous(type = "viridis",
                        limits = limits) +
  theme(legend.position = "none")
ggplot() + 
  geom_rect(data = Results$TDeltaMax[c(13,122:169),], aes(xmin = TMax-1.5,xmax = TMax+1.5,
                                                          ymin = DeltaMax-1.5, ymax=DeltaMax+1.5, fill = MSJDEpi1),
            colour  = "black") +
  geom_rect(data = Results$TDeltaMax[219:243,], aes(xmin = TMax - 0.5,
                                                    xmax = TMax + 0.5,
                                                    ymin = DeltaMax - 0.5,
                                                    ymax = DeltaMax + 0.5,
                                                    fill = MSJDEpi1), colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi1),],
            aes(xmin = TMax - .5,xmax = TMax + 0.5,
                ymin = DeltaMax - .5, ymax = DeltaMax + .5,
                fill = MSJDEpi1),
            colour = "red") +
  labs(fill = "MSJD", x = "TMax", y = "DeltaMax") +
  scale_fill_continuous(type = "viridis",
                        limits = limits)
#accept against msjd
ggplot() + 
  geom_point(data = Results$TDeltaMax,
             aes(x = MSJDEpi1,
                 y = AcceptanceEpi1),
             alpha = 0.7) + 
  labs(x = "MSJD", y = "Acceptance Probability")
#runs
ggplot(data = Results$R,
       aes(x = R,
           y = Epi1)) + 
  geom_point() + 
  geom_line() +
  labs(x = "R", y = "CNeff")

##Epidemic 2
limits <- c(0,max(Results$TDeltaMax$MSJDEpi2))
#Tile Plot
ggplot() + 
  geom_tile(data = Results$TDeltaMax[1:121,], aes(x = TMax, y = DeltaMax, fill = MSJDEpi2),
            colour = "black") + 
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi2[1:121]),],
            aes(xmin = TMax - 5,xmax = TMax + 5,
                ymin = DeltaMax - 5, ymax = DeltaMax + 5,
                fill = MSJDEpi2),
            colour = "red") +
  labs(fill = "MSJD") +
  scale_fill_continuous(type = "viridis",
                        limits = limits) +
  theme(legend.position = "none")
ggplot() + 
  geom_rect(data = Results$TDeltaMax[c(170:218),],
            aes(xmin = TMax-1.5,
                xmax = TMax+1.5,
                ymin = DeltaMax - 1.5,
                ymax=DeltaMax+1.5,
                fill = MSJDEpi2),
            colour  = "black") +
  geom_rect(data = Results$TDeltaMax[244:268,],
            aes(xmin = TMax - 0.5,
                xmax = TMax + 0.5,
                ymin = DeltaMax - 0.5,
                ymax = DeltaMax + 0.5,
                fill = MSJDEpi2),
            colour = "black") +
  geom_rect(data = Results$TDeltaMax[which.max(Results$TDeltaMax$MSJDEpi2),],
            aes(xmin = TMax - .5,xmax = TMax + 0.5,
                ymin = DeltaMax - .5, ymax = DeltaMax + .5,
                fill = MSJDEpi2),
            colour = "red") +
  labs(fill = "MSJD", x = "TMax", y = "DeltaMax") +
  scale_fill_continuous(type = "viridis",
                        limits = limits)
#accept against msjd
ggplot() + 
  geom_point(data = Results$TDeltaMax,
             aes(x = MSJDEpi2,
                 y = AcceptanceEpi2),
             alpha = 0.7) + 
  labs(x = "MSJD", y = "Acceptance Probability")
#runs
ggplot(data = Results$R,
       aes(x = R,
           y = Epi2)) + 
  geom_point() + 
  geom_line() +
  labs(x = "R", y = "CNeff")