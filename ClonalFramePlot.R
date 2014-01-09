library(ggplot2)
library(gridExtra)

cf <- read.delim("can_clonalFrame_recombProbs.txt")
for (i in (0:max(cf$Node))) {
  assign(paste("cf.", i, sep=""), cf[cf$Node == i,])
}

plot1 <- ggplot(cf.0, aes(x = Site, y = RecombinationProbability)) + geom_point() + xlim(0,4500000) 
