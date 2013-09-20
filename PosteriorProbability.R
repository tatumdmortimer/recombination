#!/usr/bin/Rscript --vanilla

library(ggplot2)

# load the Posterior Probability Data
PostProb <- read.delim("Probabilities.out")

# reshape the data
PostProbWide <- reshape(PostProb, timevar = "Group", idvar = c("Window.Number","Start", "Stop"), direction = "wide")

# replace NA with 0
PostProbWide[is.na(PostProbWide)] <- 0

# plot all posterior probabilities
p1 <- ggplot(PostProbWide) + geom_line(aes(x=Start, y=Posterior.Probability.A), colour="purple") + geom_line(aes(x=Start,y=Posterior.Probability.K), colour="green") + geom_line(aes(x=Start, y=Posterior.Probability.G), colour="blue") + scale_colour_manual("", breaks = c("A","K","G"), values = c("purple","green","blue")) + xlab("Start of Window") + ylab("Posterior Probability")
ggsave(p1, file = "PosteriorProbabilities.tiff", width=30, height=5)

# average probabilities
PostProbWide$Posterior.Probability.mean <- (PostProbWide$Posterior.Probability.K + PostProbWide$Posterior.Probability.G + PostProbWide$Posterior.Probability.A)/3

# plot average posterior probabilities
p2 <- qplot(PostProbWide$Start, PostProbWide$Posterior.Probability.mean, geom = 'line', xlab = 'Start of Window', ylab = 'Posterior Probability')
ggsave(p2, file = "AveragePosteriorProbability.tiff", width=30, height=5)

# output table of Starts, Stops, and average Posterior Probability
Keeps <- c("Start","Stop","Posterior.Probability.mean")
Average <- PostProbWide[Keeps]
write.table(Average, file="AveragePosteriorProbabilities.txt", quote=FALSE, sep="\t", row.names=FALSE)
