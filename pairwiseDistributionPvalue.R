#!/usr/bin/Rscript --vanilla

args <- commandArgs(TRUE)

data <- read.delim(args[1], header = FALSE, sep = " ")

data <- data.frame(t(data))

colnames(data) <- c("Observed", "Expected")
obs <- 2*na.omit(data$Observed)
null <- 2*na.omit(data$Expected)

print("Observed: 25-49%")
obs25 <- sum(obs > .25 & obs < .5)/length(obs)
print(obs25)
print("Observed: 50-74%")
obs50 <- sum(obs > .5 & obs < .75)/length(obs)
print(obs50)
print("Observed: >75%")
obs75 <- sum(obs > .75)/length(obs)
print(obs75)

print("Expected: 25-49%")
print(sum(null > .25 & null < .5)/length(null))
print("Expected: 50-74%")
print(sum(null > .5 & null < .75)/length(null))
print("Expected: >75%")
print(sum(null > .75)/length(null))

numComparisons <- length(obs)
null25 <- vector()
null50 <- vector()
null75 <- vector()
for (i in seq(1, length(null), by=numComparisons)) {
    sim <- null[c(seq(i, i+numComparisons-1))]
    null25 <- c(null25, (sum(sim > .25 & sim < .5)/length(sim)))
    null50 <- c(null50, (sum(sim > .5 & sim < .75)/length(sim)))
    null75 <- c(null75, (sum(sim > .75)/length(sim)))
}

print("P-value: 25-49%")
print(sum(null25 >= obs25)/length(null25))
print("P-value: 50-74%")
print(sum(null50 >= obs50)/length(null50))
print("P-value: >75%")
print(sum(null75 >= obs75)/length(null75))
