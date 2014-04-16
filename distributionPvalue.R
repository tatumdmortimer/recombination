#!/usr/bin/Rscript --vanilla

args <- commandArgs(TRUE)

nullDist <- read.table(args[1], header = FALSE)

null <- colMeans(nullDist)
summary(null)
sum(null > args[2])/length(null)
