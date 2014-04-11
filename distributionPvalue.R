#!/usr/bin/Rscript --vanilla

args <- commandArgs(TRUE)

nullDist <- read.table(args[1], header = FALSE, nrow = 1)

null <- nullDist[1,]

sum(null > args[2])/length(null)
