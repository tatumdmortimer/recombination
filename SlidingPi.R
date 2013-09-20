#!/usr/bin/Rscript --vanilla

library(PopGenome)

#get command line arguments
args <- commandArgs(trailingOnly = TRUE)

#import alignment
alignment <- readData(args[1])

#transform to sliding window
alignment.window <- sliding.window.transform(alignment, width = 300, jump = 100, type=2, whole.data=TRUE)

#calculate pi
alignment.window <- diversity.stats(alignment.window, pi = TRUE)

out <- alignment.window@Pi

write.table(out, file = args[2], quote = FALSE, sep = "\t")

