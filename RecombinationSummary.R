#!/usr/bin/Rscript --vanilla

library(ggplot2)

#get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# load the SNP data
SNPtable <- read.delim(args[1])

# load the recombination data
Pi <- read.delim("Pi.out")
phi <- read.delim("PHI.out")
brat <- read.delim("brat.out")
gene <- read.delim("geneconv.out")

# remove unwanted columns
keeps <- c("Pos","Entry","DNA")
SNPtable <- SNPtable[keeps]

# add Groups column with appropriate labels
SNPtable$Group <- "Global" 
SNPtable$Group[SNPtable$Entry == "SUM002NoPT"] <- "DS6Quebec" 
SNPtable$Group[SNPtable$Entry == "SUM003NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM004NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM005NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM006NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM007NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM008NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM009NoPT"] <- "DS6Quebec"
SNPtable$Group[SNPtable$Entry == "SUM001NoPT"] <- "C"
SNPtable$Group[SNPtable$Entry == "SUM010NoPT"] <- "C"
SNPtable$Group[SNPtable$Entry == "SUM011NoPT"] <- "C"
SNPtable$Group[SNPtable$Entry == "SUM012NoPT"] <- "C"
SNPtable$Group[SNPtable$Entry == "17_an"] <- "jm"
SNPtable$Group[SNPtable$Entry == "345_an"] <- "jm"
SNPtable$Group[SNPtable$Entry == "72"] <- "jm"
SNPtable$Group[SNPtable$Entry == "87"] <- "jm"

SNPtable$Group[SNPtable$Entry == "M19A"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "M20A"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "BC13"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "BC18"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "BC41"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "BC42"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "BC191"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "BC192"] <- "transconjugant"
SNPtable$Group[SNPtable$Entry == "NC_008596"] <- "Donor"


height = length(unique(SNPtable$Entry))
y <- vector()
x <- vector()
id <- vector()
values <- data.frame(ids = c("Pi", "phi", "brat", "gene"), value = c("Pi", "PHI", "BratNextGen", "Geneconv"))

for ( i in (1:nrow(Pi)) ) {
    id <- c(id, rep("Pi", each = 5))
    x <- c(x, Pi[i,1], Pi[i,2], Pi[i,2], Pi[i,1], Pi[i,1])
    y <- c(y, 0, 0, height, height, 0)
}
for ( i in (1:nrow(phi)) ) {
    id <- c(id, rep("phi", each = 5))
    x <- c(x, phi[i,1], phi[i,2], phi[i,2], phi[i,1], phi[i,1])
    y <- c(y, 0, 0, height, height, 0)
}
for ( i in (1:nrow(brat)) ) {
    id <- c(id, rep("brat", each = 5))
    x <- c(x, brat[i,1], brat[i,2], brat[i,2], brat[i,1], brat[i,1])
    y <- c(y, 0, 0, height, height, 0)
}
for ( i in (1:nrow(gene)) ) {
    id <- c(id, rep("gene", each =5))
    x <- c(x, gene[i,1], gene[i,2], gene[i,2], gene[i,1], gene[i,1])
    y <- c(y, 0, 0, height, height, 0)
}
positions <- data.frame(ids = id, X = x, Y = y)

datapoly <- merge(values, positions, by=c("ids"))

# graph of SNP table

p2 <- ggplot() + geom_point(aes(x=Pos, y=Entry, colour = SNPtable$Group), data=SNPtable)+ geom_polygon(aes(x=X, y=Y, group=ids, fill=value), data = datapoly, alpha = 0.3)  
pdf(file = "plot.pdf", width = 60, height = 5)
p2
dev.off()

