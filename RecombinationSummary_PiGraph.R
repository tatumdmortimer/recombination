#!/usr/bin/Rscript --vanilla 
library(ggplot2)
library(gridExtra)

#get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# load the recombination data
brat <- read.delim("brat.out")
gene <- read.delim("geneconv.out")
snp <- read.delim("SNPregions.out", header = FALSE)
RecipientFile <- paste("Pi_", args[1], "_Rec_convert.txt", sep = "")

Pi_R <- read.table(RecipientFile, sep = "\t")
colnames(Pi_R) <- c("Window","Pi")
Pi_R_Summary <- summary(Pi_R$Pi)
Pi_R_Title <- paste("Nucleotide Diversity: Recipient and ",args[1], sep = "")

plot1 <- ggplot(Pi_R, aes(x = Window, y = Pi)) + geom_line() + ggtitle(Pi_R_Title) + geom_ribbon(aes(ymin=Pi_R_Summary[2], ymax=Pi_R_Summary[5]), alpha = .2, fill = "green") + xlim(0,7500000) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(rep(0,4),"lines"))

y <- vector()
x <- vector()
id <- vector()
values <- data.frame(ids = c("snp", "brat", "gene"), value = c("SNP", "BratNextGen", "Geneconv"))

for ( i in (1:nrow(snp)) ) {
    id <- c(id, rep("snp", each = 5))
    x <- c(x, snp[i,1], snp[i,2], snp[i,2], snp[i,1], snp[i,1])
    y <- c(y, 0, 0, 1, 1, 0)
}
for ( i in (1:nrow(brat)) ) {
    id <- c(id, rep("brat", each = 5))
    x <- c(x, brat[i,1], brat[i,2], brat[i,2], brat[i,1], brat[i,1])
    y <- c(y, 1, 1, 2, 2, 1)
}
for ( i in (1:nrow(gene)) ) {
    id <- c(id, rep("gene", each = 5))
    x <- c(x, gene[i,1], gene[i,2], gene[i,2], gene[i,1], gene[i,1])
    y <- c(y, 2, 2, 3, 3, 2)
}
length(id)
length(x)
length(y)
gene[1,1]
positions <- data.frame(ids = id, X = x, Y = y)

datapoly <- merge(values, positions, by=c("ids"))

# graph of SNP table

plot2 <- ggplot() + geom_polygon(aes(x=X, y=Y, group=ids, fill=value), data = datapoly) + xlim(0,7500000) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.title = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), plot.margin = unit(rep(0,4),"lines"))
tiff(file = "plot.tiff", width = 15, height = 5, units = "in", res = 300, compression = c("jpeg"))
grid.arrange(plot1, plot2, nrow = 2, ncol = 1)
dev.off()

