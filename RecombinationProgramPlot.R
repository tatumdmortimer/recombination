#!/usr/bin/Rscript --vanilla 
library(ggplot2)
library(gtable)
library(gridExtra)

#get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# load the recombination data
brat <- read.table(paste(args[1], "_writesegments_donor.txt", sep = ""), header = FALSE)
gene <- read.delim(paste(args[1], "_DonorRef_Group3Smeg_6.10.frags", sep = ""), header=F)
snp <- read.delim("SNPregions.out", header = FALSE)
cf <- read.delim(paste(args[1], "_cf.txt", sep = ""), header = FALSE)
co <- read.delim("../clonalOriginResults2.txt", header = TRUE)
co <- co[co$Position%%100 == 0,]

plot1 <- ggplot(cf, aes(x = V2, y = V3)) + 
    geom_point(color = "black", size = .5) + 
    geom_point(data = co, aes(x = Position, y = M20A/100), color = "gray", alpha = .9, size = .5) + 
    xlim(0,7000000) + 
    ylim(0,1) +
    theme_bw(base_size = 8) + 
    theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(vjust = 1.25, size = 6),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(.5,.5,-1,0),"lines")) + 
    ylab("ClonalFrame &\nClonalOrigin\nPosterior Probability")


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
positions <- data.frame(ids = id, X = x, Y = y)

datapoly <- merge(values, positions, by=c("ids"))

# graph of SNP table

plot2 <- ggplot() + 
    geom_polygon(aes(x=X, y=Y, group=ids), data = datapoly, fill = "black") + 
    xlim(0,7000000) + 
    theme_bw(base_size = 8) +
    theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(vjust = 2),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.margin = unit(c(0,.5,.5,0),"lines")) + 
    scale_y_discrete(limits=c("True Breakpoints", "BRATNextGen", "GENECONV")) 

g1 <- ggplot_gtable(ggplot_build(plot1))
g2 <- ggplot_gtable(ggplot_build(plot2))
maxWidth = unit.pmax(g1$widths[2:3], g2$widths[2:3])

g1$widths[2:3] <- maxWidth
g2$widths[2:3] <- maxWidth
tiff(file = "plot.tiff", width = 6.5, height = 2, units = "in", res = 300, compression = c("jpeg"))
grid.arrange(g1, g2, nrow = 2, ncol = 1)
dev.off()

