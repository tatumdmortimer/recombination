#!/usr/bin/Rscript --vanilla

library(ggplot2)
library(gridExtra)
library(reshape)

betweenStrep <- read.delim('betweenGenomeBias_strep_50rep_prop.txt', header = FALSE, sep = " ")
betweenStaph <- read.delim('betweenGenomeBias_staph_50rep_prop.txt', header = FALSE, sep = " ")
betweenEnt <- read.delim('betweenGenomeBias_ent_50rep_prop.txt', header = FALSE, sep = " ")
betweenSmeg <- read.delim('betweenGenomeBias_smeg_50rep_prop.txt', header = FALSE, sep = " ")
betweenCan <- read.delim('betweenGenomeBias_can_50rep_prop.txt', header = FALSE, sep = " ")

betweenEnt <- data.frame(t(betweenEnt))
betweenSmeg <- data.frame(t(betweenSmeg))
betweenCan <- data.frame(t(betweenCan))
betweenStaph <- data.frame(t(betweenStaph))
betweenStrep <- data.frame(t(betweenStrep))

colnames(betweenEnt) <- c("Observed", "Expected")
colnames(betweenSmeg) <- c("Observed", "Expected")
colnames(betweenCan) <- c("Observed", "Expected")
colnames(betweenStrep) <- c("Observed", "Expected")
colnames(betweenStaph) <- c("Observed", "Expected")

betweenEnt2 <- melt(betweenEnt, measure=c("Observed", "Expected"))
betweenSmeg2 <- melt(betweenSmeg, measure=c("Observed", "Expected"))
betweenCan2 <- melt(betweenCan, measure=c("Observed", "Expected"))
betweenStrep2 <- melt(betweenStrep, measure=c("Observed", "Expected"))
betweenStaph2 <- melt(betweenStaph, measure=c("Observed", "Expected"))

plotEnt <- ggplot(betweenEnt2, aes(x= value*2, fill = variable, alpha = .9)) + 
  geom_histogram(aes(y=..density..*.05),position = "identity", binwidth = .05) + 
  ggtitle(expression(italic('E. faecium'))) +
  theme_bw(base_size = 8) + 
  xlim(0,1) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) 

plotSmeg <- ggplot(betweenSmeg2, aes(x= value*2, fill = variable, alpha = .9)) + 
  geom_histogram(aes(y=..density..*.05),position = "identity", binwidth = .05) + 
  ggtitle(expression(italic('M. smegmatis'))) +
  theme_bw(base_size = 8) + 
  xlim(0,1) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) 

plotCan <- ggplot(betweenCan2, aes(x= value*2, fill = variable, alpha = .9)) + 
  geom_histogram(aes(y=..density..*.05),position = "identity", binwidth = .05) + 
  ggtitle(expression(italic('M. canettii'))) +
  theme_bw(base_size = 8) + 
  xlim(0,1) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) +
  scale_alpha(guide = 'none') +
  theme(legend.direction = 'horizontal')  

plotStaph <- ggplot(betweenStaph2, aes(x= value*2, fill = variable, alpha = .9)) + 
  geom_histogram(aes(y=..density..*.05),position = "identity", binwidth = .05) + 
  ggtitle(expression(italic('S. aureus'))) +
  theme_bw(base_size = 8) + 
  xlim(0,1) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) 
plotStrep <- ggplot(betweenStrep2, aes(x= value*2, fill = variable, alpha = .9)) + 
  geom_histogram(aes(y=..density..*.05),position = "identity", binwidth = .05) + 
  ggtitle(expression(italic('S. pneumoniae'))) +
  theme_bw(base_size = 8) + 
  xlim(0,1) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) 



tiff(file = 'BetweenGenomeBias2.tiff', width=6.5, height=4,units = "in", res = 300)
grid.arrange(arrangeGrob(plotCan + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines"), axis.title = element_blank()), 
                         plotStrep + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines"), axis.title = element_blank()), 
                         plotStaph + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines"), axis.title = element_blank()), 
                         plotSmeg + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines"), axis.title = element_blank()),
                         plotEnt + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines"), axis.title = element_blank()),
                         nrow = 2, ncol = 3), 
                         sub = textGrob("       Proportion of Recombination Shared Between Pairs", gp = gpar(fontsize = 10, fontface = "bold", vjust = -.35)),
                         left = textGrob("Proportion of Pairwise Comparisons", gp = gpar(fontsize = 10, fontface = "bold"), rot = 90, vjust = .75))
dev.off()
