library(ggplot2)
library(gridExtra)
library(reshape)
betweenStrep <- read.delim('betweenGenomeBias_strep_50rep.txt', header = FALSE, sep = " ")
betweenStaph <- read.delim('betweenGenomeBias_staph_50rep.txt', header = FALSE, sep = " ")
betweenEnt <- read.delim('betweenGenomeBias_ent_50rep.txt', header = FALSE, sep = " ")
betweenSmeg <- read.delim('betweenGenomeBias_smeg_50rep.txt', header = FALSE, sep = " ")
betweenCan <- read.delim('betweenGenomeBias_cani_50rep.txt', header = FALSE, sep = " ")

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

plotEnt <- ggplot(betweenEnt2, aes(x= value, y= (..density..), fill = variable, alpha = .9)) + 
  geom_density(position = "identity", binwidth = 1) + 
  ggtitle(expression(italic('E. faecium'))) +
  theme_bw(base_size = 8) + 
  xlim(-1,50) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) +
  ylab("Proportion of Pairwise Comparisons") +
  xlab("Overlapping Events")

plotSmeg <- ggplot(betweenSmeg2, aes(x= value, y= (..density..), fill = variable, alpha = .9)) + 
  geom_density(position = "identity", binwidth = 1) + 
  ggtitle(expression(italic('M. smegmatis'))) +
  theme_bw(base_size = 8) + 
  xlim(-1,50) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) +
  ylab("Proportion of Pairwise Comparisons") +
  xlab("Overlapping Events")

plotCan <- ggplot(betweenCan2, aes(x= value, y= (..density..), fill = variable, alpha = .9)) + 
  geom_density(position = "identity", binwidth = 1) + 
  ggtitle(expression(italic('M. canettii'))) +
  theme_bw(base_size = 8) + 
  xlim(-1,50) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) +
  scale_alpha(guide = 'none') +
  theme(legend.direction = 'horizontal') + 
  ylab("Proportion of Pairwise Comparisons") +
  xlab("Overlapping Events")

plotStaph <- ggplot(betweenStaph2, aes(x= value, y= (..density..), fill = variable, alpha = .9)) + 
  geom_density(position = "identity", binwidth = 1) + 
  ggtitle(expression(italic('S. aureus'))) +
  theme_bw(base_size = 8) + 
  xlim(-1,50) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) +
  ylab("Proportion of Pairwise Comparisons") +
  xlab("Overlapping Events")
plotStrep <- ggplot(betweenStrep2, aes(x= value, y= (..density..), fill = variable, alpha = .9)) + 
  geom_density(position = "identity", binwidth = 1) + 
  ggtitle(expression(italic('S. pneumoniae'))) +
  theme_bw(base_size = 8) + 
  xlim(-1,50) + ylim(0,1) +
  scale_fill_manual("", 
                    breaks = c("Observed", "Expected"),
                    values = c("black", "gray")) +
  ylab("Proportion of Pairwise Comparisons") +
  xlab("Overlapping Events")


legend <- g_legend(plotCan)

tiff(file = 'BetweenGenomeBias.tiff', width=6.5, height=4,units = "in", res = 300)
grid.arrange(arrangeGrob(plotCan + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines")), 
                         plotStrep + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines")), 
                         plotSmeg + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines")), 
                         plotEnt + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines")),
                         plotStaph + theme(legend.position = "none", plot.margin = unit(c(.25,.25,.25,.25), "lines")),
                         nrow = 2, ncol = 3))
dev.off()
