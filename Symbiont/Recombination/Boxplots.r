library(ggplot2)
library(cowplot)
library(patchwork)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Bathymodiolus_metagenomics/Thiodubiliella")

data <- read.table("Thiodubiliella_pi.txt", header=T)

p <- ggplot(data, aes(x=site, y=intra_pi, fill=site)) + geom_boxplot(outlier.shape=NA) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) + theme_classic() + scale_fill_manual(values=c("#8F9FB4", "royalblue3", "plum3", "#FFD24D")) + labs(x="Vent", y="Intra-host nucleotide diversity") + scale_x_discrete(limits=c("Tu'i Malila", "ABE", "Tahi Moana", "Tow Cam")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("Boxplot_pi.pdf")
p
dev.off()

shapiro.test(data$intra_pi)
pairwise.t.test(data$intra_pi, data$site, p.adjust.method = "fdr", pool.sd = FALSE)


data <- read.table("rhometa_local_norm.txt", header=T)

p1 <- ggplot(data, aes(x=Species, y=RhoTheta, fill=Species)) + geom_boxplot(outlier.shape=NA) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) + geom_point(aes(y=gRhoTheta), size=7, shape=16, color="black") + theme_classic() + scale_fill_manual(values=c("darkcyan", "dodgerblue", "red3", "slateblue")) + labs(x="Symbiont", y="Rho/Theta") + scale_x_discrete(limits=c("Ca. Thiodubiliella", "Ifr-SOX", "Gamma1", "Epsilon")) + theme(legend.position = "none", text = element_text(size = 17), axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- ggplot(data, aes(x=Species, y=rm, fill=Species)) + geom_boxplot(outlier.shape=NA) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) + geom_point(aes(y=grm), size=7, shape=16, color="black") + theme_classic() + scale_fill_manual(values=c("darkcyan", "dodgerblue", "red3", "slateblue")) + labs(x="Symbiont", y="r/m") + scale_x_discrete(limits=c("Ca. Thiodubiliella", "Ifr-SOX", "Gamma1", "Epsilon")) + theme(legend.position = "none", text = element_text(size = 17), axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Boxplot_rhometa_norm.pdf", width=18, height=10)
plot_grid(p1, p2, labels = c('A', 'B'), align="hv", label_size = 30)
dev.off()

shapiro.test(data$RhoTheta)
pairwise.t.test(data$RhoTheta, data$Species, p.adjust.method = "fdr", pool.sd = FALSE)

shapiro.test(data$rm)
pairwise.t.test(data$rm, data$Species, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(data$rm, data$Species, p.adjust.method = "fdr", pool.sd = FALSE)
