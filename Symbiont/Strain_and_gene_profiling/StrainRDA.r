library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(psych)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Bathymodiolus_metagenomics/Thiodubiliella/STRONG")

data <- t(read.table("mag_and_haplo_cov.tsv", header=T, row.names=1))

env <- read.table("Thiodubiliella.RDA.txt", header=T)

identical(rownames(data), env[,1]) 

pdf("Variable_correlation.pdf")
pairs.panels(env[,3:7], scale=T)
dev.off()

pred <- subset(env, select=-c(Latitude, Longitude))

rda <- capscale(data ~ Fluid + Depth + Year, data = pred, na.action = na.exclude, add = "cailliez")

RsquareAdj(rda)
vif.cca(rda) #Should be below 10

signif.full <- anova.cca(rda, parallel=getOption("mc.cores"))
signif.full
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
signif.margin <- anova.cca(rda, by="margin", parallel=getOption("mc.cores"))
signif.margin

env$Vent = factor(env$Vent, levels = c("Tow Cam", "Tahi Moana", "ABE 16", "ABE 09", "Tu'i Malila 16", "Tu'i Malila 09"))
col <- c("Tow Cam" = "plum3", "Tahi Moana" = "royalblue3", "ABE 16" = "lightsteelblue4", "ABE 09" = "lightsteelblue", "Tu'i Malila 16" = "darkgoldenrod1", "Tu'i Malila 09" = "lightgoldenrod1")

pdf("Thiodubiliella_RDA.pdf")
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)
points(rda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=col[env$Vent])
text(rda, scaling=3, display="bp", col="black", cex=1)
text(rda, scaling=3, display="species", col="black", cex=.5)
legend("bottomleft", legend=levels(env$Vent), bty="n", col="gray32", pch=21, cex=1, pt.bg=col)
dev.off()
