library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(ggrepel)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Bathymodiolus_metagenomics/Thiodubiliella")

data <- read.table("Thiodubiliella.Freebayes.FINAL.CT.FORMAT", header=T, row.names=1)
com <- t(apply(data, 2, function(i) i/sum(i, na.rm=TRUE)))
data2 <- t(read.table("Thiodubiliella.Freebayes.FINAL.GT.FORMAT", header=T, row.names=1))
pop <- read.table("Thiodubiliella.clst", header=T)

dist <- vegdist(com, method="bray", binary=F, na.rm=TRUE)
matrix <- df2genind(data2, ploidy = 1, ncode = 1)
dist2 <- vegdist(matrix, method="euclidean", binary=F, na.rm=TRUE)
dist2[is.na(dist2)] = 0
pcoa <- pcoa(dist, correction = "cailliez")
pcoa2 <- pcoa(dist2, correction = "cailliez")

Axis1 <- pcoa$vectors[,1]
Axis2 <- pcoa$vectors[,2]
var1 <- round(pcoa$values[1,2], digits=4)*100
var2 <- round(pcoa$values[2,2], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
Axis1b <- pcoa2$vectors[,1]
Axis2b <- pcoa2$vectors[,2]
var1b <- round(pcoa2$values[1,2], digits=4)*100
var2b <- round(pcoa2$values[2,2], digits=4)*100
x_axisb <- paste("Axis1"," ","[",var1b,"%]",sep="")
y_axisb <- paste("Axis2"," ","[",var2b,"%]",sep="")

plot1 <- cbind(com, Axis1, Axis2, pop)
plot2 <- cbind(matrix, Axis1b, Axis2b, pop)

plot1$Vent <- factor(plot1$Vent, levels = c("Tow Cam", "Tahi Moana", "ABE 16", "ABE 09", "Tu'i Malila 16", "Tu'i Malila 09"))
plot2$Vent <- factor(plot2$Vent, levels = c("Tow Cam", "Tahi Moana", "ABE 16", "ABE 09", "Tu'i Malila 16", "Tu'i Malila 09"))

col <- c("Tow Cam" = "plum3", "Tahi Moana" = "royalblue3", "ABE 16" = "lightsteelblue4", "ABE 09" = "lightsteelblue", "Tu'i Malila 16" = "darkgoldenrod1", "Tu'i Malila 09" = "lightgoldenrod1")

# Colored by vent
pdf("Thiodubiliella_PCoA.12.pdf")
p <- ggplot(plot1, aes(Axis1, Axis2, color=Vent)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, name="Vent") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("Thiodubiliella_PCoA.geno.12.pdf")
p <- ggplot(plot2, aes(Axis1b, Axis2b, color=Vent)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, name="Vent") + labs(x = x_axisb, y = y_axisb)
p
dev.off()