library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Bathymodiolus_metagenomics/B_septemdierum/Transcriptome")

#IBS: Identical by descent/state, average genotype distance between individuals, ignores allele frequencies, used in MDS
data <- as.matrix(read.table("B_septemdierum_66.ibsMat"))
data[is.na(data)] = 0 
#COV: Covariance matrix, weighting scheme based on allele frequencies, used in PCA
data2 <- as.matrix(read.table("B_septemdierum_66.covMat"))
data2[is.na(data2)] = 0
pop <- read.table("B_septemdierum.clst", header=T)
pcoa <- pcoa(data, correction = "cailliez")
pca <- prcomp(data2, scale=TRUE, center=TRUE)

Axis1 <- pcoa$vectors.cor[,1]
Axis2 <- pcoa$vectors.cor[,2]
var1 <- round(pcoa$values[1,3], digits=4)*100
var2 <- round(pcoa$values[2,3], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
PC1 <- pca$rotation[,1]
PC2 <- pca$rotation[,2]
eigs <- pca$sdev^2
pc1_var <- round(eigs[1]/sum(eigs), digits=4)*100
pc2_var <- round(eigs[2]/sum(eigs), digits=4)*100
pc1_axis <- paste("PC1"," ","[",pc1_var,"%]",sep="")
pc2_axis <- paste("PC2"," ","[",pc2_var,"%]",sep="")
plot1 <- cbind(data, Axis1, Axis2, pop)
plot2 <- cbind(data2, PC1, PC2, pop)

plot1$Vent <- factor(plot1$Vent, levels = c("Tow Cam", "Tahi Moana", "ABE 16", "ABE 09", "Tu'i Malila 16", "Tu'i Malila 09"))
plot2$Vent <- factor(plot2$Vent, levels = c("Tow Cam", "Tahi Moana", "ABE 16", "ABE 09", "Tu'i Malila 16", "Tu'i Malila 09"))

col <- c("Tow Cam" = "plum3", "Tahi Moana" = "royalblue3", "ABE 16" = "lightsteelblue4", "ABE 09" = "lightsteelblue", "Tu'i Malila 16" = "darkgoldenrod1", "Tu'i Malila 09" = "lightgoldenrod1")

pdf("B_septemdierum_PCoA66.12.pdf")
p <- ggplot(plot1, aes(Axis1, Axis2, color=Vent)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, name="Vent") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("B_septemdierum_PCA66.12.pdf")
p <- ggplot(plot2, aes(PC1, PC2, color=Vent)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Vent)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, name="Vent") + labs(x = pc1_axis, y = pc2_axis)
p
dev.off()
