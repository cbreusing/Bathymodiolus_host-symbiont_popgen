library(ggplot2)
library(vegan)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(phangorn)
library(GUniFrac)
library(ape)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Bathymodiolus_metagenomics/Thiodubiliella")

#data <- read.table("Thiodubiliella_fst.txt", header=T)
#newmat <- xtabs(fst ~ sample1 + sample2, data=data)
#write.table(newmat,"Thiodubiliella_fst.mat", sep="\t", col.names = TRUE, row.names = TRUE)

#data <- t(read.table("STRONG/mag_and_haplo_cov.tsv", header=T, row.names=1))
#tree <- read.tree(file="STRONG/haplotypes_tree.nwk")
#rooted_tree <- midpoint(tree)
#beta <- GUniFrac(data, rooted_tree)
#dw <- beta$unifracs[, , "d_1"]

data <- as.matrix(read.table("Thiodubiliella_fst.matrix", header=T, row.names=1))
pop <- read.table("Thiodubiliella.clst3", header=T)
colRamp = colorRamp2(c(-1,0,1), c("darkcyan", "oldlace", "darkslateblue"))
sample_colors <- list(Vent = c("Tow Cam" = "plum3", "Tahi Moana" = "royalblue3", "ABE" = "#8F9FB4", "Tu'i Malila" = "#FFD24D"))

ha <- HeatmapAnnotation(df = pop, annotation_height = unit(4, "mm"), col = list(Vent = c("Tow Cam" = "plum3", "Tahi Moana" = "royalblue3", "ABE" = "#8F9FB4", "Tu'i Malila" = "#FFD24D")), show_legend=FALSE)
hb <- rowAnnotation(df = pop, annotation_height = unit(4, "mm"), col = list(Vent = c("Tow Cam" = "plum3", "Tahi Moana" = "royalblue3", "ABE" = "#8F9FB4", "Tu'i Malila" = "#FFD24D")), show_legend=FALSE)

pdf("Fst-Pst_between_samples.pdf")
map <- Heatmap(data, top_annotation = ha, left_annotation = hb, km = 1, name = "Fst+Pst", col = colRamp, show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_column_dend = FALSE, column_dend_reorder = FALSE, heatmap_legend_param = list(legend_direction = "horizontal"))
draw(map, heatmap_legend_side = "top")
dev.off()
