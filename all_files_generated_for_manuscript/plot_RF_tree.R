library("ggplot2")
library("phyloseq")
library("ape")
library(microbiome)

#Set the working directory to be the folder containing this script
wd = getwd()
setwd <- (wd)

#Import all of the data
asv_table <- read.csv("random_forest_R/random_forest.csv", sep=",", row.names=1)
asv_table = as.matrix(asv_table)
phy_tree <- read_tree("agglom_1/reduced_tree.tree")

#Convert these to phyloseq objects
ASV = otu_table(asv_table, taxa_are_rows = TRUE)
physeq = phyloseq(ASV,phy_tree)

#The following two commands will show you how many taxa there were before and after removing low abundance ASVs. You could choose a lower number for your remove low abundance, but the agglomeration may struggle if this is above ~10,000
physeq
png("random_forest_R/this_RF_tree.png", width=500, height=2500)
plot_tree(physeq, "treeonly", nodeplotblank)
#plot_tree(physeq, label.tips="taxa_names")
dev.off()
tree = phy_tree(physeq)
labels = c(tree$tip.label)
plot_tree(tree)
write.table(labels, "random_forest_R/this_RF_labels.csv", sep=",")
ape::write.tree(tree, paste(wd, "/random_forest_R/this_RF_tree.tree", sep=""))
