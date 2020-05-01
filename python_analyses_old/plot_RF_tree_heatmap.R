library("ggplot2")
library("phyloseq")
library("ape")
#library(microbiome)

wd = getwd()
setwd <- (wd)

asv_table <- read.csv("random_forest_R/random_forest.csv", sep=",", row.names=1)
asv_table = as.matrix(asv_table)
phy_tree <- read_tree("agglom_1/reduced_tree.tree")
taxonomy <- read.csv("taxonomy_name_only.csv", sep=",", row.names=1)
taxonomy = as.matrix(taxonomy)
asv_df = as.data.frame(asv_table)
tax_df = as.data.frame(taxonomy, stringsAsFactors=F)

ASV = otu_table(asv_table, taxa_are_rows = TRUE)
physeq = phyloseq(ASV,phy_tree)
tree = phy_tree(physeq)

detach("package:phyloseq", unload=TRUE)
library("ggtree")

asv_tax_full <- merge(asv_df, tax_df, by="row.names")
asv_tax <- asv_tax_full[,c("Row.names", "Species.name")]

colnames(asv_tax) <- c("label", "Species.name")

labels = tree$tip.label
species_name = c()

for (i in 1:length(labels)) {
	for (j in 1:length(asv_tax$label)) {
		if (labels[i] == asv_tax$label[j]) {
			
			species_name <- c(species_name, asv_tax$Species.name[j])
		}
	}
}
d <- data.frame(label=tree$tip.label, species=species_name)

p <- ggtree(tree) %<+% d + xlim(NA, 6) + geom_tiplab(aes(label=species), parse=F, size=1, align=T, linesize=.12)
heat <- gheatmap(p, ASV, legend_title='Relative abundance', offset=1, colnames=TRUE, colnames_angle=90, font.size=1, hjust=0, colnames_position='top') + scale_fill_viridis_c()
pdf("random_forest_R/tree_and_heatmap.pdf")
heat
dev.off()