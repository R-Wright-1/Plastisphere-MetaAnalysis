library("ggplot2")
library("phyloseq")
library("ape")
library(ggnewscale)

wd = getwd()
setwd <- (wd)

asv_table <- read.csv("ancom/ancom_significant.csv", sep=",", row.names=1, check.names=FALSE)
asv_table = as.matrix(asv_table)
colors = asv_table['Colors', ]
asv_table <- asv_table[!rownames(asv_table) %in% c('Colors'), ]
class(asv_table) <- "numeric"

phy_tree <- read_tree("agglom_1/reduced_tree.tree")
taxonomy <- read.csv("ancom/taxonomy_name_only.csv", sep=",", row.names=1)
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

p <- ggtree(tree, layout="fan", open.angle=90)
offset = 0
for (n in 1:length(colnames(ASV))) {
	name = colnames(ASV)[n]
	this_df = ASV[, name]
	cols = colors[[name]]
	cols = strsplit(cols, ", ")
	p <- p + new_scale_fill()
	p <- gheatmap(p, this_df, offset=offset, width=.15, colnames_angle=90, font.size=2, hjust=1, color="black") + scale_fill_gradient2(low=cols[[1]][1], mid="white", high=cols[[1]][2])
	offset = offset+0.4
}
offset = offset+0.5
p <- p %<+% d + geom_tiplab(aes(label=species), parse=F, size=2, align=T, linesize=.12, offset=offset)
p <- p + theme(legend.position="none")
pdf("ancom/tree_and_heatmap.pdf")
p
dev.off()