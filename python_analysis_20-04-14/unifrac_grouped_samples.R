library("phyloseq")

wd = getwd()
setwd <- (wd)

#Import all of the data (these are large files so this may take a while)
asv_table <- read.csv("grouped_samples_for_unifrac.csv", sep=",", row.names=1)
asv_table = as.matrix(asv_table)
phy_tree <- read_tree("tree.nwk")

#Convert these to phyloseq objects
ASV = otu_table(asv_table, taxa_are_rows = TRUE)

#You can just run the following commands to ensure that all samples and OTUs have the same names in each matrix - you will see that phy_tree also contains other names (these are from the reference database), but don't worry about this - it will sort itself out when you combine the objects
#sample_names(ASV)[1:100]
#sample_names(META)[1:100]

#Now we will add together the OTU data and the tree data and remove those ASVs that are below 1% relative abundance
physeq = phyloseq(ASV,phy_tree)

w_uf <- UniFrac(physeq, weighted=TRUE, normalized=FALSE, fast=TRUE)

#And finally we can write the data to file so that we can continue our analysis in Python
w_uf_df <- as.data.frame(as.matrix(w_uf))
write.csv(w_uf_df, paste(wd, "/weighted_unifrac.csv", sep=""))
