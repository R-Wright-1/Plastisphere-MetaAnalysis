library(metacoder)

args = commandArgs(trailingOnly=TRUE)

color1 = args[1]
color2 = args[2]

wd = getwd()
setwd <- (wd)
print(wd)

asvs <- read.csv("metacoder/metacoder.csv", header=TRUE, sep = ",")

obj <- parse_tax_data(asvs, class_cols = "lineage", class_sep = ";") 
obj$data$tax_data <- calc_obs_props(obj, "tax_data") 
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")
set.seed(1)
middle = '#C3C3C3'

cn = colnames(obj$data$tax_abund)
cn = cn[2:length(cn)]
treats = c()
for (a in 1:length(cn)) {
	str = cn[a]
	treats = c(treats, substr(str, 1, 6))
}

obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund", cols = cn, groups = treats, combinations=list(c('Treat1', 'Treat2')))
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method="fdr")
range(obj$data$diff_table$wilcox_p_value, finite=TRUE)
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0 
heat_tree(obj, node_label=taxon_names, node_size=n_obs, node_color=log2_median_ratio, node_color_interval=c(-3,3), node_color_range=c(color2, middle, color1), layout="davidson-harel", initial_layout="reingold-tilford", make_node_legend=FALSE, make_edge_legend=FALSE, output_file="metacoder/metacoder.pdf", node_label_max=0) 
heat_tree(obj, node_label=taxon_names, node_size=n_obs, node_color=log2_median_ratio, node_color_interval=c(-3,3), node_color_range=c(color2, middle, color1), layout="davidson-harel", initial_layout="reingold-tilford", make_node_legend=FALSE, make_edge_legend=FALSE, output_file="metacoder/metacoder_labels.pdf") 