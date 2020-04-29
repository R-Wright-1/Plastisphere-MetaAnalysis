#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:58:36 2020

@author: robynwright

Before running, you should check that:
- QIIME2 is installed (this uses version 2019.10)
- PICRUSt2 is installed (this uses version XX) and the environment is called picrust2
- R is installed
- You have all documents from: github_address
- You have checked all packages from README.rtf are installed
"""

import all_functions as af
import os
import pickle

"""
Input files from this point are: 
- feature-table_w_tax.txt
- tree.nwk
- metadata.txt
- dna-sequences.fasta
"""

basedir = os.getcwd()
ft_tax = 'qiime_output/feature-table_w_tax.txt'
meta_fn = 'metadata.txt'
seqs = 'qiime_output/dna-sequences.fasta'
study_locs, study_dates = 'Study_location.csv', 'Study_dates.csv'

"""
Set up folders
"""
folder_names = ["agglom_1", "agglom_2", "agglom_05", "figures", "ancom", "figures/ancom", "metacoder", "figures/metacoder", "random_forest", "figures/random_forest", "random_forest_R"]
for fn in folder_names:
    if os.path.exists(basedir+"/"+fn):
       print(fn+" folder already exists. Skipping making, but the files in this folder will be over-written if they have the same names as the ones saved by this script.") 
    else:
        os.system("mkdir "+basedir+"/"+fn)

"""
FIRST PART
Format for R
"""
#ft, tax_dict = af.format_R(ft_tax)
meta, meta_names, meta_dict = af.get_meta(meta_fn)

"""
SECOND PART
Run R phyloseq and unifrac
"""

#os.system("/usr/local/bin/Rscript agglomerate_and_unifrac.R")
ft_agglom = os.getcwd()+'/agglom_1/otu_table.csv'
tree_agglom = os.getcwd()+'/agglom_1/reduced_tree.tree'

"""
THIRD PART
Format and run PICRUSt?
"""

#seqs_agglom_fn, ft_agglom_fn, ft_agglom = af.filter_seqs(ft_agglom, seqs)
#os.system("source activate picrust2")
#os.system("picrust2_pipeline.py -s "+seqs_agglom_fn+" -i "+ft_agglom_fn+" -o picrust/picrust_out --custom_trait_tables ko_all.txt --stratified --no_pathways")

"""
FOURTH PART
Make all plots
Some parts of this will take a long time to run
Random Forests: (on a regular laptop, approximately overnight) and calls a lot of different functions within it, as well as an R script for plotting. Giving a lower value for est (default is 10000) will reduce running time as well as robustness
"""

w_uf, uw_uf = 'agglom_1/weighted_unifrac.csv', 'agglom_1/unweighted_unifrac.csv'
w_uf_not_agglom, uw_uf_not_agglom = 'qiime_output/weighted_unifrac_not_agglom.tsv', 'qiime_output/unweighted_unifrac_not_agglom.tsv'
#af.study_map(study_dates, study_locs, basedir) #get map and cumulative number of studies (we don't use the metadata file for this as it also contains data for studies for which no sequencing data were accessible/made available by the authors)
#af.nmds_plot(w_uf, meta_dict, basedir) #get nmds plot for agglomerated weighted unifrac similarity matrix
#af.nmds_plot(uw_uf, meta_dict, basedir) #get nmds plot for agglomerated unweighted unifrac similarity matrix
#af.nmds_plot(w_uf_not_agglom, meta_dict, basedir) #get nmds plot for not agglomerated weighted unifrac similarity matrix
#af.nmds_plot(uw_uf_not_agglom, meta_dict, basedir) #get nmds plot for not agglomerated unweighted unifrac similarity matrix
#af.similarity_heatmap(w_uf, basedir) #get similarity heatmap plot for agglomerated weighted unifrac similarity matrix
#af.similarity_heatmap(uw_uf, basedir) #get similarity heatmap plot for agglomerated unweighted unifrac similarity matrix
#af.similarity_heatmap(w_uf_not_agglom, basedir) #get similarity heatmap plot for not agglomerated weighted unifrac similarity matrix
#af.similarity_heatmap(uw_uf_not_agglom, basedir) #get similarity heatmap plot for not agglomerated unweighted unifrac similarity matrix
#af.random_forests(ft_agglom, tax_dict, basedir) #get details and plots for which ASVs are contributing to differences between which meta-categories
with open('tax_dict.dictionary', 'rb') as tax_dict:
    tax_dict = pickle.load(tax_dict)
#af.bar_dendro_venn(ft_agglom, meta_dict, basedir, tax_dict) #get the dendrogram, stacked bar, heatmap, diversity and venn diagrams of overlapping ASVs
#for al in [1, 2, 3, 4, 5, 7]: #level corresponds to phylogenetic level that we want to do comparisons on, with 0 being kingdom, 1 being phylum, etc. 7 is the lowest classification made for each. We don't include 6 here because this is species only, while 7 is the lowest taxonomic classification
#    af.tree_heatmap(ft_agglom, meta_dict, basedir, tax_dict, level=al) #get the circular phylogenetic tree and heatmap for significant differences found by ancom 
af.metacoder(ft_agglom, tax_dict, meta_dict, basedir)
"""
FIFTH PART
Run Metacoder in R
"""