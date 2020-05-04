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

import all_functions_new as af
import os
import pandas as pd
import pickle
from datetime import datetime
start_time = datetime.now()

"""
Input files from this point are: 
- feature-table_w_tax.txt
- tree.nwk
- metadata.txt
- dna-sequences.fasta
"""

basedir = '/Users/robynwright/Documents/OneDrive/ACU_meta-analysis/Data_2020/paper_data_20-04-14/'
#os.chdir(basedir)
ft_tax, meta_fn, seqs, study_locs, study_dates = basedir+'qiime_output/feature-table_w_tax.txt', 'metadata.txt', basedir+'qiime_output/dna-sequences.fasta', 'Study_location.csv', 'Study_dates.csv'
n_jobs=10

"""
Set up folders
Add folders with all of the names that we will use to save either files from the analysis or figures
"""
#folder_names = ["agglom", "picrust", "figures", "ancom", "figures/ancom", "metacoder", "figures/metacoder", "random_forest", "random_forest/single_environment", "figures/random_forest", "figures/random_forest/single_environment", "random_forest_R", "random_forest/leave_one_dataset_out", "figues/random_forest/leave_one_dataset_out"]
#for fn in folder_names:
#    if os.path.exists(basedir+"/"+fn):
#       print(fn+" folder already exists\nSkipping making, but the files in this folder will be over-written if they have the same names as the ones saved by this script\n") 
#    else:
#        os.system("mkdir "+basedir+"/"+fn)

"""
FIRST PART
Format for R
"""
#ft, tax_dict = af.format_R(basedir+'/'+ft_tax, basedir)
tax_dict = af.open_pickle(basedir+'tax_dict.dictionary')

"""
SECOND PART
Run R phyloseq and unifrac and open these files afterwards
"""

#os.system("/usr/local/bin/Rscript agglomerate_and_unifrac.R")
ft_df, tree_agglom = af.open_and_sort(basedir+'/agglom/otu_table.csv'), basedir+'/agglom/reduced_tree.tree'
w_uf, uw_uf = af.open_and_sort(basedir+'agglom/weighted_unifrac.csv'), af.open_and_sort(basedir+'agglom/unweighted_unifrac.csv') #file names for unifrac distances based on agglomerated data
meta, meta_names, meta_dict = af.get_meta(meta_fn)
meta_df = af.get_meta_df(meta, meta_names, list(ft_df.columns))
ft_full = basedir+'feature_table.csv'
#ASV_dict = af.get_ASV_dict(ft_df, basedir+seqs)
#af.write_pickle(basedir+'ASV_dict.dictionary', ASV_dict)
ASV_dict = af.open_pickle(basedir+'ASV_dict.dictionary')

"""
THIRD PART
Format files for PICRUSt and run
"""

#seqs_agglom_fn, ft_agglom_fn = af.filter_seqs(ft_df, seqs) #filter the sequences file to have only the sequences that are in the agglomerated feature table
#os.system("source activate picrust2") #activate the picrust2 environment
#os.system("picrust2_pipeline.py -s "+seqs_agglom_fn+" -i "+ft_agglom_fn+" -o picrust/picrust_out --custom_trait_tables ko_all.txt --stratified --no_pathways") #run picrust with your files
#picrust_file = 'picrust/picrust_out/ko_all_metagenome_out/pred_metagenome_unstrat.tsv.gz' #picrust file name within your picrust folder
#ko_file = 'picrust/kegg_list.csv' #list of kegg orthologs including information on whether these are genes for xenobiotic degradation, pathogenesis or antimicrobial resistance
#KO_dict, KO_dict_full = af.make_KO_dict(ko_file) #make a dictionary with these kegg orthologs
#picrust = pd.read_csv(basedir+'/'+picrust_file, sep='\t', header=0, index_col=0) #now read in the picrust file
#picrust, KO_dict = af.filter_picrust(basedir+picrust_file, KO_dict, KO_dict_full) #and filter the picrust file based on the kegg orthologs that we want to keep

"""
FOURTH PART
Do all analyses and make all plots
Some parts of this will take a long time to run
Random Forests: (on a regular laptop, approximately overnight) and calls a lot of different functions within it, as well as an R script for plotting. 
Giving a lower value for est (estimators; default is 10000) will reduce running time as well as robustness
"""

#af.study_map(study_dates, study_locs, basedir) #get map and cumulative number of studies (we don't use the metadata file for this as it also contains data for studies for which no sequencing data were accessible/made available by the authors)
#af.study_metrics(meta_df, basedir) #get some basic information on number of samples in each environment and plastic types of all samples

#af.nmds_plot_study_env(w_uf, uw_uf, meta_dict, basedir) #get the nmds plot for weighted and unweighted unifrac distances for environment and study
#af.similarity_heatmap_combined(w_uf, uw_uf, basedir) #and get heatmaps showing average sample similarity between studies

#af.bar_dendro_venn(ft_df, ft_full, meta_dict, basedir, tax_dict) #get the dendrogram, stacked bar, heatmap, diversity and venn diagrams of overlapping ASVs

#[af.tree_heatmap(ft_df, meta_dict, basedir, tax_dict, level=al) for al in [1, 2, 3, 4, 5, 6]] #get the circular phylogenetic tree and heatmap for significant differences found by ancom. #level corresponds to phylogenetic level that we want to do comparisons on, with 0 being kingdom, 1 being phylum, etc. 7 is the lowest classification made for each. 
#af.metacoder(ft_df, tax_dict, meta_dict, basedir) #get the metacoder plots summarising differences between early and late incubation times, within and between plastic types

#af.get_random_forests(ft_df, tax_dict, meta_df, basedir, est=10000, n_jobs=n_jobs) #get all random forest models for metadata categories across all samples
#af.get_random_forest_plots(ft_df, tax_dict, ASV_dict, meta_dict, basedir) #make the summary heatmaps for each of these
#af.get_random_forests_leave_one_dataset_out(ft_df, tax_dict, meta_df, basedir, meta_dict, est=10000, n_jobs=n_jobs) #get all random forest models for metadata categories across all samples
#af.leave_one_out_plots(basedir, meta_df)

#af.get_environment_random_forest(ft_df, tax_dict, meta_df, meta_dict, basedir, est=10000, n_jobs=10) #get the separate random forests for only the general plastic type
#af.get_environment_random_forest_plots(ft_df, meta_df, tax_dict, ASV_dict, meta_dict, basedir) #and make the basic heatmaps
#af.make_env_rf_plot(ft_df, tax_dict, basedir, ASV_dict, meta_dict, mx=0.01)

#af.get_overall_random_forest_plot(ft_df, meta_df, tax_dict, ASV_dict, meta_dict, basedir)

#picrust = af.open_pickle(basedir+'picrust.dataframe') #and the picrust data
#KO_dict = af.open_pickle(basedir+'KO_dict.dictionary') #load the picrust kegg ortholog dictionary 
#af.picrust_plots(picrust, KO_dict, meta_dict, basedir) #now make the picrust volcano and bar plot

#af.plots_per_study(ft_df, meta_df, meta_dict, w_uf, uw_uf, tax_dict, ASV_dict, basedir, est=10000, n_jobs=10) #get the plots summarising the results per study

print('Running time was', datetime.now()-start_time)