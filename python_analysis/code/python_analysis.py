#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:58:36 2020

@author: robynwright

"""

import all_functions as af
from Bio import SeqIO
import csv
from itertools import chain
import IPython
from lifelines.utils import concordance_index
import math
import matplotlib.colors as colors
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, TransformedBbox, BboxPatch, BboxConnector
import numpy as np
import os
import pandas as pd
from pdf2image import convert_from_path, convert_from_bytes
import pickle
from pip._internal.operations import freeze
from scipy.cluster import hierarchy
from scipy.spatial import distance
import scipy.spatial.distance as ssd
from sinfo import sinfo
from skbio.stats.composition import ancom
from sklearn import manifold, preprocessing
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split

''' Basic setup '''
basedir = '/Users/robynwright/Documents/OneDrive/Papers_writing/Plastisphere_Meta-analysis/test_recreate_analyses/paper_data/' 
meta_fn, seqs, seqs_rare, study_locs, study_dates, map_img = basedir+'metadata.txt', basedir+'qiime_output/dna-sequences_nr.fasta', basedir+'qiime_output/dna-sequences_rare.fasta', basedir+'Study_location.csv', basedir+'Study_dates.csv', basedir+'world_map.jpg'
ft_tax_rare = basedir+'qiime_output/feature-table_w_tax_rare.txt'
ft_tax_nr = basedir+'qiime_output/feature-table_w_tax_not_rare.txt'
ft_tax_nr_filt = basedir+'qiime_output/feature-table_w_tax_not_rare_filtered.txt'
n_jobs, est = 10, 10000

''' Make empty folders '''
folder_names = ["agglom", "picrust", "figures", "ancom", "separate_studies", "figures/ancom", "figures/metacoder", "random_forest", "random_forest/rare", "random_forest/log", "random_forest/clr", "random_forest/rel_abun", "random_forest/rare/single_environment", "random_forest/log/single_environment", "random_forest/clr/single_environment", "random_forest/rel_abun/single_environment", "figures/random_forest", "figures/random_forest/single_environment", "figures/random_forest/single_category"]
for fn in folder_names:
    if not os.path.exists(basedir+"/"+fn):
       os.system("mkdir "+basedir+"/"+fn)
       
''' Reformat QIIME2 files for R '''
if os.path.exists(basedir+'tax_dict_rare.dictionary'):
    tax_dict_rare = af.open_pickle(basedir+'tax_dict_rare.dictionary')
else:
    ft_rare, tax_dict_rare = af.format_R(ft_tax_rare, basedir, 'rare')

if os.path.exists(basedir+'tax_dict_nr.dictionary'):
    tax_dict_nr = af.open_pickle(basedir+'tax_dict_nr.dictionary')
else:
    ft_nr, tax_dict_nr = af.format_R(ft_tax_nr_filt, basedir, 'nr')

''' Perform agglomeration '''
'''
Call R script
'''

''' Make a taxonomy matrix for R '''
ft_nr = af.open_and_sort(basedir+'agglom/otu_table_nr.csv')
tax_matrix = pd.read_csv(basedir+'tax_dict_nr.csv', header=0, index_col=0)
tax_matrix = tax_matrix.loc[ft_nr.index.values, :]
tax_matrix.drop(['Species name'], axis=1, inplace=True)
tax_matrix = tax_matrix.reset_index()

''' Calculate distances '''
'''
Call R script
'''

''' Read the files into Python again and open the information that we need '''
# (Note here that R seems to sometimes save things with a different encoding than is expected by pandas. If you get errors related to NaNs here, just try opening the .csv files in excel and resaving them as .csv files (not actually changing anything, the default for excel is to save in utf-8, which is what pandas wants) should fix this. It otherwise reads in the dataframes as all NaNs sometimes. )
ft_df_rare, ft_df_log, ft_df_clr, ft_df_rel_abun, ft_df_nr = af.open_and_sort(basedir+'/agglom/otu_table_rare.csv'), af.open_and_sort(basedir+'/agglom/otu_table_log.csv'), af.open_and_sort(basedir+'/agglom/otu_table_clr.csv'), af.open_and_sort(basedir+'/agglom/otu_table_relabun.csv'), af.open_and_sort(basedir+'/agglom/otu_table_nr.csv')
tree_agglom_nr, tree_agglom_rare = basedir+'/agglom/reduced_tree_nr.tree', basedir+'/agglom/reduced_tree_rare.tree'
tax_dict_rare, tax_dict_nr = af.open_pickle(basedir+'tax_dict_rare.dictionary'), af.open_pickle(basedir+'tax_dict_nr.dictionary')
tax_matrix_nr, tax_matrix_rare = pd.read_csv(basedir+'tax_dict_nr.csv'), pd.read_csv(basedir+'tax_dict_rare.csv')

ft_full_nr, ft_full_rare = basedir+'/feature_table_nr.csv', basedir+'/feature_table_rare.csv'
meta, meta_names, meta_dict = af.get_meta(meta_fn)
meta_df = af.get_meta_df(meta, meta_names, list(ft_df_rare.columns))

''' Read in distance matrices '''
w_uf_rare, uw_uf_rare, w_uf_rel_abun, uw_uf_rel_abun, w_uf_log, uw_uf_log, philr_clr, eucl_dist = af.open_and_sort(basedir+'/agglom/weighted_unifrac_rare.csv'), af.open_and_sort(basedir+ '/agglom/unweighted_unifrac_rare.csv'), af.open_and_sort(basedir+'/agglom/weighted_unifrac_rel_abun.csv'), af.open_and_sort(basedir+ '/agglom/unweighted_unifrac_rel_abun.csv'), af.open_and_sort(basedir+'/agglom/weighted_unifrac_log.csv'), af.open_and_sort(basedir+ '/agglom/unweighted_unifrac_log.csv'), af.open_and_sort(basedir+'/agglom/philr_distance.csv'), af.get_eucl_dist(ft_df_clr.transpose())
eucl_dist_df = pd.DataFrame(eucl_dist, index=ft_df_clr.columns, columns=ft_df_clr.columns)
eucl_dist_df.to_csv(basedir+'agglom/euclidean_distance.csv')
eucl_clr = af.open_and_sort(basedir+'/agglom/euclidean_distance.csv')

w_uf_rare.to_csv(basedir+'agglom/sorted_weighted_unifrac_rare.csv'), uw_uf_rare.to_csv(basedir+'agglom/sorted_unweighted_unifrac_rare.csv'), w_uf_log.to_csv(basedir+'agglom/sorted_weighted_unifrac_log.csv'), uw_uf_log.to_csv(basedir+'agglom/sorted_unweighted_unifrac_log.csv'), w_uf_rel_abun.to_csv(basedir+'agglom/sorted_weighted_unifrac_rel_abun.csv'), uw_uf_rel_abun.to_csv(basedir+'agglom/sorted_unweighted_unifrac_rel_abun.csv'), philr_clr.to_csv(basedir+'agglom/sorted_philr_distance.csv'), eucl_clr.to_csv(basedir+'agglom/sorted_euclidean_distance.csv')

if os.path.exists(basedir+'ASV_dict_nr.dictionary'):
    ASV_dict_nr = af.open_pickle(basedir+'ASV_dict_nr.dictionary')
else:
    ASV_dict_nr = af.get_ASV_dict(ft_df_nr, seqs, basedir)
    af.write_pickle(basedir+'ASV_dict_nr.dictionary', ASV_dict_nr)
if os.path.exists(basedir+'ASV_dict_rare.dictionary'):
    ASV_dict_rare = af.open_pickle(basedir+'ASV_dict_rare.dictionary')
else:
    ASV_dict_rare = af.get_ASV_dict(ft_df_rare, seqs_rare, basedir)
    af.write_pickle(basedir+'ASV_dict_rare.dictionary', ASV_dict_rare)

''' Get Figure 1 (study map and metrics) '''
if not os.path.exists(basedir+'/figures/Fig1_map_metrics'+af.ext):
  af.study_map_metrics(study_dates, study_locs, basedir, map_img, meta_df)

''' Get Figure S1 (richness and number of reads for all studies) '''
if os.path.exists(basedir+'feature_table_raw.csv'):
    ft_df_raw = af.open_and_sort(basedir+'feature_table_raw.csv')
else:
    ft_df_raw, tax_dict_raw = af.format_R(ft_tax_nr_filt, basedir, 'raw')

if not os.path.exists(basedir+'figures/FigS1_richness'+af.ext):
  af.calculate_richness(ft_df_rare, ft_df_nr, ft_df_raw, meta_df, basedir)

''' Get Figure 2 (NMDS plots) '''
# First get the groupings:
if os.path.exists(basedir+'anosim_permanova_all.list'):
  stats = af.open_pickle(basedir+'anosim_permanova_all.list')
  run_tests = False
else:
  group_env, group_study = af.get_stats_groups(w_uf_rare, meta_dict)
  run_tests = True

# Carry out the stats tests
'''
Call R script
'''

# Make the NMDS plots
'''
if run_tests:
  stats_rare = r.rare_stats
  stats_rel_abun = r.rel_abun_stats
  stats_log = r.log_stats
  stats_clr = r.clr_stats
  stats = [stats_rare, stats_rel_abun, stats_log, stats_clr]
  write_pickle(basedir+'anosim_permanova_all.list', stats)
'''

count = 0
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  if not os.path.exists(basedir+'/figures/Fig2_nmds_overall_'+norm+af.ext):
    this_stats = [[stats[count][:4], stats[count][4:8]], [stats[count][8:12], stats[count][12:]]]
    if norm != 'clr':
      af.nmds_plot_study_env('agglom/sorted_weighted_unifrac_'+norm+'.csv', 'agglom/sorted_unweighted_unifrac_'+norm+'.csv', meta_dict, basedir, stats=this_stats, norm_name=norm)
    else:
      af.nmds_plot_study_env('agglom/sorted_euclidean_distance.csv', 'agglom/sorted_philr_distance.csv', meta_dict, basedir, stats=this_stats, norm_name=norm)
  count += 1

''' Get Figure 3 (distance heatmaps) '''
dist_mats = [[w_uf_rare, uw_uf_rare], [w_uf_rel_abun, uw_uf_rel_abun], [w_uf_log, uw_uf_log], [eucl_clr, philr_clr]]
count = 0
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  if not os.path.exists(basedir+'/figures/Fig3_heatmap_combined_'+norm+af.ext):
    af.similarity_heatmap_combined(dist_mats[count][0], dist_mats[count][1], basedir, norm_name=norm) 
  count += 1

''' Get Figure 4 (summarising sample types, groupings and the number of taxa shared within environments and sample types) '''
count = 0
ft_full_nr = basedir+'/feature_table_nr.csv'
ft_full_rare = basedir+'/feature_table_rare.csv'
full_fts = [ft_full_rare, ft_full_nr, ft_full_nr, ft_full_nr]
norm_fts = [ft_df_rare, ft_df_rel_abun, ft_df_log, ft_df_clr]
tax_dicts = [tax_dict_rare, tax_dict_nr, tax_dict_nr, tax_dict_nr]
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  if not os.path.exists(basedir+'/figures/Fig4_dendro_venn_'+norm+af.ext):
    af.bar_dendro_venn(norm_fts[count], full_fts[count], meta_dict, basedir, tax_dicts[count], str_norm=norm)
  count += 1

''' Get the overall random forests for each normalisation methods '''
af.get_random_forests(ft_df_rare, tax_dict_rare, meta_df, basedir, est=est, n_jobs=n_jobs, norm='rare')
af.get_random_forests(ft_df_log, tax_dict_nr, meta_df, basedir, est=est, n_jobs=n_jobs, norm='log')
af.get_random_forests(ft_df_clr, tax_dict_nr, meta_df, basedir, est=est, n_jobs=n_jobs, norm='clr')
af.get_random_forests(ft_df_rel_abun, tax_dict_nr, meta_df, basedir, est=est, n_jobs=n_jobs, norm='rel_abun')

''' Get the random forests for each environment for general plastic type '''
af.get_environment_random_forest(ft_df_rare, tax_dict_rare, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='rare') 
af.get_environment_random_forest(ft_df_log, tax_dict_nr, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='log') 
af.get_environment_random_forest(ft_df_clr, tax_dict_nr, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='clr') 
af.get_environment_random_forest(ft_df_rel_abun, tax_dict_nr, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='rel_abun') 

''' Get Figure S2 (comparison of different normalisation methods on overall random forest accuracy and concordance) '''
if not os.path.exists(basedir+'figures/FigS2_RF_compare'+af.ext):
  af.get_rf_comparison(basedir)

''' Get Figure S3 (comparison of different normalisation methods on random forest accuracy and concordance within environments for general plastic type) '''
if not os.path.exists(basedir+'/figures/FigS3_RF_env_compare'+af.ext):
  af.get_rf_comparison_env(basedir)

''' Get Figure 5, Figure S4 and figures for Supplementary Sections 2 and 3 (the overall random forest plots) '''
fts = [ft_df_rare, ft_df_rel_abun, ft_df_log, ft_df_clr]
tax_dicts = [tax_dict_rare, tax_dict_nr, tax_dict_nr, tax_dict_nr]
ASV_dicts = [ASV_dict_rare, ASV_dict_nr, ASV_dict_nr, ASV_dict_nr]
count = 0
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  af.get_overall_random_forest_plot(fts[count], meta_df, tax_dicts[count], ASV_dicts[count], meta_dict, basedir, norm=norm)
  count += 1
  
af.get_random_forest_plots(ft_df_rare, tax_dict_rare, ASV_dict_rare, meta_dict, basedir, norm='rare')
af.get_random_forest_plots(ft_df_clr, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='clr')
af.get_random_forest_plots(ft_df_log, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='log')
af.get_random_forest_plots(ft_df_rel_abun, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='rel_abun')

os.system('cp '+basedir+'/figures/random_forest/ASVs_overall_forest_rare.png '+basedir+'/figures/FigS4_rare_ASV_accuracy.png')

''' Get Figure S5 and more for Supplementary Section 3 (random forest plots for general plastic type within each environment) '''
af.get_environment_random_forest_plots(ft_df_rare, meta_df, tax_dict_rare, ASV_dict_rare, meta_dict, basedir, norm='rare', skip_individual=False)
af.get_environment_random_forest_plots(ft_df_log, meta_df, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='log', skip_individual=False)
af.get_environment_random_forest_plots(ft_df_rel_abun, meta_df, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='rel_abun', skip_individual=False)
af.get_environment_random_forest_plots(ft_df_clr, meta_df, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='clr', skip_individual=False)
if not os.path.exists(basedir+'figures/FigS5_environment_random_forest_0.005_rare.pdf'):
  if not os.path.exists('figures/FigS5_environment_random_forest_0.005_rare'+af.ext):
    af.make_env_rf_plot(ft_df_rare, tax_dict_rare, basedir, ASV_dict_rare, meta_dict, norm='rare')
if not os.path.exists(basedir+'figures/FigS5_environment_random_forest_0.01_rare.pdf'):
  if not os.path.exists('figures/FigS5_environment_random_forest_0.01_rare'+af.ext):
    af.make_env_rf_plot(ft_df_rare, tax_dict_rare, basedir, ASV_dict_rare, meta_dict, mx=0.01, norm='rare')

''' Get Figure 6 and figures for Supplementary Section 4 (metacoder plots) '''
count = 0
ft_dfs = [ft_df_rare, ft_df_rel_abun, ft_df_log, ft_df_clr]
tax_dicts = [tax_dict_rare, tax_dict_nr, tax_dict_nr, tax_dict_nr]
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  af.metacoder_py(ft_dfs[count], tax_dicts[count], meta_dict, basedir, norm=norm)
  ft_dfs[count].drop(['lineage'], axis=1, inplace=True)
  count += 1

# Convert these PDFs to .png (or whatever extension we are using)
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  path_to_pdf = basedir+'figures/metacoder/'+norm+'/'
  af.convert_pdf(path_to_pdf)
  
# Make a colorbar
af.make_colorbar_fig(basedir)

''' Get figures for Supplementary Section 5 (ANCOM trees and heatmaps) '''
[af.tree_heatmap(ft_df_rare, meta_dict, basedir, tax_dict_rare, level=al) for al in [1, 2, 3, 4, 5, 6]]

# Convert these PDFs to .png
path_to_pdf = basedir+'figures/ancom/'
af.convert_pdf(path_to_pdf)

''' Get figures for Supplementary Section 1 (separate plots for each study) '''
if af.fs_main == 10:
  af.fs_main -= 2
  af.fs_small -= 2
  af.fs_title -= 2
af.plots_per_study(ft_df_rare, meta_df, meta_dict, w_uf_rare, uw_uf_rare, tax_dict_rare, ASV_dict_rare, basedir, est=est, n_jobs=n_jobs) #get the plots summarising the results per study
