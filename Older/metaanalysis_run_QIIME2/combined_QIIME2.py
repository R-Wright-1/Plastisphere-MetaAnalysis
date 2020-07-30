#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:30:03 2020

@author: robynwright
"""

import os
from datetime import datetime
from optparse import OptionParser
parser = OptionParser()

start_time = datetime.now()

parser.add_option("--cores", dest="cores", default=1, help="Numeric value of how many cores to use. Default=1")
parser.add_option("--max-depth", dest="max_depth", default=995391, help="Numeric value of the maximum sequencing depth. Default=995391")
parser.add_option("--rare-depth", dest="rare_depth", default=2000, help="Numeric value of the depth to rarefy to and remove samples with below this number of sequences. Default=2000")
(options, args) = parser.parse_args()

cores, max_depth, rare_depth = options.cores, options.max_depth, options.rare_depth
if cores == 1 and max_depth == 995391 and rare_depth == 2000:
    print("Using all default values as you haven't changed any")
    print("This might take a long time to run with only one core")

if not os.path.exists("merged_table.qza"):
    print("merged_table.qza doesn't exist")
    print("You haven't merged your files, or renamed the original files")
    exit()

str1 = "qiime feature-table summarize \
		     --i-table merged_table.qza  \
		     --o-visualization merged_table_summary.qzv"
str2 = "qiime feature-classifier classify-sklearn \
		      --i-reads merged_representative_sequences.qza \
		      --i-classifier ref_alignments/classifier_silva_132_99_16S.qza \
		      --p-n-jobs "+str(cores)+" \
		      --output-dir taxa"
str3 = "qiime tools export \
		      --input-path taxa/classification.qza \
		      --output-path taxa"
str4 = "qiime feature-table filter-features \
		      --i-table merged_table.qza \
		      --p-min-frequency 10 \
		      --p-min-samples 1 \
		      --o-filtered-table merged_table_filtered.qza"
str5 = "qiime taxa filter-table \
		      --i-table merged_table_filt.qza \
		      --i-taxonomy taxa/classification.qza \
		      --p-include D_1__ \
		      --exclude mitochondria, chloroplast \
		      --o-filtered-table merged_table_filtered_contamination.qza"
str6 = "qiime feature-table summarize \
		      --i-table merged_table_filtered_contamination.qza \
		      --o-visualization merged_table_filtered_contamination_summary.qzv"
str7 = "qiime diversity alpha-rarefaction \
		      --i-table merged_table_filtered_contamination.qza \
		      --p-max-depth "+str(max_depth)+" \
	 	      --p-steps 20 \
		      --p-metrics 'observed_otus' \
		      --o-visualization merged_rarefaction_curves.qzv"
str8 = "qiime feature-table filter-samples \
		      --i-table merged_table_filtered_contamination.qza \
		      --p-min-frequency "+str(rare_depth)+" \
		      --o-filtered-table  merged_table_final.qza"
str9 = "qiime feature-table rarefy \
		       --i-table merged_table_final.qza \
		       --p-sampling-depth "+str(rare_depth)+" \
		       --o-rarefied-table merged_table_final_rarefied.qza"
str10 = "qiime feature-table filter-seqs \
		      --i-data representative_sequences.qza \
		      --i-table merged_table_final_rarefied.qza \
		      --o-filtered-data  representative_sequences_final_rarefied.qza"
str11 = "qiime tools export \
		      --input-path representative_sequences_final_rarefied.qza \
		      --output-path exports"
str12 = "sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv"
str13 = "qiime tools export \
		      --input-path merged_table_final_rarefied.qza \
		      --output-path exports"
str14 = "biom add-metadata \
		      --i exports/feature-table.biom \
		      --o exports/feature-table_w_tax.biom \
		      --observation-metadata-fp taxa/taxonomy.tsv \
		      --sc-separated taxonomy"
str15 = "biom convert \
		      --i exports/feature-table_w_tax.biom \
		      --o exports/feature-table_w_tax.txt \
		      --to-tsv \
		      --header-key taxonomy"
str16 = "qiime fragment-insertion sepp \
		      --i-representative-sequences representative_sequences_final_rarefied.qza \
		      --i-reference-database ref_alignments/sepp-refs-silva-128.qza \
		      --o-tree insertion_tree_rarefied.qza \
		      --placements insertion_placements_rarefied.qza \
		      --p-threads "+str(cores)
str17 = "qiime diversity core-metrics-phylogenetic \
		      --i-table merged_table_final_rarefied.qza \
		      --i-phylogeny insertion_tree_rarefied.qza \
		      --p-sampling-depth 2000 \
		      --m-metadata-file metadata.txt \
		      --p-n-jobs "+str(cores)+" \
		      --output-dir diversity"
str18 = "qiime tools export \
		       --input-path insertion_tree_rarefied.qza \
		       --output-path exports"
str19 = "qiime tools export \
		       --input-path diversity/weighted_unifrac_distance_matrix.qza \
		       --output-path diversity"
str20 = "mv diversity/distance-matrix.tsv exports/weighted_unifrac_not_agglom.tsv"
str21 = "qiime tools export \
		       --input-path diversity/unweighted_unifrac_distance_matrix.qza \
		       --output-path diversity"
str22 = "mv diversity/distance-matrix.tsv exports/unweighted_unifrac_not_agglom.tsv"

os.system(str1)
os.system(str2)
os.system(str3)
os.system(str4)
os.system(str5)
os.system(str6)
os.system(str7)
os.system(str8)
os.system(str9)
os.system(str10)
os.system(str11)
os.system(str12)
os.system(str13)
os.system(str14)
os.system(str15)
os.system(str16)
os.system(str17)
os.system(str18)
os.system(str19)
os.system(str20)
os.system(str21)
print('Finished running, running time was: ', datetime.now() - start_time)