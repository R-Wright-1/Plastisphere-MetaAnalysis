# Food or just a free ride? A meta-analysis reveals the global diversity of the Plastisphere


This document contains all of the commands necessary for the analyses performed in:
**Food or just a free ride? A meta-analysis reveals the global diversity of the Plastisphere**.

Please contact [Robyn Wright](mailto:robyn.wright@dal.ca) with any questions.

You can see the version of this file that contains example figures at [this FigShare file](https://doi.org/10.6084/m9.figshare.12923855).

## 1. Studies included

This provides links to all studies included in this meta-analysis, listed by the names that I have used to refer to them throughout. 

### Studies
<li>[AmaralZettler2015](doi.org/10.1890/150017)</li>
<li>[AriasAndres2018](doi.org/10.1016/j.envpol.2018.02.058)</li>
<li>[Canada2020](doi.org/10.1016/j.aquaculture.2019.734540)</li>
<li>[Curren2019](doi.org/10.1016/j.scitotenv.2018.11.250)</li>
<li>[Delacuvellerie2019](doi.org/10.1016/j.jhazmat.2019.120899)</li>
<li>[DeTender2015](doi.org/10.1021/acs.est.5b01093)</li>
<li>[DeTender2017](doi.org/10.1021/acs.est.7b00697)</li>
<li>[DussudHudec2018](doi.org/10.3389/fmicb.2018.01571)</li>
<li>[DussudMeistertzheim2018](doi.org/10.1016/j.envpol.2017.12.027)</li>
<li>[ErniCassola2019](doi.org/10.1007/s00248-019-01424-5)</li>
<li>[Esan2019](doi.org/10.1371/journal.pone.0214376)</li>
<li>[Frere2018](doi.org/10.1016/j.envpol.2018.07.023)</li>
<li>[Hoellein2014](doi.org/10.1371/journal.pone.0098485)</li>
<li>[Hoellein2017](doi.org/10.1086/693012)</li>
<li>[Jiang2018](doi.org/10.1016/j.scitotenv.2017.12.105)</li>
<li>[Kesy2019](doi.org/10.3389/fmicb.2019.01665)</li>
<li>[Kirstein2018](doi.org/10.1016/j.marenvres.2018.09.028)</li>
<li>[Kirstein2019](doi.org/10.1371/journal.pone.0215859)</li>
<li>[McCormick2014](doi.org/10.1021/es503610r)</li>
<li>[McCormick2016](doi.org/10.1002/ecs2.1556)</li>
<li>[Oberbeckmann2016](doi.org/10.1371/journal.pone.0159289)</li>
<li>[Oberbeckmann2018](doi.org/10.3389/fmicb.2017.02709)</li>
<li>[Ogonowski2018](doi.org/10.1111/1462-2920.14120)</li>
<li>[Parrish2019](doi.org/10.1039/c8ew00712h)</li>
<li>[Pinto2019](doi.org/10.1371/journal.pone.0217165)</li>
<li>[Pollet2018](doi.org/10.1093/femsec/fiy083)</li>
<li>[Rosato2020](doi.org/10.1016/j.scitotenv.2019.135790)</li>
<li>[Syranidou2019](doi.org/10.1016/j.jhazmat.2019.04.078)</li>
<li>[SyranidouPE2017](doi.org/10.1371/journal.pone.0183984)</li>
<li>[SyranidouPS2017](doi.org/10.1038/s41598-017-18366-y)</li>
<li>[Tagg2019](doi.org/10.1016/j.marpolbul.2019.06.013)</li>
<li>[Woodall2018](doi.org/10.1371/journal.pone.0206220)</li>
<li>[Wu2019](doi.org/10.1016/j.watres.2019.114979)</li>
<li>[Xu2019](doi.org/10.1016/j.marpolbul.2019.05.036)</li>
<li>[Zhang2019](doi.org/10.1016/j.scitotenv.2019.06.108)</li>
<br/>
If you are adding a study then you should add details of this (following the existing formats) to the files in the 'python_analysis_20-04-14' folder:<br/>
<ul>
<li>Study_dates.csv</li>
<li>Study_location.csv</li>
<li>metadata.txt</li>
</ul>

## 2. Setup environment 

I have made this notebook using R studio and Python and used python in a conda environment named r-reticulate. You can see the session information below for the versions of these and of packages that I use.

If you do not have any of these packages already installed then you will need to install them. 

### Setup reticulate

If you have problems with R finding Python then it might be worth explicitly telling R where to find the Python version you want to use, as described [here](https://rstudio.github.io/reticulate/reference/use_python.html) or [here](https://rstudio.github.io/reticulate/articles/versions.html). What worked for me was changing it in my Rprofile - this is supposedly in the folder given by running R.home() in the R studio console:
```{R, cache=TRUE}
R.home()
```

However, when I looked here, I couldn't find any Rprofile and found (from searching in finder) that the only Rprofile was actually in:
```{}
/Library/Frameworks/R.framework/Versions/3.6/Resources/library/base/R/Rprofile
```

I then opened this in a text editor and added this line to the end of the file:
```{}
Sys.setenv(RETICULATE_PYTHON = "/Users/robynwright/opt/miniconda3/envs/r-reticulate/bin/python")
```

### Setup R

```{R, setup_r, results='hide', message=FALSE, warning=FALSE}
load_all_R_packages <- function() {
  library(reticulate)
  library(kableExtra)
  library(knitr)
  library(exactRankTests)
  library(nlme)
  library(dplyr)
  library(ggplot2)
  library(vegan)
  library(phyloseq)
  library(ape)
  library(microbiome)
  library(ggnewscale)
  library(metacoder)
  library(philr)
  library(coda.base)
  library(ALDEx2)
  library(compositions)
}
load_all_R_packages()
```

### Setup Python

I have found conda-forge to be most successful for installing the majority of these packages.

```{python, setup_python, results='hide', message=FALSE}
import os
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
from datetime import datetime
import importlib
from itertools import chain
import IPython
from lifelines.utils import concordance_index
import math
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, TransformedBbox, BboxPatch, BboxConnector
import numpy as np
from optparse import OptionParser
import os
import pandas as pd
from pdf2image import convert_from_path, convert_from_bytes
import pickle
from pip._internal.operations import freeze
import random
from scipy.cluster import hierarchy
from scipy.stats import pearsonr
from scipy.spatial import distance
import scipy.spatial.distance as ssd
import scipy.stats as stats
from sinfo import sinfo
from skbio import DistanceMatrix
from skbio.stats.distance import anosim
from skbio.stats.distance import permanova
from skbio.stats.composition import ancom
from sklearn.cluster import AgglomerativeClustering
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
```

### Print session information 

#### R
```{R, r_session_info}
sessionInfo()
```

#### Python
Python doesn't really have a great option for session info like R does, but here I have printed first all of the packages that I use, and then all of the packages that are installed within my environment.
```{python, python_session_info}
sinfo()
print('\n\n')
x = freeze.freeze()
for p in x:
    print(p)
```

## 3. Getting data 

I have provided all of the data for this study in several different ways:
<ul>
<li>[The raw reads files in QIIME2 .qza format](link)</li>
<li>[QIIME2 files separately for each study](https://doi.org/10.6084/m9.figshare.12217682) - these include the deblur outputs for each study as well as the QIIME2 representative sequences and feature tables.</li>
<li>[QIIME2 merged files and output](https://doi.org/10.6084/m9.figshare.12227522) - including the merged QIIME2 representative sequences and feature table, as well as the taxonomic classifications, filtered tables, SEPP insertion tree, etc.</li>
</ul>
If you want to either add a study to this or recreate the analyses used here, then you can skip ahead to section 3 and use the merged_representative_sequences.qza and merged_table.qza files for this. 

Otherwise, details are given below on how to download the reads for all studies for which data is in the NCBI SRA. Data was provided directly by the authors for Curren2019, Kirstein2018, Kirstein2019 (this is available via NCBI for both Kirstein studies, but I asked for the data without primers removed) and Pollet2018. Frere2018 is available from: VAMPS portal (www.vamps. mbl.edu) under the project name LQM_MPLA_Bv4v5.

### Downloading NCBI reads

Here I use AmaralZettler2015 with accession SRP026054 as an example.

#### (I) Get accession list
e.g. download from [this link](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP026054&o=acc_s%3Aa)

#### (II) Use this to download the reads 
For this you will need the [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit):
```{bash, eval=FALSE}
prefetch --option-file SRR_Acc_List.txt
```
These files by default get added to an ncbi/sra/ folder in your home directory and can then be moved wherever you like.

#### (III) Convert this to fastq R1 and R2 files:
```{bash, eval=FALSE}
for i in folder_name/* ; do fastq-dump -split-files --gzip $i ; done
```

#### (IV) I then renamed the files using a Python script. 

If you are not adding many files then this is probably easier to just do manually, but get in touch if you want details of how I did this. I renamed all files to follow this format:
AZ2015s001_S001_L001_R1_001.fastq.gz	AZ2015s001_S001_L001_R2_001.fastq.gz
AZ2015s002_S002_L001_R1_001.fastq.gz	AZ2015s002_S002_L001_R2_001.fastq.gz
*e.g.* these are the R1 and R2 files for samples 1 and 2 from AmaralZettler2015 - I don't have details of the run anything was performed on or the lane samples were run in, so I have just edited the file names to follow this format so as they are easily compatible with QIIME2. 

#### (V) And I added all of the first parts of these names (e.g. AZ2015s001, AZ2015s002) as the names for these samples in my metadata file. 
**The most current version of this is [here](https://github.com/R-Wright-1/Plastisphere-MetaAnalysis/blob/master/python_analysis_20-04-14/metadata.txt). As long as your sample names follow this format (i.e. sample name before the first underscore) then the subsequent parts of this analysis shouldn't struggle even if your naming is different.**

## 4. QIIME2 analysis 

I used 12 threads on a server for most of my analyses, but you can change this in these code chunks accordingly if you have less available. I think some parts will probably struggle with so many samples if you try to do this locally. 
This follows the Microbiome helper tutorial [here](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2020.2)).
You can view all of the summary files (when you use the summarise commands) [here](https://view.qiime2.org/).

Activate the QIIME2 environment (if you do not already have this installed then follow the instructions [here](https://docs.qiime2.org/2020.6/install/):
```{bash, eval=FALSE}
conda activate qiime2-2019.10
```

### A. For each individual study

Note that this uses [Deblur](https://github.com/biocore/deblur). [DADA2](https://benjjneb.github.io/dada2/) could also be used, but given that we didn't know whether all samples from each study came from the same sequencing run, we chose the per sample denoising approach of Deblur. 

#### (I) Run quality checks on the data:
```{bash, eval=FALSE}
mkdir fastqc_out
fastqc -t 12 raw_data/*.fastq.gz -o fastqc_out
multiqc fastqc_out/
```
You should now look at the multiqc_report.html to ensure that everything is high enough quality to proceed.

#### (II) Import your data to the QIIME2 format:
```{bash, eval=FALSE}
qiime tools import \
            --type SampleData[PairedEndSequencesWithQuality] \
            --input-path raw_data/ \
            --output-path reads.qza \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
```

#### (III) Trim primers (if present) with cutadapt. 
The primer sequences shown here are for 341F and 802R - these will need changing if you have used different primers:
```{bash, eval=FALSE}
qiime cutadapt trim-paired \
            --i-demultiplexed-sequences reads.qza \
            --p-cores 12 \
            --p-front-f CCTACGGGNGGCWGCAG \
            --p-front-r GACTACHVGGGTATCTAATCC \
            --p-discard-untrimmed \
            --p-no-indels \
            --o-trimmed-sequences reads_trimmed.qza
```

#### (IV) Summarize the trimmed files:
```{bash, eval=FALSE}
qiime demux summarize \
            --i-data reads_trimmed.qza \
            --o-visualization reads_trimmed_summary.qzv
```

#### (V) Join paired ends 
If the reads were already trimmed then just use reads.qza as the input here:
```{bash, eval=FALSE}
qiime vsearch join-pairs \
            --i-demultiplexed-seqs reads_trimmed.qza \
            --o-joined-sequences reads_joined.qza
```

#### (VI) Summarize the joined pairs 
If too many reads were removed then you may need to play around with some of the other options at [here](https://docs.qiime2.org/2020.2/plugins/available/vsearch/join-pairs/):
```{bash, eval=FALSE}
qiime demux summarize \
            --i-data reads_joined.qza \
            --o-visualization reads_joined_summary.qzv
```

#### (VII) Filter out low quality reads:
```{bash, eval=FALSE}
qiime quality-filter q-score-joined \
            --i-demux reads_joined.qza \
            --o-filter-stats filt_stats.qza \
            --o-filtered-sequences reads_joined_filtered.qza
```

#### (VIII) Summarize these reads (and look at where to trim) 
You should look at the positions where the quality starts to drop below 30 and use these as trim lengths:
```{bash, eval=FALSE}
qiime demux summarize \
            --i-data reads_joined_filtered.qza \
            --o-visualization reads_joined_filtered_summary.qzv
```

#### (IX) Run deblur 
You can remove the --p-left-trim-len if you don't need to remove any from this end:
```{bash, eval=FALSE}
qiime deblur denoise-16S \
            --i-demultiplexed-seqs reads_joined_filtered.qza \
            --p-trim-length 402 \
            --p-left-trim-len 0 \
            --p-sample-stats \
            --p-jobs-to-start 12 \
            --p-min-reads 1 \
            --output-dir deblur_output_quality
```

#### (X) Summarize the feature table to see how many reads we now have:
```{bash, eval=FALSE}
qiime feature-table summarize \
            --i-table deblur_output_quality/table.qza  \
            --o-visualization deblur_table_summary.qzv
```

### B. Merge studies

#### (I) Rename the current representative sequences and merged tables:
```{bash, eval=FALSE}
mv merged_representative_sequences.qza previous_merged_representative_sequences.qza
mv merged_table.qza previous_merged_table.qza
```

#### (II) Combine feature tables. 
You will need to replace 'your_folder_name' with the folder that contains your tables to be added (when I did this for all studies, I just added additional --i-tables table_name.qza lines):
```{bash, eval=FALSE}
qiime feature-table merge \
            --i-tables your_folder_name/deblur_output_quality/table.qza \
            --i-tables previous_merged_table.qza \
            --o-merged-table merged_table.qza
```

#### (III) Combine the sequence objects. 
You will again need to replace 'your_folder_name' with the folder that contains your sequences to be added (when I did this for all studies, I just added additional --i-data representative_sequences_name.qza lines):
```{bash, eval=FALSE}
qiime feature-table merge-seqs \
            --i-data your_folder_name/deblur_output/representative_sequences.qza \
            --i-data previous_merged_representative_sequences.qza \
            --o-merged-data merged_representative_sequences.qza
```

### C. Combined processing

Now that all of the samples that we are looking at are combined into the merged sequences and table files, we can classify and analyze them.

#### (I) Summarize the combined feature tables (this is to check that everything looks OK after the merges, and can be skipped if not necessary):
```{bash, eval=FALSE}
qiime feature-table summarize \
            --i-table merged_table.qza  \
            --o-visualization merged_table_summary.qzv
```

#### (II) Classify the features (this part will probably take the longest - it may take at least a day or so and is the part that may not be possible on a local computer):
```{bash, eval=FALSE}
qiime feature-classifier classify-sklearn \
            --i-reads merged_representative_sequences.qza \
            --i-classifier ref_alignments/classifier_silva_132_99_16S.qza \
            --p-n-jobs 12 \
            --output-dir taxa
```
As these sequences come from different 16S regions, I downloaded the full length 16S classifier from [here](https://docs.qiime2.org/2020.6/data-resources/). There is now an updated SILVA version, but I used the Silva 132 classifier (this can only improve upon classification accuracy, so I recommend using the latest one).

#### (III) Export this file to look at the classifications:
```{bash, eval=FALSE}
qiime tools export \
            --input-path taxa/classification.qza \
            --output-path taxa
```

#### (IV) Filter low abundance features:
```{bash, eval=FALSE}
qiime feature-table filter-features \
            --i-table merged_table.qza \
            --p-min-frequency 10 \
            --p-min-samples 1 \
            --o-filtered-table merged_table_filtered.qza
```

#### (V) Filter potential contaminants and those not classified at the kingdom level:
```{bash, eval=FALSE}
qiime taxa filter-table \
            --i-table merged_table_filtered.qza \
            --i-taxonomy taxa/classification.qza \
            --p-include D_1__ \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table merged_table_filtered_contamination.qza
```

#### (VI) Summarize the filtered table:
```{bash, eval=FALSE}
qiime feature-table summarize \
            --i-table merged_table_filtered_contamination.qza \
            --o-visualization merged_table_filtered_contamination_summary.qzv
```
Now find out how many features you have as well as the maximum sample depth (this is the "Maximum Frequency" in the "Frequency per sample" section).

#### (VII) Obtain rarefaction curves for samples:
```{bash, eval=FALSE}
qiime diversity alpha-rarefaction \
            --i-table merged_table_filtered_contamination.qza \
            --p-max-depth 995391 \
            --p-steps 20 \
            --p-metrics 'observed_otus' \
            --o-visualization merged_rarefaction_curves.qzv
```

#### (VIII) Filter samples that have below 2000 reads:
```{bash, eval=FALSE}
qiime feature-table filter-samples \
            --i-table merged_table_filtered_contamination.qza \
            --p-min-frequency 2000 \
            --o-filtered-table  merged_table_final.qza
```

#### (IX) Rarefy remaining samples to 2000:
```{bash, eval=FALSE}
qiime feature-table rarefy \
            --i-table merged_table_final.qza \
            --p-sampling-depth 2000 \
            --o-rarefied-table merged_table_final_rarefied.qza
```

#### (X) Filter the sequences to contain only those that are in the rarefied feature table:
```{bash, eval=FALSE}
qiime feature-table filter-seqs \
            --i-data merged_representative_sequences.qza \
            --i-table merged_table_final_rarefied.qza \
            --o-filtered-data  representative_sequences_final_rarefied.qza
```

#### (XI) Export feature table and sequences for rarefied data:
```{bash, eval=FALSE}
qiime tools export \
            --input-path representative_sequences_final_rarefied.qza \
            --output-path exports
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv
qiime tools export \
            --input-path merged_table_final_rarefied.qza \
            --output-path exports
biom add-metadata \
            -i exports/feature-table.biom \
            -o exports/feature-table_w_tax.biom \
            --observation-metadata-fp taxa/taxonomy.tsv \
            --sc-separated taxonomy
biom convert \
            -i exports/feature-table_w_tax.biom \
            -o exports/feature-table_w_tax.txt \
            --to-tsv \
            --header-key taxonomy
```

#### (XII) Obtain phylogenetic tree using SEPP fragment insertion and the silva reference database for rarefied data:
```{bash, eval=FALSE}
qiime fragment-insertion sepp \
            --i-representative-sequences representative_sequences_final_rarefied.qza \
            --i-reference-database ref_alignments/sepp-refs-silva-128.qza \
            --o-tree insertion_tree_rarefied.qza \
            --o-placements insertion_placements_rarefied.qza \
            --p-threads 12
```
You can download the reference file [here](https://docs.qiime2.org/2020.6/data-resources/). At the time of writing, this still used Silva 128, but I would recommend using an updated version if there is one.

#### (XIII) Export the resulting insertion tree:
```{bash, eval=FALSE}
qiime tools export \
            --input-path insertion_tree_rarefied.qza \
            --output-path exports
```

#### (XIV) The files inside the exports folder should then be copied to the folder that the subsequent analyses will be carried out in, e.g.:
```{bash, eval=FALSE}
for i in exports/* ; cp $i paper_data/qiime_output/; done
```

And rename the files:
```{bash, eval=FALSE}
mv paper_data/qiime_output/feature-table_w_tax.txt paper_data/qiime_output/feature-table_w_tax_rare.txt
mv paper_data/qiime_output/dna-sequences.fasta paper_data/qiime_output/dna-sequences_rare.fasta
mv paper_data/qiime_output/tree.nwk paper_data/qiime_output/tree_rare.nwk
```

#### Now do the same for the non-rarefied data

#### (XV) Move the data to a new folder and summarize the number of features etc:
```{bash, eval=FALSE}
mkdir not_rarefied
cp merged_table_final.qza not_rarefied/
cp representative_sequences_final.qza not_rarefied/
cd not_rarefied
qiime feature-table summarize \
            --i-table merged_table_final.qza \
            --o-visualization merged_table_final_summary.qzv
```
Now go to [QIIME2 view](https://view.qiime2.org/) and look at this merged_table_final_summary.qzv
This tells me that I have 2,056 samples, 162,656 features and a median of 21,654 reads per sample. So I am dividing this median number by 100 and using this as the min-frequency here:
```{bash, eval=FALSE}
qiime feature-table filter-features \
            --i-table merged_table_final.qza \
            --p-min-frequency 217 \
            --p-min-samples 1 \
            --o-filtered-table merged_table_final_filtered.qza

qiime feature-table summarize \
            --i-table merged_table_final_filtered.qza \
            --o-visualization merged_table_final_filtered_summary.qzv

qiime feature-table filter-seqs \
            --i-data representative_sequences_final.qza \
            --i-table merged_table_final_filtered.qza \
            --o-filtered-data  representative_sequences_final_filtered.qza
```
For the 2056 samples I include, this now gives a median of 18,962 and 24,873 features (much faster to make a tree with)

#### (XVI) Export feature table and sequences for not rarefied data (before and after filtering, because we still want to know how many sequences in each sample to calculate relative abundance):
```{bash, eval=FALSE}
qiime tools export \
           --input-path representative_sequences_final_filtered.qza \
           --output-path exports
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxonomy.tsv
qiime tools export \
           --input-path merged_table_final_filtered.qza \
           --output-path exports
biom add-metadata \
           -i exports/feature-table.biom \
           -o exports/feature-table_w_tax_filtered.biom \
           --observation-metadata-fp taxonomy.tsv \
           --sc-separated taxonomy
biom convert \
           -i exports/feature-table_w_tax_filtered.biom \
           -o exports/feature-table_w_tax_filtered.txt \
           --to-tsv \
           --header-key taxonomy

qiime tools export \
           --input-path merged_table_final.qza \
           --output-path exports
biom add-metadata \
           -i exports/feature-table.biom \
           -o exports/feature-table_w_tax.biom \
           --observation-metadata-fp taxonomy.tsv \
           --sc-separated taxonomy
biom convert \
           -i exports/feature-table_w_tax.biom \
           -o exports/feature-table_w_tax.txt \
           --to-tsv \
           --header-key taxonomy
```

#### (XVI) Now put these sequences into the tree:
```{bash, eval=FALSE}
qiime fragment-insertion sepp \
            --i-representative-sequences representative_sequences_final_filtered.qza \
            --i-reference-database /home/robyn/other/ref_alignments/sepp-refs-silva-128.qza \
            --o-tree insertion_tree_not_norm.qza \
            --o-placements insertion_placements_not_norm.qza \
            --p-threads 24
```

#### (XVII) Export the tree:
```{bash, eval=FALSE}
qiime tools export \
           --input-path insertion_tree_not_norm.qza \
           --output-path exports
```

#### (XVIII) The files inside the exports folder should then be copied to the folder that the subsequent analyses will be carried out in, e.g.:
```{bash, eval=FALSE}
for i in exports/* ; cp $i paper_data/qiime_output/; done
```

And rename the files:
```{bash, eval=FALSE}
mv paper_data/qiime_output/feature-table_w_tax.txt paper_data/qiime_output/feature-table_w_tax_not_rare.txt
mv paper_data/qiime_output/feature-table_w_tax_filtered.txt paper_data/qiime_output/feature-table_w_tax_not_rare_filtered.txt
mv paper_data/qiime_output/dna-sequences.fasta paper_data/qiime_output/dna-sequences_not_rare.fasta
mv paper_data/qiime_output/tree.nwk paper_data/qiime_output/tree_not_rare.nwk
```

## 5. Data and statistical analysis 

### A. Introduction

This can be used to recreate all of the analyses found in the Plastisphere meta-analysis paper, or alternatively to re-do these analyses while also incorporating additional data.

### B. Other files needed

These can be found [here](https://github.com/R-Wright-1/Plastisphere-MetaAnalysis/tree/master/files_needed):
<ul>
<li>metadata.txt</li>
<li>Study_dates.csv</li>
<li>Study_location.csv</li>
<li>world_map.jpg</li>
</ul>

This goes through all analyses run, but you can download all of the analysis files used and created from [this Figshare file](https://doi.org/10.6084/m9.figshare.12227303). In particular, taking the random_forest and ancom files will make a very big difference to how long this takes to run. If you just want to re-make the figures (possibly with changes) then you can add all files but remove the figures folder.

It is expected that the base directory contains, at a minimum:
<ul>
<li>metadata.txt</li>
<li>Study_dates.csv</li>
<li>Study_location.csv</li>
<li>world_map.jpg</li>
<li>qiime_output
<ul>
<li>dna-sequences_rare.fasta</li>
<li>dna-sequences_nr.fasta</li>
<li>feature-table_w_tax_not_rare_filtered.txt</li>
<li>feature-table_w_tax_not_rare.txt</li>
<li>feature-table_w_tax_rare.txt</li>
<li>tree_nr.nwk</li>
<li>tree_rare.nwk</li>
</ul>
</li>

### C. Analysis

The majority of these code chunks contain statements to ensure that they aren't run if the output already exists, but if you want to run them anyway then you should change these/move the files already created somewhere else. Many steps will rely on the output of previous steps, although I have tried to ensure that once run they will save their output so that they won't need to be re-run multiple times (they are often very time/computationally intensive). 
The sections with an asterisk will need running each time you want to run anything (and after the initial run should only take a few seconds to run), but the others are only if you are running the analyses for that section.

#### (*I) Set the base directory for all inputs and outputs and the basic files and define all functions that we will call later on.

Python:
```{python, give_input, results='hide', fig.keep='all', message=FALSE}
basedir = '/Users/robynwright/Documents/OneDrive/Papers_writing/Plastisphere_Meta-analysis/test_recreate_analyses/paper_data/' 
meta_fn, seqs, seqs_rare, study_locs, study_dates, map_img = basedir+'metadata.txt', basedir+'qiime_output/dna-sequences_nr.fasta', basedir+'qiime_output/dna-sequences_rare.fasta', basedir+'Study_location.csv', basedir+'Study_dates.csv', basedir+'world_map.jpg'
ft_tax_rare = basedir+'qiime_output/feature-table_w_tax_rare.txt'
ft_tax_nr = basedir+'qiime_output/feature-table_w_tax_not_rare.txt'
ft_tax_nr_filt = basedir+'qiime_output/feature-table_w_tax_not_rare_filtered.txt'
n_jobs, est = 10, 10000

seed = 3 #random seed
fs_title, fs_main, fs_small = 14, 10, 8 #font sizes to use in figures
color_env = {'marine':'#03A498', 'freshwater':'#03B1FC', 'aquatic':'#9A75FC', 'terrestrial':'#FD9A64'} #colors for plotting each environment
label_env = ['Marine', 'Freshwater', 'Aquatic', 'Terrestrial'] #labels for each environment
ext, dpi = '.png', 600 #extension and dots per inch to save figures  with
color_source = {'aliphatic':'#5F9EA0', 'biofilm':'#FCA94A', 'other plastic':'#8B008B', 'unknown plastic':'#3593FC', 'planktonic':'#F9E79F', 'blank':'gray', 'early':'#FBC704', 'late':'#900C3F', 'collection':'gray'}
phylo_level_names = {0:'kingdoms', 1:'phyla', 2:'classes', 3:'orders', 4:'families', 5:'genera', 6:'species', 7:'ASVs'}
rename_plots = {'AmaralZettler':'Amaral Zettler $et$ $al$. 2015', 'AriasAndres':'Arias-Andres $et$ $al$. 2018', 'Canada':'Canada $et$ $al$. 2020', 'Curren':'Curren & Leong 2019', 'Delacuvellerie':'Delacuvellerie $et$ $al$. 2019', 'DeTender':'De Tender $et$ $al$. 2015', 'DeTenderB':'De Tender $et$ $al$. 2017', 'DussudHudec':'Dussud, Hudec $et$ $al$. 2018', 'DussudMeistertzheim':'Dussud, Meistertzheim $et$ $al$. 2018', 'ErniCassola':'Erni-Cassola $et$ $al$. 2019', 'Esan':'Esan $et$ $al$. 2019', 'Frere':'Frere $et$ $al$. 2018', 'Hoellein':'Hoellein $et$ $al$. 2014', 'HoelleinB':'Hoellein $et$ $al$. 2017', 'Jiang':'Jiang $et$ $al$. 2018', 'Kesy':'Kesy $et$ $al$. 2019', 'Kirstein':'Kirstein $et$ $al$. 2018', 'KirsteinB':'Kirstein $et$ $al$. 2019', 'McCormick':'McCormick $et$ $al$. 2014', 'McCormickB':'McCormick $et$ $al$. 2016', 'Oberbeckmann':'Oberbeckmann $et$ $al$. 2016', 'OberbeckmannB':'Oberbeckmann $et$ $al$. 2018', 'Ogonowski':'Ogonowski $et$ $al$. 2018', 'Parrish':'Parrish $et$ $al$. 2019', 'Pinto':'Pinto $et$ $al$. 2019', 'Pollet':'Pollet $et$ $al$. 2018', 'Rosato':'Rosato $et$ $al$. 2020', 'Syranidou':'Syranidou $et$ $al$. 2019', 'SyranidouPE':'Syranidou $et$ $al$. 2017a', 'SyranidouPS':'Syranidou $et$ $al$. 2017b', 'Tagg':'Tagg $et$ $al$. 2019', 'Woodall':'Woodall $et$ $al$. 2018', 'Wu':'Wu $et$ $al$. 2019', 'Xu':'Xu $et$ $al$. 2019', 'Zhang':'Zhang $et$ $al$. 2019', 'WaterOrSediment':'Water or Sediment', 'LabOrField':'Laboratory or Field', 'IncubationOrCollection':'Incubation or Collection', 'MaterialType':'Material type', 'PlasticTypeSpecific':'Plastic type (specific)', 'PlasticTypeGeneral':'Plastic type (general)', 'DEPTH':'Depth', 'IncubationTime':'Incubation time (specific)', 'IncubationGeneral':'Incubation time (general)', 'PrimerPair':'Primer pair', 'DNAExtraction':'DNA extraction method', 'lab':'Laboratory', 'not_plastic':'Not plastic', 'aged_oxope':'Aged Oxo-PE', 'freeliving':'Free living', 'particleassociated':'Particle associated', 'oxope':'Oxo-PE', 'rinse_pe':'PE rinse water', 'rinse_ps':'PS rinse water', 'rinse_wood':'Wood rinse water', 'bhet':'BHET', 'hdpe':'HDPE', 'ldpe':'LDPE', 'na':'NA', 'pa':'PA', 'pe':'PE', 'pes':'PES', 'pestur':'PESTUR', 'pet':'PET', 'phbv':'PHBV', 'pla':'PLA', 'pp':'PP', 'ps':'PS', 'pvc':'PVC', 'san':'SAN', 'w_pe':'Weathered PE', '10:14':'10:14 light:dark', '12:12':'12:12 light:dark', '16:08':'16:08 light:dark', '27F_519R':'27F-519R', '319F_806R':'319F-806R', '338F_806R':'338F-806R', '341F_785R':'341F-785R', '341F_802R':'341F-802R', '341F_806R':'341F-806R', '515F_806R':'515F-806R', '515FY_926R':'515FY-926R', '518F_926R':'518F-926R', '543F_783R':'543F-783R', '967F_1064R':'967F-1064R', 'B969F_BA1406R':'B969F-BA1406R', 'redextract_sigma':'REDExtract-$N$-AmpTM', 'gentra_puregene':'Gentra Puregene', 'purelink':'PureLink', 'powersoil':'PowerSoil', 'phenol_chloroform':'Phenol-Chloroform', 'powerbiofilm':'PowerBiofilm', 'ultraclean_soil':'UltraClean soil', 'fastdna_soil':'FastDNA soil', 'orders':'Order', 'classes':'Class', 'phyla':'Phylum', 'genera':'Genera', 'families':'Family', 'species':'Species', 'ASVs':'ASV', 'kingdoms':'Kingdom', 'PlasticOnly':'Plastic only', '534f_783r':'534F-783R', 'Phenol_chloroform':'Phenol-chloroform'}
meta_name_ind = {'Study':0, 'Latitude':1, 'Longitude':2, 'Environment':3, 'WaterOrSediment':4, 'LabOrField':5, 'IncubationOrCollection':6, 'Source':7, 'MaterialType':8, 'PlasticTypeSpecific':9, 'PlasticTypeGeneral':10, 'DEPTH':11, 'IncubationTime':12, 'IncubationGeneral':13, 'Temperature':14, 'Salinity':15, 'Light':16, 'Season':17, 'PrimerPair':18, 'DNAExtraction':19, 'PlasticOnly':20}
name_dict = {'La2020':'Latva\n$et$ $al$. 2020', 'AZ2015':'Amaral-Zettler\n$et$ $al$. 2015', 'AA2018':'Arias-Andres\n$et$ $al$. 2018', 'Ca2020':'Canada\n$et$ $al$. 2020', 'Cu2019':'Curren & Leong\n2019', 'De2019':'Delacuvellerie\n$et$ $al$. 2019', 'DT2015':'De Tender\n$et$ $al$. 2015', 'DT2017':'De Tender\n$et$ $al$. 2017', 'DH2018':'Dussud \n$et$ $al$. 2018a', 'DM2018':'Dussud\n$et$ $al$. 2018b', 'EC2019':'Erni-Cassola\n$et$ $al$. 2019', 'Es2019':'Esan\n$et$ $al$. 2019', 'Fr2018':'Frere\n$et$ $al$. 2018', 'Ho2014':'Hoellein\n$et$ $al$. 2014', 'Ho2017':'Hoellein\n$et$ $al$. 2017', 'Ji2018':'Jiang\n$et$ $al$. 2018', 'Ke2019':'Kesy\n$et$ $al$. 2019', 'Ki2018':'Kirstein\n$et$ $al$. 2018', 'Ki2019':'Kirstein\n$et$ $al$. 2019', 'MC2014':'McCormick\n$et$ $al$. 2014', 'MC2016':'McCormick\n$et$ $al$. 2016', 'Ob2016':'Oberbeckmann\n$et$ $al$. 2016', 'Ob2018':'Oberbeckmann\n$et$ $al$. 2018', 'Og2018':'Ogonowski\n$et$ $al$. 2018', 'Pi2019':'Pinto\n$et$ $al$. 2019', 'Po2018':'Pollet\n$et$ $al$. 2018', 'Ro2020':'Rosato\n$et$ $al$. 2020', 'Sy2019':'Syranidou\n$et$ $al$. 2019', 'SyPE20':'Syranidou\n$et$ $al$. 2017a', 'SyPS20':'Syranidou\n$et$ $al$. 2017b', 'Ta2019':'Tagg\n$et$ $al$. 2019', 'Wo2018':'Woodall\n$et$ $al$. 2018', 'Wr2019':'Wright\n$et$ $al$. 2020', 'Wu2019':'Wu\n$et$ $al$. 2019', 'Xu2019':'Xu\n$et$ $al$. 2019', 'Zh2019':'Zhang\n$et$ $al$. 2019', 'Br2016':'Bryant\n$et$ $al$. 2016', 'Pin201':'Pinnell\n$et$ $al$. 2019', 'Pa2019':'Parrish\n$et$ $al$. 2019'}
name_dict_2 = {'La2020':'Latva $et$ $al$. 2020', 'AZ2015':'Amaral-Zettler $et$ $al$. 2015', 'AA2018':'Arias-Andres $et$ $al$. 2018', 'Ca2020':'Canada $et$ $al$. 2020', 'Cu2019':'Curren & Leong 2019', 'De2019':'Delacuvellerie $et$ $al$. 2019', 'DT2015':'De Tender $et$ $al$. 2015', 'DT2017':'De Tender $et$ $al$. 2017', 'DH2018':'Dussud $et$ $al$.  2018a', 'DM2018':'Dussud $et$ $al$. 2018b', 'EC2019':'Erni-Cassola $et$ $al$. 2019', 'Es2019':'Esan $et$ $al$. 2019', 'Fr2018':'Frere $et$ $al$. 2018', 'Ho2014':'Hoellein $et$ $al$. 2014', 'Ho2017':'Hoellein $et$ $al$. 2017', 'Ji2018':'Jiang $et$ $al$. 2018', 'Ke2019':'Kesy $et$ $al$. 2019', 'Ki2018':'Kirstein $et$ $al$. 2018', 'Ki2019':'Kirstein $et$ $al$. 2019', 'MC2014':'McCormick $et$ $al$. 2014', 'MC2016':'McCormick $et$ $al$. 2016', 'Ob2016':'Oberbeckmann $et$ $al$. 2016', 'Ob2018':'Oberbeckmann $et$ $al$. 2018', 'Og2018':'Ogonowski $et$ $al$. 2018', 'Pi2019':'Pinto $et$ $al$. 2019', 'Po2018':'Pollet $et$ $al$. 2018', 'Ro2020':'Rosato $et$ $al$. 2020', 'Sy2019':'Syranidou $et$ $al$. 2019', 'SyPE20':'Syranidou $et$ $al$. 2017a', 'SyPS20':'Syranidou $et$ $al$. 2017b', 'Ta2019':'Tagg $et$ $al$. 2019', 'Wo2018':'Woodall $et$ $al$. 2018', 'Wr2019':'Wright $et$ $al$. 2020', 'Wu2019':'Wu $et$ $al$. 2019', 'Xu2019':'Xu $et$ $al$. 2019', 'Zh2019':'Zhang $et$ $al$. 2019', 'Br2016':'Bryant $et$ $al$. 2016', 'Pin201':'Pinnell $et$ $al$. 2019', 'Pa2019':'Parrish $et$ $al$. 2019'}
name_env = {'AZ2015':'marine', 'AA2018':'freshwater', 'Ca2020':'aquatic', 'Cu2019':'marine', 'De2019':'marine', 'DT2015':'marine', 'DT2017':'marine', 'DH2018':'marine', 'DM2018':'marine', 'EC2019':'marine', 'Es2019':'terrestrial', 'Fr2018':'marine', 'Ho2014':'freshwater', 'Ho2017':'freshwater', 'Ke2019':'aquatic', 'Ki2018':'marine', 'Ki2019':'marine', 'MC2014':'freshwater', 'MC2016':'freshwater', 'Ob2016':'marine', 'Ob2018':'aquatic', 'Og2018':'marine', 'Pi2019':'marine', 'Po2018':'marine', 'Ro2020':'marine', 'Sy2019':'marine', 'SyPE20':'marine', 'SyPS20':'marine', 'Ta2019':'aquatic', 'Wo2018':'marine', 'Wr2019':'marine', 'Wu2019':'freshwater', 'Xu2019':'marine', 'Zh2019':'terrestrial', 'Pa2019':'aquatic'}
norm_names = {'rare':'Rarefied', 'rel_abun':'Relative\nabundance', 'log':'Log', 'clr':'CLR'}

def open_csv(fn): #open csv with file name fn, returning the data as rows
    '''
    Opens csv file and returns the data as rows
    Takes as input:
        - fn (.csv file name)
    Returns:
        - rows of the file as a list
    '''
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    return rows

def open_txt(fn):
    '''
    Opens text file and returns the data as rows
    Takes as input:
        - fn (.txt file name - assumes that data is tab delimited and lines terminated with \n)
    Returns:
        - rows of the file as a list
    '''
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f, delimiter='\t', lineterminator='\n'):
            rows.append(row)
    return rows

def open_pickle(fn):
    '''
    Opens python data object file and returns the data as a python object
    Takes as input:
        - fn (file name for python object)
    Returns:
        - python object
    '''
    with open(fn, 'rb') as pick: 
        pick = pickle.load(pick)
    return pick

def write_csv(fn, data):
    '''
    Save a .csv file of the data given
    Takes as input:
       - fn (the name with which to save the csv file)
       - data (a list of lists - the data to go in the csv file)
    '''
    with open(fn, 'w') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)
    return

def write_txt(fn, data):
    '''
    Save a .txt file (tab delimited with \n to break lines) of the data given
    Takes as input:
       - fn (the name with which to save the txt file)
       - data (a list of lists - the data to go in the txt file)
    '''
    with open(fn, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for row in data:
            writer.writerow(row)
    return

def write_pickle(fn, data): #write pickle python object with name fn and the data in data
    '''
    Write a pickle python object with the data given
    Takes as input:
        - fn (the name with which to save the data)
        - data (the python object to save)
    '''
    with open(fn, 'wb') as f:
        pickle.dump(data, f)
    return
    
def open_and_sort(fn):
    '''
    This is a function to sort the columns of a .csv file and return the resulting sorted dataframe
    Takes as input:
        - fn (the name of the .csv file)
    If 'unifrac' is in the name then it will sort the rows as well as the columns. Otherwise only the columns will be sorted
    Returns:
        - df (a pandas dataframe of the sorted .csv file)
    '''
    df = pd.read_csv(fn, header=0, index_col=0) #open the file
    df.sort_index(inplace=True) #sort the column names
    if 'unifrac' in fn or 'distance' in fn: #if it is a distance matrix
        df.sort_index(axis=1, inplace=True) #also sort the rows
    return df

def get_meta(meta): #get the information contained in the meta file as a dictionary
    '''
    Function to get the information contained in the metadata file as a dataframe and a dictionary
    Takes as input:
        - meta (.txt tab delimited file with sample names as rows and metadata categories as columns)
    Returns:
        - meta (the metadata file as a list)
        - meta_names (the names of the metadata categories)
        - meta_dict (a dictionary containing sample names as keys and all metadata stored)
    '''
    meta = open_txt(meta) #open the file
    meta_names = meta[0] #save the column names
    del meta[0] #delete the column names
    meta_dict = {} #create a dictionary of all information contained, with sample names as dictionary keys
    for a in range(len(meta)):
        meta_dict[meta[a][0]] = meta[a][1:]
        for b in range(len(meta[a])):
            if meta[a][b] not in rename_plots:
                try:
                    float(meta[a][b])
                except:
                    rename_plots[meta[a][b]] = meta[a][b].capitalize()
    return meta, meta_names, meta_dict

def get_meta_df(meta, meta_names, ft_samples):
    '''
    Function to check that only the samples that were included in the QIIME analysis (after filtering, rarefaction, etc) remain in the metadata dataframe
    Takes as input:
        - meta (the metadata file as a list (one of the outputs of the get_meta function)
        - meta_names (list of metadata categories (also output from the get_meta function))
        - ft_samples (the sample names from the feature table)
    Returns:
        - meta_df (a dataframe containing only the samples that are in the feature table)
    '''
    meta_df = pd.DataFrame(meta, columns=meta_names) #open the meta file
    meta_df.index = meta_df['#SampleID'] #change the index to be the sample ID column
    meta_df.drop('#SampleID', axis=1, inplace=True) #remove this column
    sn_meta = list(meta_df.index.values) #get a list of all sample names
    for a in sn_meta: #for each of the samples
        if a not in ft_samples: #if it's not also in the samples list
            meta_df.drop(a, axis=0, inplace=True) #remove it from the meta dataframe
    return meta_df
    
def get_eucl_dist(ft, dist_met='euclidean'):
  X = ft.iloc[0:].values
  y = ft.iloc[:,0].values
  X_true = X
  similarities = distance.cdist(X_true, X_true, dist_met)
  return similarities
    
def get_ASV_dict(ft, seqs, basedir):
    '''
    Function to rename the ASVs with a number (rather than the random sequence of letters/numbers that is standard for QIIME2 output)
    Takes as input:
        - ft (feature table pandas dataframe with sample names as columns and ASVs as rows)
        - seqs (name of the sequences fasta file that contains sequences for all of the ASVs)
    Returns:
        - ASV_dict (a dictionary containing the previous ASV names as keys and the new ASV names as data)
    '''
    asvs = list(ft.index.values)
    ASV_dict = {}
    for a in range(len(asvs)):
        ASV_dict[asvs[a]] = 'ASV'+str(a).zfill(4)
    seqs_rename = [] #create a list to add the sequence records to
    for record in SeqIO.parse(seqs, "fasta"): #open the sequences file and loop through it
        if record.id in asvs: #if the record id matches one found in the asv list
            record.id = ASV_dict[record.id]
            seqs_rename.append(record) #add the record to the list
    SeqIO.write(seqs_rename, basedir+"sequences_agglom_renamed.fasta", "fasta") #save the new list of sequence records
    return ASV_dict
    
def mark_inset(parent_axes, inset_axes, loc1a=1, loc1b=1, loc2a=2, loc2b=2, **kwargs):
    rect = TransformedBbox(inset_axes.viewLim, parent_axes.transData)

    pp = BboxPatch(rect, fill=False, **kwargs)
    parent_axes.add_patch(pp)

    p1 = BboxConnector(inset_axes.bbox, rect, loc1=loc1a, loc2=loc1b, **kwargs)
    inset_axes.add_patch(p1)
    p1.set_clip_on(False)
    p2 = BboxConnector(inset_axes.bbox, rect, loc1=loc2a, loc2=loc2b, **kwargs)
    inset_axes.add_patch(p2)
    p2.set_clip_on(False)

    return pp, p1, p2
    
def transform_for_NMDS(similarities, n_jobs=1): #transform the similarity matrix to 2D space (n_components=2)
    '''
    For a similarity matrix, calculate the best fit of points in 2D 
    Takes as input:
        - similarities (a distance matrix with no row or sample names)
        - n_jobs (number of processors to use for constraining - doesn't seem to work so we leave this as the default 1)
    '''
    #similarities = np.nan_to_num(similarities)
    X_true = similarities.iloc[0:].values
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity="precomputed", n_jobs=n_jobs)
    similarities = np.nan_to_num(similarities)
    pos = mds.fit(similarities).embedding_
    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12, dissimilarity="precomputed", random_state=seed, n_jobs=n_jobs, n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)
    npos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((npos ** 2).sum())
    clf = PCA()
    npos = clf.fit_transform(npos)
    return npos, nmds.stress_

def format_R(ft, basedir, name): 
    '''
    Function that transforms the feature table output by QIIME2 to a version that can be used by the R agglomerate script (phyloseq and unifrac)
    Takes as input:
        - feature table.txt with samples as columns and ASVs as rows
    Gives as output:
        - ft as a pandas dataframe with sorted sample names
        - tax_dict - a dictionary containing the full taxonomy for each ASV
        - saves various files: feature_table.csv, taxonomy.csv, taxonomy_name_only.csv, random_forest/taxonomy_name_only.csv, tax_dict.dictionary
    '''
    ft = open_txt(ft) #open the text file
    del ft[0] #delete the first line of the file that tells you it was created from a .biom file
    tax_dict = {} #create a dictionary that will have all of the taxon information
    tax_file = [['OTUID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Species name']] #create a list that will have all of the taxon information
    tax_name_only = [['OTUID', 'Species name']]
    for a in range(len(ft)): #for each row of the feature table
        if a > 0: #if it's not the header row
            string = ft[a][-1].split(';') #split this taxon information based on semi-colons (the default format from QIIME2)
            new_string = []
            for b in range(len(string)):
                if b < 5 and 'Alteromonas' in string[b]:
                    break
                if 'uncultured' in string[b]:
                    continue
                elif 'Ambiguous' in string[b]:
                    continue
                elif 'Uncultured' in string[b]:
                    continue
                elif 'ambiguous' in string[b]:
                    continue
                else:
                    if string[b][:3] == ' D_':
                        string[b] = string[b][6:]
                    elif string[b][:2] == 'D_':
                        string[b] = string[b][5:]
                    new_string.append(string[b])
            sp_name = new_string[-1]
            while len(new_string) < 7:
                new_string.append('')
            new_string.append(sp_name)
            tax_name_only.append([ft[a][0], sp_name])
            tax_dict[ft[a][0]] = new_string #creates the taxon in the taxa dictionary
            tax_file.append([ft[a][0]]+new_string) #add this asv name and taxon information
        ft[a] = ft[a][:-1] #remove the taxon information from the feature table
    ft[0][0] = 'OTUID' #Rename the first cell to match the expected input in R
    write_csv(basedir+'feature_table_'+name+'.csv', ft) #write the new feature table to .csv
    ft = pd.read_csv(basedir+'feature_table_'+name+'.csv', header=0, index_col=0)
    ft.sort_index(inplace=True)
    ft.sort_index(axis=1, inplace=True)
    write_csv(basedir+'tax_dict_'+name+'.csv', tax_file) #write the taxonomy to .csv
    write_csv(basedir+'tax_dict_'+name+'_name_only.csv', tax_name_only) #write the taxonomy to .csv
    write_pickle(basedir+'tax_dict_'+name+'.dictionary', tax_dict) #write the taxon information to a python object
    return ft, tax_dict
    
def study_map_metrics(dates, locs, basedir, map_img, meta_df):
    '''
    Function to plot and save Figure 1 map (cumulative number of studies at the beginning of 2020 and the world map showing all of these studies) and metrics (number of samples included for each environment, plastic type and whether these were from incubations or collections, i.e. Figure 1C, 1D, 1E)
    It uses files that are not the metadata file in order to also include those studies with data that wasn't accessible
    Takes as input:
        - dates (file containing details of the numbers of studies published in each environment in each year, and whether the data for these studies are accessible, named Study_dates.csv in the data to be downloaded for this study)
        - locs (file containing details of study locations for all studies published up to 2020, including whether the data was accessible for each, named Study_location.csv in the data to be downloaded for this study)
        - basedir (location of the folder to save the figure to)
        - map_img (path to the map figure used as the background)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
    Returns:
        - nothing, but saves figure 'study_map' to the figures folder
    '''
    dates = open_csv(dates) #open csv file with study dates
    locs = open_csv(locs) #open csv file with study locations
    plt.figure(figsize=(15,8)) #set up figure
    ax1 = plt.subplot2grid((80,4), (0,0), rowspan=36) #add axis for number of publications
    ax2 = plt.subplot2grid((2,4), (0,1), colspan=2) #add axis for map
    ax1.set_title('A', loc='left', fontsize=fs_main, fontweight='bold') #add axis label
    ax2.set_title('B', loc='left', fontsize=fs_main, fontweight='bold') #add axis label
    plt.sca(ax2)
    img = plt.imread(map_img) #get the background map picture and set the main and inset axis to show only the regions of interest
    ax2.imshow(img, extent=[-180, 180, -90, 90], alpha=0.6)
    years, marine, freshwater, terrestrial, aquatic, marine_unav, freshwater_unav, terrestrial_unav, aquatic_unav = [], [], [], [], [], [], [], [], [] #set up lists to add data to
    mar, fw, terr, aqu, mar_unav, fw_unav, terr_unav, aqu_unav = 0, 0, 0, 0, 0, 0, 0, 0 #set all numbers to zero to start

    for a in range(1, len(dates)): #go through and see if there is data for each cell
        for b in range(len(dates[a])):
            if dates[a][b] == '': dates[a][b] = 0 #if there's nothing in the cell, make this equal to 0
            else: dates[a][b] = int(dates[a][b]) #if there is a number, make this into an integer
        years.append(dates[a][0]) #add year to list
        #add all cumulative numbers on
        mar += dates[a][1]
        mar_unav += dates[a][2]
        fw += dates[a][3]
        fw_unav += dates[a][4]
        aqu += dates[a][5]
        aqu_unav += dates[a][6]
        terr += dates[a][7]
        terr_unav += dates[a][8]
        #append the new cumulative numbers to the respective lists for each environment
        marine.append(mar), marine_unav.append(mar_unav), freshwater.append(fw), freshwater_unav.append(fw_unav), aquatic.append(aqu), aquatic_unav.append(aqu_unav), terrestrial.append(terr), terrestrial_unav.append(terr_unav)
    adding = [terrestrial, terrestrial_unav, aquatic, aquatic_unav, freshwater, freshwater_unav, marine, marine_unav]
    #now add all of the environments together
    for a in range(1, len(adding)):
        for b in range(len(adding[a])):
            adding[a][b] += adding[a-1][b]
    #set up the colors that we will use and the transparency for the studies with no data
    colors_2 = [color_env['terrestrial'], color_env['terrestrial'], color_env['aquatic'], color_env['aquatic'], color_env['freshwater'], color_env['freshwater'], color_env['marine'], color_env['marine']]
    a1, a2 = 1, 0.3
    alphas = [a1, a2, a1, a2, a1, a2, a1, a2]
    for a in range(len(adding)): #for each line that we have set up
        ax1.plot(years, adding[a], 'grey', lw=1) #plot the line in grey
        if a > 0: #and fill between that line and either the previous line or zero
            ax1.fill_between(years, adding[a-1], adding[a], facecolor=colors_2[a], alpha=alphas[a])
        else:
            ax1.fill_between(years, 0, adding[a], facecolor=colors_2[a], alpha=alphas[a])
    #show only the parts of the axis that we are interested in
    ax1.set_xlim([2012, 2020])
    ax1.set_ylim([0, 55])
    #add axis labels and change the font sizes to be consistent throughout
    ax1.set_ylabel('Cumulative number of\nPlastisphere publications', fontsize=fs_small)
    ax1.set_xlabel('Year', fontsize=fs_main)
    ax1.set_xticks([2012, 2014, 2016, 2018, 2020])
    ax1.tick_params(axis='both', which='major', labelsize=fs_small)
    ax1.tick_params(axis='both', which='minor', labelsize=fs_small)
    #make a custom legend (this is so we have boxes for colors rather than lines)
    marine = mlines.Line2D([], [], color=color_env['marine'], marker='s', markersize=8, markeredgecolor='k', label=label_env[0], linestyle=' ')
    freshwater = mlines.Line2D([], [], color=color_env['freshwater'], marker='s', markersize=8, markeredgecolor='k', label=label_env[1], linestyle=' ')
    aquatic = mlines.Line2D([], [], color=color_env['aquatic'], marker='s', markersize=8, markeredgecolor='k', label=label_env[2], linestyle=' ')
    terrestrial = mlines.Line2D([], [], color=color_env['terrestrial'], marker='s', markersize=8, markeredgecolor='k', label=label_env[3], linestyle=' ')
    plt.sca(ax1)
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial], loc='upper left', fontsize=fs_small)  
    
    #make all locations into floats (rather than the default strings)
    for a in range(1, len(locs)):
        for b in range(len(locs[a])):
            if b == 1 or b == 2 or b == 3:
                locs[a][b] = float(locs[a][b])

    axins1 = zoomed_inset_axes(ax2, 2.5, loc='upper left', bbox_to_anchor=(0.25, 0.5), bbox_transform=ax2.transAxes)
    axins1.imshow(img, extent=[-180, 180, -90, 90], alpha=0.6)
    axins1.set_xlim([-22, 30]), axins1.set_ylim([25, 65])
    mark_inset(ax2, axins1, loc1a=1, loc1b=4, loc2a=2, loc2b=3, fc="none", ec="k")
    
    axins2 = zoomed_inset_axes(ax2, 2.1, loc='upper right', bbox_to_anchor=(0.99, 0.45), bbox_transform=ax2.transAxes)
    axins2.imshow(img, extent=[-180, 180, -90, 90], alpha=0.6)
    axins2.set_xlim([105, 130]), axins2.set_ylim([15, 45])
    mark_inset(ax2, axins2, loc1a=1, loc1b=4, loc2a=2, loc2b=3, fc="none", ec="k")
    
    for ax in [ax2, axins1, axins2]:
      plt.sca(ax)
      plt.yticks([]), plt.xticks([])

    #plot each location, checking which color and marker these need based on the values in the other columns of the file
    for i in range(1, len(locs)):
        if locs[i][5] == 'Yes': alpha, edge = 1, 'k'
        else: alpha, edge = 1, 'w'
        if locs[i][6] == 'Lab': marker='^'
        elif locs[i][6] == 'Field': marker = 'o'
        elif locs[i][6] == 'Both': marker = '*'
        ax2.scatter(locs[i][2], locs[i][1], color=color_env[locs[i][4]], edgecolor=edge, marker=marker, s=15, linewidth=0.4, alpha=alpha)
        axins1.scatter(locs[i][2], locs[i][1], color=color_env[locs[i][4]], edgecolor=edge, marker=marker, s=25, linewidth=0.4, alpha=alpha)
        axins2.scatter(locs[i][2], locs[i][1], color=color_env[locs[i][4]], edgecolor=edge, marker=marker, s=20, linewidth=0.4, alpha=alpha)
    
    #again add a custom legend
    plt.sca(ax2)
    lab = mlines.Line2D([], [], color='w', marker='^', markersize=8, markeredgecolor='k', label='Laboratory', linestyle=' ')
    field = mlines.Line2D([], [], color='w', marker='o', markersize=8, markeredgecolor='k', label='Field', linestyle=' ')
    both = mlines.Line2D([], [], color='w', marker='*', markersize=8, markeredgecolor='k', label='Both', linestyle=' ')
    gap = mlines.Line2D([], [], color='w', marker='o', markersize=2, markeredgecolor='w', alpha=0, label='', linestyle=' ')
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial, gap, lab, field, both], loc='lower left', fontsize=fs_small)  

    #now save the figure (using the file extension and dpi specified at the top of the file) and close it

    #set up the figure and plots
    #plt.figure(figsize=(10,3))
    ax1 = plt.subplot2grid((2,4), (1,0))
    ax2 = plt.subplot2grid((2,4), (1,1))
    ax3 = plt.subplot2grid((2,4), (1,2))
    samples = list(meta_df.index.values) #get all sample names
    envs, env_caps = ['terrestrial', 'aquatic', 'freshwater', 'marine'], ['Terrestrial', 'Aquatic', 'Freshwater', 'Marine'] #set up lists of the environment names
    for a in range(len(envs)): #for each of the environments
        bottom = 0 #reset the bottom
        for b in ['lab', 'field']: #for each of lab and field
            count = 0 #reset the count
            for c in range(len(samples)): #for each of the sample names
                if meta_df.loc[samples[c], 'Environment'] == envs[a] and meta_df.loc[samples[c], 'LabOrField'] == b: #if the metadata shows the sample is in the environment and in the lab/field category we are looking at
                    count += 1 #add it to the count
            if b == 'lab': #if we are looking at lab samples
                ax1.barh(a, count, left=bottom, color=color_env[envs[a]], hatch='//', edgecolor='k') #plot a bar with the correct colors and hatching
            else: #otherwise
                ax1.barh(a, count, left=bottom, color=color_env[envs[a]], edgecolor='k') #use the right color, but not hatched
            if count > 100: #if the count is above 100 (i.e. the bar is large enough to show text within)
                ax1.text((count/2)+bottom, a, str(count), color='w', ha='center', va='center', fontsize=fs_small) #add text saying the number of samples in this category
            bottom = count+bottom #now add this count to the bottom
            if b == 'field': #if we are looking at field samples (i.e. we have now plotted both bars for this environmet)
                ax1.text(bottom+20, a, str(bottom), color='k', ha='left', va='center', fontsize=fs_small) #print the total number of samples in this environment next to the bar
        pie, inc = [], [] #set up new lists
        for c in range(len(samples)): #for each of the samples
            if meta_df.loc[samples[c], 'Environment'] == envs[a]: #if this is the environment we are looking at
                pie.append(meta_df.loc[samples[c], 'PlasticTypeGeneral']) #add the plastic type to the list
                inc.append(meta_df.loc[samples[c], 'IncubationGeneral']) #add the incubation time to the other list
        unique = list(set(pie)) #get all of the unique values for plastic type in this environment
        pie_color = [color_source[x] for x in unique] #and get the colors that match each plastic type from the dictionary at the top
        pie_count = [pie.count(x) for x in unique] #and count how many samples of each plastic type
        pie_count = [(x/sum(pie_count))*100 for x in pie_count] #and get these counts as relative abundance
        left=0 #reset the left count
        for d in range(len(pie_count)): #for each sample type
            ax2.barh(a, pie_count[d], left=left, color=pie_color[d], label=unique[d], edgecolor='k') #plot the stacked bar corresponding to sample type
            fc = 'w' #set the fontcolor as white
            if unique[d] == 'planktonic' or unique[d] == 'biofilm': #unless it's planktonic or biofilm
                fc = 'k' #in which case, the fontcolor should be black
            if pie_count[d] > 5: #if this category has more than 5% of samples (i.e. is big enought to add text)
                ax2.text((pie_count[d]/2)+left, a, str(int(pie_count[d])), color=fc, ha='center', va='center', fontsize=fs_small) #add text showing the percentage of samples with this plastic type
            left += pie_count[d] #now add this number to the left count (to make the bar stacked)
        unique = list(set(inc)) #now we are doing the same for incubation time, so get the unique values for incubation time in this environment
        if len(unique) == 3: #if this list is equal to 3
            unique = ['early', 'late', 'collection'] #change the order to be this
        pie_color = [color_source[x] for x in unique] #and get the colors that match these incubation times
        pie_count = [inc.count(x) for x in unique] #count how many samples for each incubation time
        pie_count = [(x/sum(pie_count))*100 for x in pie_count] #and convert these counts to relative abundance
        left = 0 #reset the left count
        for d in range(len(pie_count)): #for each sample type
            ax3.barh(a, pie_count[d], left=left, color=pie_color[d], label=unique[d], edgecolor='k') #plot the stacked bar corresponding to sample type
            fc = 'w' #set the fontcolor as white
            if unique[d] == 'early': #unless it's an early incubation time
                fc = 'k' #in which case the fontcolor should be black
            if pie_count[d] > 5: #if this accounts for more than 5% of relative abundance (i.e. is big enough to add text)
                ax3.text((pie_count[d]/2)+left, a, str(int(pie_count[d])), color=fc, ha='center', va='center', fontsize=fs_small) #add text showing the percentage of samples with this incubation time
            left += pie_count[d] #now add the number to the left count 
    #set x and y labels, ticks, limits, etc.  and make patches and legends for each plot
    ax1.set_xlabel('Number of samples', fontsize=fs_small)
    plt.sca(ax1)
    plt.yticks([0, 1, 2, 3], env_caps, fontsize=fs_small)
    plt.xticks(fontsize=fs_small)
    plt.xlim([0, 1400])
    lab = patches.Patch(facecolor='w', edgecolor='k', hatch='//', label='Laboratory')
    field = patches.Patch(facecolor='w', edgecolor='k', label='Field')
    ax1.legend(handles=[lab, field], fontsize=fs_main, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=2)
    sources = ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic', 'blank']
    handles = [patches.Patch(facecolor=color_source[x], edgecolor='k', label=x.capitalize()) for x in sources]
    ax2.legend(handles=handles, fontsize=fs_small, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=2)
    plt.sca(ax2)
    plt.yticks([])
    plt.xticks(fontsize=fs_small)
    plt.xlabel('Relative abundance (%)', fontsize=fs_small)
    plt.xlim([0, 100])
    times = ['early', 'late', 'collection']
    handles = [patches.Patch(facecolor=color_source[x], edgecolor='k', label=x.capitalize()) for x in times]
    ax3.legend(handles=handles, fontsize=fs_small, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=3)
    plt.sca(ax3)
    plt.yticks([])
    plt.xticks(fontsize=fs_small)
    plt.xlabel('Relative abundance (%)', fontsize=fs_small)
    plt.xlim([0, 100])
    plt.subplots_adjust(wspace=0.1)
    ax1.set_title('C', loc='left', fontsize=fs_main, fontweight='bold')
    ax2.set_title('D', loc='left', fontsize=fs_main, fontweight='bold')
    ax3.set_title('E', loc='left', fontsize=fs_main, fontweight='bold')
    plt.savefig(basedir+'/figures/'+'Fig1_map_metrics'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def get_single_nmds(rows, filter_on, filt_ind, color_on, color_ind, ax, leg, colors, names, meta_dict, save_name, second_filter='', second_filter_ind='', npos='', n_jobs=1, get_stats=False):
    '''
    Function to get a single nmds plot on an axis
    Takes as input:
        - rows (a distance matrix with no row or column names, i.e. containing only numbers)
        - filter_on (if we want to filter out samples from the matrix, then this contains what we are filtering on)
        - filter_ind (and the index in the metadata that we are filtering on)
        - color_on (the category that we are coloring samples on)
        - color_ind (the index for this category in the metadata)
        - ax (the axis we are plotting on)
        - leg (the position of the index we are plotting)
        - colors (the colors to use for plotting categories)
        - names (the sample names, in the same order as the samples occur in rows, the distance matrix)
        - meta_dict (a dictionary containing all metadata for each sample, with sample names as keys)
        - second_filter (if we want to filter samples based on a second category, then this contains what we are filtering on. Default is '')
        - second_filter_ind (and the index in the metadata that we are filtering on. Default is '')
        - npos (the nmds positions - this saves calculating the nmds points for a second time if we are plotting the same data again, e.g. just with different colored points)
        - n_jobs (the number of processors to use for calculating the nmds plot - this doesn't seem to work if we put more than 1, so isn't actually used)
        - get_stats (this tells us whether to get stats based on our categories - e.g. permanova and anosim - default is False)
    '''
    plt.sca(ax)
    color = []
    #If we are filtering, go through and save only the samples that are included into new lists
    if filter_on != 'none':
        new_names, new_ind = [], []
        for a in range(len(names)):
            if meta_dict[names[a]][filt_ind] == filter_on:
                if second_filter != '':
                    if meta_dict[names[a]][second_filter_ind] != second_filter:
                        continue
                new_names.append(names[a])
                new_ind.append(a)
        new_rows = []
        if new_names == []:
            return
        for b in range(len(rows)):
            if b in new_ind:
                this_row = []
                for c in range(len(rows[b])):
                    if c in new_ind:
                        this_row.append(rows[b][c])
                new_rows.append(this_row)
        names = new_names
        rows = new_rows
    #get the appropriate colors for all included samples
    groups = []
    for a in range(len(names)):
        added = False
        groups.append(meta_dict[names[a]][color_ind])
        for b in range(len(color_on)):
            if meta_dict[names[a]][color_ind] == color_on[b]:
                color.append(colors[b])
                added = True
        if not added:
            print('Didnt find ', color_on[b], 'in the metadata file so no color for this or in the list of colors')
    #if we don't already have the values for npos, then transform the similarity matrix for this
    if npos == '':
        s = pd.DataFrame(rows).fillna(value=0)
        if os.path.exists(save_name):
          npos = open_pickle(save_name)
        else:
          npos, stress = transform_for_NMDS(s, n_jobs=n_jobs)
          write_pickle(save_name, npos)
    else:
        s = pd.DataFrame(rows).fillna(value=0)
    if get_stats != False:
        string = 'ANOSIM: R='+str(round(get_stats[0], 3))+r', $p$='+str(round(get_stats[1], 3))
        string += '\nPERMANOVA: F='+str(round(get_stats[2], 3))+r', $p$='+str(round(get_stats[3], 3))
        t = plt.text(0.02, 0.98, string, ha='left', va='top', transform = ax.transAxes, fontsize=fs_main, fontweight='bold')
        t.set_bbox(dict(facecolor='white', alpha=0.8))
    size = 20
    #plot the values for nmds 1 (npos[a,0]) and nmds 2 (npos[a,1]) for all samples, giving the color as determined by sample type
    for a in range(len(rows[0])):
        plt.scatter(npos[a,0], npos[a,1], color=color[a], marker='o', s=size, edgecolor='k')
    #get the legend handles incase these are being plotted
    
    handles = []
    for a in range(len(color_on)):
        if color_on[a] in rename_plots:
            this_color = rename_plots[color_on[a]]
        else:
            this_color = color_on[a]
        handles.append(mlines.Line2D([], [], color=colors[a], marker='o', markersize=fs_main, markeredgecolor='k', label=this_color, linestyle=' '))    
    if filter_on != 'none':
        return handles
    else:
        return npos, handles

def nmds_plot_study_env(dist_matr_fn_w, dist_matr_fn_uw, meta_dict, basedir, n_jobs=1, save_name='/figures/Fig2_nmds_overall_', stats=False, norm_name='rare'):
    '''
    Function to make nmds plots for weighted and unweighted unifrac distances for all samples, either colored by environment or by study
    Takes as input:
        - dist_matr_fn_w (file name of the weighted unifrac - can be .csv or .txt - assumes this is a file with sample names as rows and columns with values showing their similarity in the matrix)
        - dist_matr_fn_uw (file name of the unweighted unifrac - as above)
        - meta_dict (dictionary containing metadata for all samples, with sample names as keys)
        - basedir (name of the directory to save the figures to)
        - n_jobs (number of processors to use for calculating the best ordination - doesn't seem to work so we leave with the default of 1)
    Returns:
        - nothing, but saves figure 'nmds_overall' to the figures folder
    '''
    plt.figure(figsize=(15,15))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (0,1))
    ax3 = plt.subplot2grid((2,2), (1,0))
    ax4 = plt.subplot2grid((2,2), (1,1))
    
    dist = [dist_matr_fn_w, dist_matr_fn_uw]
    ax = [[ax1, ax2], [ax3, ax4]]
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    color_env = ['#03A498', '#03B1FC', '#9A75FC', '#FD9A64']
    study, study_index = ['AmaralZettler','AriasAndres','Canada','Curren','Delacuvellerie','DeTender','DeTenderB','DussudHudec','DussudMeistertzheim','ErniCassola','Esan','Frere','Hoellein','HoelleinB','Jiang','Kesy','Kirstein','KirsteinB','McCormick','McCormickB','Oberbeckmann','OberbeckmannB','Ogonowski','Parrish', 'Pinto','Pollet','Rosato','Syranidou','SyranidouPE','SyranidouPS','Tagg','Woodall','Wu','Xu','Zhang'], 0
    colormap_40b, colormap_40c, colormap_40a = mpl.cm.get_cmap('tab20b', 256), mpl.cm.get_cmap('tab20c', 256), mpl.cm.get_cmap('tab20', 256)
    norm, norm2, norm3 = mpl.colors.Normalize(vmin=0, vmax=19), mpl.colors.Normalize(vmin=20, vmax=39), mpl.colors.Normalize(vmin=40, vmax=59)
    m, m2, m3 = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40b), mpl.cm.ScalarMappable(norm=norm2, cmap=colormap_40c), mpl.cm.ScalarMappable(norm=norm3, cmap=colormap_40a)
               
    for z in [0, 1]:
        dist_matr = pd.read_csv(basedir+dist[z], header=0, index_col=0) #read in the distance matrix
        dist_matr = dist_matr.astype('float').fillna(value=0) #turn all values into floats
        names = list(dist_matr.columns) #get a list of column names
        dist_matr = dist_matr.rename_axis('ID').values #and turn it into a matrix with no row or column names
        color_study = []
        for a in range(len(study)):
            if a < 20:
                color_study.append(m.to_rgba(a))
            elif a < 40:
                color_study.append(m2.to_rgba(a))
            else:
                color_study.append(m3.to_rgba(a))
        filter_on, filter_index = 'none', 'none'
        second_filter, second_filter_ind = '', ''
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, envs, env_index, ax[z][0], 'upper left', color_env, names, meta_dict, basedir+'nmds_'+norm_name+'.df'+str(z), second_filter, second_filter_ind, '', n_jobs=n_jobs, get_stats=stats[z][0])
        if z == 0:
            plt.sca(ax[z][0])
            plt.legend(handles=handles, loc='upper right', fontsize=fs_main, edgecolor='k')
        else:
            plt.sca(ax[z][0])
            plt.legend(handles=handles, loc='upper right', fontsize=fs_main, edgecolor='k')
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, study, study_index, ax[z][1], 'upper right', color_study, names, meta_dict, basedir+'nmds_'+norm_name+'.df'+str(z), second_filter, second_filter_ind, npos, n_jobs=n_jobs, get_stats=stats[z][1])
        if z == 0:
            plt.sca(ax[z][1])
            plt.legend(handles=handles, bbox_to_anchor=(1.05,1.025), fontsize=fs_main, edgecolor='k')
    titles = [r'$\bf{Environment}$', r'$\bf{Study}$', '', '']
    title_letter = ['A', 'B', 'C', 'D']
    axes = [ax1, ax2, ax3, ax4]
    if norm_name == 'clr':
      ylabs, xlabs = [r'$\bf{Aitchison}$ $\bf{distance}$'+'\nnMDS2', '', r'$\bf{PHILR}$ $\bf{distance}$'+'\nnMDS2', ''], ['', '', 'nMDS1', 'nMDS1']
    else:
      ylabs, xlabs = [r'$\bf{Weighted unifrac}$ $\bf{distance}$'+'\nnMDS2', '', r'$\bf{Unweighted unifrac}$ $\bf{distance}$'+'\nnMDS2', ''], ['', '', 'nMDS1', 'nMDS1']
    for a in range(len(titles)):
        axes[a].set_title(titles[a], fontsize=fs_title)
        axes[a].set_title(title_letter[a], fontsize=fs_title, loc='left')
        axes[a].set_xlabel(xlabs[a], fontsize=fs_title)
        axes[a].set_ylabel(ylabs[a], fontsize=fs_title)
    plt.savefig(basedir+save_name+norm_name+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def plot_box(ax, l, r, b, t, line_col):
    '''
    Function to plot a box on a given axes
    Takes as input:
        - ax (the axis to plot on)
        - l, r, b, t (left, right, bottom and top coordinates)
        - line_col (the color to use to plot the line)
    '''
    plt.sca(ax)
    plt.plot([l, r], [b, b], color=line_col, lw=2)
    plt.plot([l, r], [t, t], color=line_col, lw=2)
    plt.plot([l, l], [t, b], color=line_col, lw=2)
    plt.plot([r, r], [t, b], color=line_col, lw=2)
    return

def similarity_heatmap_combined(dist_matr_fn_w, dist_matr_fn_uw, basedir, save_name='/figures/Fig3_unifrac_heatmap_combined_', norm_name='rare'):
    '''
    Function to plot the heatmap showing both weighted and unweighted unifrac distances for all studies
    Takes as input:
        - dist_matr_fn_w (file name of the weighted unifrac - can be .csv or .txt - assumes this is a file with sample names as rows and columns with values showing their similarity in the matrix)
        - dist_matr_fn_uw (file name of the unweighted unifrac - as above)
        - basedir (name of the directory to save the figures to)
    Returns:
        - nothing, but saves figure 'unifrac_heatmap_combined' to the figures folder
    '''
    dist_matr_files = [dist_matr_fn_w, dist_matr_fn_uw]
    fig = plt.figure(figsize=(12,12))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    ax = [ax1, ax2]
    for z in range(len(dist_matr_files)):
        plt.sca(ax[z])
        dist_matr = dist_matr_files[z]
        dist_matr = dist_matr.reset_index().drop('index', axis=1)
        dist_matr = [list(dist_matr.columns)]+dist_matr.values.tolist()
        studies, first = [], [] #set up lists to store the names of all studies as well as the index for when these start within the distance matrices
        for a in range(len(dist_matr[0])): #loop through the study names
            if dist_matr[0][a] == '': #if the cell is empty, continue
                continue
            if dist_matr[0][a][:6] not in studies: #if the study part of the name (excluding sample number) isn't in the list, add it
                studies.append(dist_matr[0][a][:6])
                first.append(a) #and also add the index for where this study starts
        marine, freshwater, aquatic, terrestrial = [], [], [], [] #set up lists to store which samples are from which environment
        for a in range(len(studies)): #loop through all study names
            #for each, get the environment from the dictionary at the top
            #append the study name to the correct environment if it matches one of the environment types
            if name_env[studies[a]] == 'marine':
                marine.append(a)
            elif name_env[studies[a]] == 'freshwater':
                freshwater.append(a)
            elif name_env[studies[a]] == 'aquatic':
                aquatic.append(a)
            elif name_env[studies[a]] == 'terrestrial':
                terrestrial.append(a)
        order = marine+aquatic+freshwater+terrestrial #the order for the studies to be plotted (so that environments are grouped together)
        lens = [len(marine), len(aquatic), len(freshwater), len(terrestrial)] #get the length (/number of studies) of each of the environments
        text_cols = [] #set up a list for the text colors to be added to so these match the environment type
        for a in range(len(marine)):
            text_cols.append(color_env['marine'])
        for a in range(len(aquatic)):
            text_cols.append(color_env['aquatic'])
        for a in range(len(freshwater)):
            text_cols.append(color_env['freshwater'])
        for a in range(len(terrestrial)):
            text_cols.append(color_env['terrestrial'])
        new_text_cols = [] #now rearrange them so they are in the order they will be plotted
        for a in range(len(order)):
            for b in range(len(order)):
                if a == order[b]:
                    new_text_cols.append(text_cols[b])
        each_study, each_row = [], []
        #go  through the distance matrix and calculate means for the distance between samples within and between different studies
        for a in range(1, len(dist_matr)):
            this_num, this_row = [], []
            for b in range(1, len(dist_matr[a])):
                if b in first or b == len(dist_matr[a])-1:
                    if b == len(dist_matr[a])-1:
                        this_num.append(float(dist_matr[a][b]))
                    if this_num == []:
                        continue
                    this_row.append(np.mean(this_num))
                    this_num = []
                else:
                    if float(dist_matr[a][b]) != 0:
                        this_num.append(float(dist_matr[a][b]))
            if a in first or a == len(dist_matr)-1:
                if a == len(dist_matr)-1:
                    each_row.append(this_row)
                if each_row == []:
                    continue
                this_study = []
                for c in range(len(each_row[0])):
                    this_num = []
                    for d in range(len(each_row)):
                        this_num.append(each_row[d][c])
                    this_study.append(np.mean(this_num))
                each_study.append(this_study)
                each_row = []
            actual_row = []
            for d in range(len(order)):
                actual_row.append(this_row[order[d]])
            each_row.append(actual_row)
        #now change the order of the studies so that they are grouped by environment
        each_study_new, studies_new = [], []
        for a in range(len(order)):
            each_study_new.append(each_study[order[a]])
            studies_new.append(studies[order[a]])
        #reverse these as they will be plotted from the bottom not the top
        each_study_new.reverse()
        studies_new.reverse()
        each_study, studies = each_study_new, studies_new
        
        #get the minimum and maximum distances between studies so that colors can be normalised to this range
        colormap = mpl.cm.get_cmap('inferno_r', 256)                
        ma, mi = [], []
        for a in range(len(each_study)):
            ma.append(max(each_study[a]))
            mi.append(min(each_study[a]))
        norm = mpl.colors.Normalize(vmin=min(mi), vmax=max(ma))
        m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
        #set up figure and axes
        ally = []
        #now go through and get the equivalent color for each mean distance and plot these as a bar plot
        means = []
        for a in range(len(each_study)):
            colors = []
            x, y, bottom = [], [], []
            mean = []
            for b in range(len(each_study[a])):
                if b != a:
                    mean.append(each_study[a][b])
                color = m.to_rgba(each_study[a][b])
                colors.append(color)
                x.append(b+1)
                y.append(1)
                bottom.append(a)
                plt.bar(b+1, 1, bottom=a, color=color, edgecolor='k', width=1)
            means.append(np.mean(mean))
            ally.append(a+0.5)
        cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colormap),ax=ax[z], shrink=0.4, pad=0.02, orientation='vertical')
        cb.ax.tick_params(labelsize=fs_small)
        if z == 0:
            if norm_name == 'clr':
              cb.set_label('Aitchison distance', fontsize=fs_main)
              il = 'Aitchison distance'
            else:
              cb.set_label('Weighted unifrac distance', fontsize=fs_main)
              il = 'Weighted unifrac distance'
        else:
            if norm == 'clr':
              cb.set_label('PHILR distance', fontsize=fs_main)
              il = 'PHILR distance'
            else:
              cb.set_label('Unweighted unifrac distance', fontsize=fs_main)
              il = 'Unweighted unifrac distance'
        #change x and y limits
        plt.xlim([0.5,ally[-1]+1])
        plt.ylim([0,ally[-1]+0.5])
        colors = ['#03A498', '#03B1FC', '#9A75FC', '#FD9A64']
        line_col = 'w'
    
        #plot white boxes around samples of the same environments
        l, r, b, t = 0.54, lens[0]+0.43, ally[-1]-lens[0]+0.57, ally[-1]+0.48#left, right, bottom, top
        plot_box(ax[z], l, r, b, t, line_col)
        l, r, b, t = l+lens[0], r+lens[1], b-lens[1], t-lens[0]
        plot_box(ax[z], l, r, b, t, line_col)
        l, r, b, t = l+lens[1], r+lens[2], b-lens[2], t-lens[1]
        plot_box(ax[z], l, r, b, t, line_col)
        l, r, b, t = l+lens[2], r+lens[3], b-lens[3], t-lens[2]
        plot_box(ax[z], l, r, b, t, line_col)
        for a in range(len(studies)): #rename the studies to a full study name rather than the short code they were given
            studies[a] = name_dict_2[studies[a]]
        new_means = pd.DataFrame([means], columns=studies, index=[il])
        if z == 0:
            means_df = new_means
        else:
            means_df = pd.concat([means_df, new_means])
        plt.sca(ax[z])
        empty = []
        text_cols.reverse()
        for a in range(len(ally)): #plot the study names with the correct colors  on the y axis
            plt.text(0, ally[a], studies[a], color=text_cols[a], va='center', ha='right', fontsize=fs_small)
            empty.append('')
        plt.yticks(ally, empty, fontsize=fs_small) #remove the y ticks
        if z == 0:
            studies.reverse()
            for a in range(len(x)): #not plot the study names on the x axis (the list is reversed as we plotted from the bottom up)
                plt.text(x[a], ally[-1]+1.5, studies[a], rotation=90, color=text_cols[-(a+1)], va='bottom', ha='center', fontsize=fs_small)
            plt.xticks(x, empty, fontsize=fs_small) #remove the x ticks
            ax[z].xaxis.tick_top()
        else:
            plt.xticks([])
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(basedir+save_name+norm_name+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return
    
def get_venn_labels(data): 
    '''
    Calculate the number of overlapping ASVs between different samples
    It takes as input:
        - data (which should be a list of lists, containing the taxa present in different samples)
    '''
    N = len(data)
    sets_data = [set(data[i]) for i in range(N)]
    s_all = set(chain(*data))
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value
    labels = []
    for k in set_collections:
        labels.append(str(len(set_collections[k])))
    return labels

def get_venn(labels, ax, colors, names=['A', 'B', 'C', 'D', 'E'], plotting=[True, True, True, True, True]):
    '''
    This is a function to get a venn diagram with up to 5 ellipses (only plotted if that sample type is present, but the layout remains the same)
    It takes as input:
        - labels (a list of labels as output by the get_venn_labels function)
        - ax (axis to plot the venn diagram on)
        - colors (list of colors to use for each ellipse, expected that this will have length 5)
        - names (list of names to use for each ellipse, expected that this will have length 5. Default=['A', 'B', 'C', 'D', 'E'])
        - plotting (list of True/False for whether to plot each ellipse, expected that this will have length 5. Default is plotting=[True, True, True, True, True])
    '''
    #list of locations for each point of the ellipses
    el1, el2, el3, el4, el5 = [0.428, 0.469, 0.558, 0.578, 0.489], [0.449, 0.543, 0.523, 0.432, 0.383], [0.87, 0.87, 0.87, 0.87, 0.87], [0.50, 0.50, 0.50, 0.50, 0.50], [155.0, 82.0, 10.0, 118.0, 46.0]
    for a in range(len(el1)): #plot all ellipses 
        if plotting[a]: #only if this value is true for each
            e = patches.Ellipse(xy=(el1[a], el2[a]), width=el3[a], height=el4[a], angle=el5[a], color=colors[a], alpha=0.4)
            ax.add_patch(e)
    #list of text locations
    x = [0.27, 0.72, 0.55, 0.91, 0.78, 0.84, 0.76, 0.51, 0.39, 0.42, 0.50, 0.67, 0.70, 0.51, 0.64, 0.10, 0.20, 0.76, 0.65, 0.18, 0.21, 0.81, 0.74, 0.27, 0.34, 0.33, 0.51, 0.25, 0.28, 0.36, 0.51]
    y = [0.11, 0.11, 0.13, 0.58, 0.64, 0.41, 0.55, 0.90, 0.15, 0.78, 0.15, 0.76, 0.71, 0.74, 0.67, 0.61, 0.31, 0.25, 0.23, 0.50, 0.37, 0.37, 0.40, 0.70, 0.25, 0.72, 0.22, 0.58, 0.39, 0.66, 0.47]
    for a in range(len(x)):
        if labels[a] != '0': #only plot the overlap number if it is bigger than 0 (it will probably be 0 only if that sample type wasn't present and we aren't plotting the ellipse anyway)
            ax.text(x[a], y[a], labels[a], horizontalalignment='center',verticalalignment='center',fontsize=fs_small,color="black")
    #list of text locations
    x, y = [0.02, 0.72, 0.97, 0.88, 0.12], [0.72, 0.94, 0.74, 0.05, 0.05]
    for a in range(len(x)):
        if plotting[a]: #add the sample type name if the value for plotting is true
            ax.text(x[a], y[a], names[a], horizontalalignment='center',verticalalignment='center',fontsize=fs_main,color="black")
    return


def get_diversity(diversity, sample):
    '''
    function to calculate a range of different diversity metrics
    It takes as input:
        - diversity (the name of the diversity metric we want, can be 'Simpsons', 'Shannon', 'Richness', 'Evenness', 'Maximum' (Maximum is not a diversity metric, the function will just return the maximum abundance value given in sample)
        - sample (a list of abundance values that should correspond to one sample)
    Returns:
        - The diversity index for the individual sample
    '''
    for a in range(len(sample)):
        sample[a] = float(sample[a])
    total = sum(sample)
    if diversity == 'Simpsons':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)**2
        simpsons = 1-(sum(sample))
        return simpsons
    elif diversity == 'Shannon':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)
            if sample[b] != 0:
                sample[b] = -(sample[b] * (np.log(sample[b])))
        shannon = sum(sample)
        return shannon
    elif diversity == 'Richness':
        rich = 0
        for b in range(len(sample)):
            if sample[b] != 0:
                rich += 1
        return rich
    elif diversity == 'Evenness':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)
            if sample[b] != 0:
                sample[b] = -(sample[b] * (np.log(sample[b])))
        shannon = sum(sample)
        rich = 0
        for b in range(len(sample)):
            if sample[b] != 0:
                rich += 1
        even = shannon/(np.log(rich))
        return even
    elif diversity == 'Maximum':
        ma = (max(sample)/total)*100
        return ma
    return

def get_cols(num):
    '''
    This is a function to get a list of either 20 or 40 colors depending on the input 'num' (assumed that this is a number between 0 and 40)
    '''
    colormap_20 = mpl.cm.get_cmap('tab20', 256)
    norm = mpl.colors.Normalize(vmin=0, vmax=19)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_20)
    colors_20 = []
    for a in range(20):
        colors_20.append(m.to_rgba(a))
    colormap_40b = mpl.cm.get_cmap('tab20b', 256)
    norm = mpl.colors.Normalize(vmin=0, vmax=19)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40b)
    colors_40 = []
    for a in range(20):
        colors_40.append(m.to_rgba(a))
    colormap_40c = mpl.cm.get_cmap('tab20c', 256)
    norm = mpl.colors.Normalize(vmin=20, vmax=39)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40c)
    for a in range(20):
        a += 20
        colors_40.append(m.to_rgba(a))
    if num < 21:
        return colors_20
    else:
        return colors_40

def bar_dendro_venn(ft, ft_full, meta_dict, basedir, tax_dict, str_norm='rare'):
    '''
    Function to make the dendrogram, stacked bar, heatmap, simpsons diversity and venn diagrams of shared ASVs
    Takes as input:
        - ft (dataframe with samples as columns and ASVs as rows)
        - meta_dict (dictionary containing metadata for all samples, with sample names as keys)
        - basedir (name of the directory to save the figures to)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
    Returns:
        - nothing, but saves figure to the figures folder
    '''
    #set up the figure and all axes
    fig = plt.figure(figsize=(12,19))
    ax1 = plt.subplot2grid((2,7), (0,0), frameon=False)
    ax2 = plt.subplot2grid((2,7), (0,2), colspan=3)
    ax3 = plt.subplot2grid((2,7), (0,5))
    ax4 = plt.subplot2grid((2,7), (0,6))
    ax3_colbar = plt.subplot2grid((50, 7), (24,6))
    ax5 = plt.subplot2grid((27,2), (16,0), rowspan=4, frameon=False)
    ax6 = plt.subplot2grid((27,2), (16,1), rowspan=4, frameon=False)
    ax7 = plt.subplot2grid((27,2), (21,0), rowspan=4, frameon=False)
    ax8 = plt.subplot2grid((27,2), (21,1), rowspan=4, frameon=False)
    ft_full = pd.read_csv(ft_full, header=0, index_col=0)
    ft_full = ft_full.div(ft_full.sum(axis=0)).multiply(100) #convert the full feature table to % relative abundance
    #ft = ft*100 #convert the feature table to % relative abundance
    colnames = list(ft.columns) #get the column names (sample names)
    env, source, env_source_dict, env_source = [], [], {}, []
    for a in range(len(colnames)): #use the sample names to get the metadata categories that these samples belong to (environment and source)
        env.append(meta_dict[colnames[a]][3])
        source.append(meta_dict[colnames[a]][10])
        if meta_dict[colnames[a]][5] == 'lab':
            meta_dict[colnames[a]][5] = 'laboratory'
        env_source_dict[colnames[a]] = meta_dict[colnames[a]][3]+' '+'\n'+meta_dict[colnames[a]][10]
        env_source.append(meta_dict[colnames[a]][3]+' '+'\n'+meta_dict[colnames[a]][10])
    envs, sources = ['marine', 'freshwater', 'aquatic', 'terrestrial'], ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic']
    colors = [color_source[sources[0]], color_source[sources[1]], color_source[sources[2]], color_source[sources[3]], color_source[sources[4]]]
    venn_ax = [ax5, ax6, ax7, ax8]
    #sort the data and make venn diagrams with ASVs that overlap between treatments
    for a in range(len(envs)):
        asv_present, all_unique_asv, plotting = [], [], []
        names = []
        for b in range(len(sources)):
            keeping = []
            for c in range(len(env)):
                if env[c] == envs[a] and source[c] == sources[b]:
                    keeping.append(True)
                else:
                    keeping.append(False)
            source_ft = ft.loc[:, keeping] #this should now only keep the samples that belong to envs[a] and sources[b]
            ma = list(source_ft.max(axis=1))
            all_asv = list(source_ft.index.values) #get all ASVs
            asvs_here = []
            #only keep the ASVs if they are present (i.e. have abundance above 0)
            for d in range(len(all_asv)):
                if ma[d] == 'nan':
                    continue
                if ma[d] > 0:
                    asvs_here.append(all_asv[d])
            asv_present.append(asvs_here)
            all_unique_asv += asvs_here
            #tell it to only plot if this treatment exists
            if len(asvs_here) == 0:
                names.append(sources[b].capitalize())
                plotting.append(False)
            else:
                names.append(sources[b].capitalize())
                plotting.append(True)
        all_unique_asv = set(all_unique_asv) #get all ASVs that are present in this environment
        labels = get_venn_labels(asv_present) #get the labels based on which ASVs are present in each treatment (source)
        get_venn(labels, venn_ax[a], colors, names=names, plotting=plotting) #get the venn diagram for this environment
        venn_ax[a].set_title(envs[a].capitalize()+'\nTotal ASVs = '+str(len(all_unique_asv)), fontweight="bold", fontsize=fs_main) #change the title for this axes
        plt.sca(venn_ax[a])
        plt.xticks([]), plt.yticks([]) #remove the x and y ticks
    #get simpsons diversity for each sample
    indiv_div = ['Simpsons']
    for cn in colnames:
        indiv_div.append(get_diversity('Simpsons', list(ft_full.loc[:, cn])))
    #rename all sample names by a combination of the environment and source
    ft_full.rename(columns=env_source_dict, inplace=True)
    env_source = list(set(env_source)) #now get a list of all unique environments and sources
    ft_full = ft_full.reset_index()
    ft_full.rename(columns={'OTUID':'index'}, inplace=True)
    ft_full = ft_full.append(pd.Series(indiv_div, index=ft_full.columns), ignore_index=True) #add the simpsons diversity to the overall dataframe
    ft_full.index = ft_full['index']
    ft_full.drop('index', axis=1, inplace=True)
    ft_dendro = ft_full.drop('Simpsons', axis=0) #make a new dataframe without the simpsons index of diversity
    ft_dendro = ft_dendro.groupby(by=ft_dendro.columns, axis=1).mean() #group this by the environment/source name, taking a mean rather than a sum for each ASV
    ft_dendro.to_csv(basedir+'/grouped_samples_for_unifrac.csv')
    
    r.r_unifrac('agglom/reduced_tree_nr.tree') #run the R function to perform unifrac on these groupings
    
    w_uf = pd.read_csv(basedir+'weighted_unifrac_grouped_samples.csv', header=0, index_col=0)
    snames = list(w_uf.columns)
    r_name_dict = {}
    for s in snames:
        r_name_dict[s] = s.replace('..', ' \n').replace('.', ' ')
    w_uf.rename(columns=r_name_dict, index=r_name_dict, inplace=True)
    plt.sca(ax1)
    #ft_dendro_T = ft_dendro.transpose() #transpose the dataframe
    Z = ssd.squareform(w_uf)
    Z = hierarchy.linkage(Z, "ward")
    mpl.rcParams['lines.linewidth'] = 2
    hierarchy.set_link_color_palette(['k'])
    dn = hierarchy.dendrogram(Z, above_threshold_color='k', orientation='left') #plot this dendrogram of sample groupings
    y_labels, locs, xlocs, labels = list(ax1.get_yticklabels()), list(ax1.get_yticks()), list(ax1.get_xticks()), [] #get all labels so that we can work out which order samples are plotted in (and apply this to the other plots)
    for y in y_labels:
        labels.append(y.get_text())
    plot_labels, pl, colors = [], [], []
    grouping = list(w_uf.columns)
    sources.append('blank')
    #get the appropriate color and label for the order of samples
    for a in range(len(labels)):
        for b in range(len(envs)):
            if envs[b] in grouping[int(labels[a])-1]:
                colors.append(color_env[envs[b]])
        plot_labels.append(grouping[int(labels[a])-1].capitalize())
        pl.append(grouping[int(labels[a])-1])
    #add the yticks (with color corresponding to the environment that they came from)
    plt.yticks([])
    for a in range(len(locs)):
        ax1.text(xlocs[0]-1.5, locs[a], plot_labels[a], fontsize=fs_main, bbox=dict(facecolor=colors[a], alpha=0.4, pad=0.5, edgecolor='w'), ha='center', va='center')
    plt.sca(ax1)
    plt.xticks([])
    asv = list(ft_dendro.index.values)
    phyla, clss, spec = {}, {}, {}
    #for all ASVs, get the phylum, class and species that they belong to (although we don't actually use the species information)
    for a in range(len(asv)):
        tax = tax_dict[asv[a]]
        phyla[asv[a]] = tax[1]
        clss[asv[a]] = tax[3]
        spec[asv[a]] = tax[-1]
    #make new feature tables where ASVs are named by phyla/class and then group all of these phylums/classes (summing rather than taking means now)
    ft_T = ft_dendro.transpose()
    ft_phyla = ft_T.rename(columns=phyla)
    ft_phyla = ft_phyla.groupby(by=ft_phyla.columns, axis=1).sum()
    ft_class = ft_T.rename(columns=clss)
    ft_class = ft_class.groupby(by=ft_class.columns, axis=1).sum()
    ft_phyla = ft_phyla[ft_phyla.columns[ft_phyla.max() > 1]] #only keep those phyla that have a maximum abundance above 1%
    yt = []
    plot_names = ['Alteromonadales', 'Bacillales', 'Flavobacteriales', 'Oceanospirillales', 'Rhodobacterales', 'Rhodospirillales', 'Sphingomonadales', 'Vibrionales']
    colormap = mpl.cm.get_cmap('hot', 256)
    norm = mpl.colors.Normalize(vmin=0, vmax=15)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
    for a in range(len(pl)): #for every environment/source combination
        lst = list(ft_phyla.loc[pl[a], :]) #get a list of the relative abundance values for each phylum
        left = 0
        colors = get_cols(len(lst)) #get the colors for plotting (based on number of phyla above 1%)
        if len(lst) > 40: #if we have more than 40 phyla, print a warning
            print('Warning: there are not enough distinct colors for the stacked bar of the dendro and venn plot')
        for b in range(len(lst)): #plot the stacked bar chart for each phyla
            ax2.barh(a, lst[b], left=left, height=0.8, edgecolor='k', color=colors[b])
            left += lst[b]
        yt.append(a)
        ax2.barh(a, 100-left, left=left, height=0.8, edgecolor='k', color='k') #make the sample up to 100% (accounting for those phyla removed due to being <1% either above here or in the R script that werent agglomerated due to low abundance)
        left = 0
        x3 = []
        for c in range(len(plot_names)): #for each of the families that we are specifically interested in (due to them having been reported in a large number of Plastisphere studies)
            num = ft_class.loc[pl[a], plot_names[c]]
            ax3.barh(a, 1, left=left, height=1, edgecolor='k', color=m.to_rgba(num)) #now make the heatmap square, and scaled to this abundance
            x3.append(left+0.5) #add the y value
            left += 1
    #make the legend for the stacked bar of phyla, with the correct colors for each phylum
    handles = []
    phy = list(ft_phyla.columns)
    for a in range(len(phy)):
        handles.append(mlines.Line2D([], [], color=colors[a], marker='o', markersize=8, markeredgecolor='k', label=phy[a], linestyle=' '))
    handles.append(mlines.Line2D([], [], color='k', marker='o', markersize=8, markeredgecolor='k', label='Other', linestyle=' '))
    ax2.legend(handles=handles, bbox_to_anchor=(1.03, -0.06), ncol=round(len(lst)/4), fontsize=fs_small)
    #change the x and y limits and ticks for these plots
    ax2.set_ylim([-0.5, len(pl)-0.5])
    ax3.set_ylim([-0.5, len(pl)-0.5])
    ax2.set_xlim([0, 100])
    ax3.set_xlim([0,8])
    plt.sca(ax2)
    plt.yticks(yt, [])
    plt.sca(ax3)
    #add the ticks for each of the families of interest
    for a in range(len(plot_names)):
        plot_names[a] = r'$'+plot_names[a]+'$'
    plt.xticks(x3, plot_names, rotation=90, fontsize=fs_main)
    plt.yticks(yt, [])
    #now get only the simpsons diversity index for each (kept separately because we don't want to group it by sample as we want to be able to plot the outliers etc in the boxplot)
    ft_simps = ft_full.loc[['Simpsons'], :]
    for a in range(len(pl)): #for each sample type, get all samples and then make a box plot
        this_simp = ft_simps.loc[:, [pl[a]]]
        ax4.boxplot(this_simp, positions=[a], vert=False, showfliers=False)
    #add x ticks etc and some titles to our plots
    plt.sca(ax4)
    plt.yticks(yt, [])
    ax2.set_xlabel('Relative abundance (%)', fontsize=fs_main)
    ax2.set_title('Mean relative abundance of phyla', fontweight='bold', fontsize=fs_title-2)
    ax3.set_title('Mean relative\nabundance of\nkey orders', fontweight='bold', fontsize=fs_title-2)
    ax4.set_title('Simpsons\nindex\nof diversity', fontweight='bold', fontsize=fs_title-2)
    cb1 = mpl.colorbar.ColorbarBase(ax3_colbar, cmap=colormap, norm=norm, orientation='horizontal') #get the colorbar for the heatmap
    cb1.set_label('Relative abundance (%)', fontsize=fs_main)
    cb1.set_ticks([0, 5, 10, 15])
    cb1.ax.tick_params(labelsize=fs_small)
    plt.savefig(basedir+'/figures/Fig4_dendro_venn_'+str_norm+ext, dpi=dpi, bbox_inches='tight') #save the figure
    plt.close()
    return
    
def metacoder_py(ft, tax_dict, meta_dict, basedir, norm='rare'):
    '''
    This will separate the feature table to environments and perform pairwise comparisons between the same plastic types at early and late incubation times
    as well as different samples at the same incubation times.
    It takes as input:
        - ft (a feature table containing the abundance of ASVs as rows and samples as columns)
        - tax_dict (a dictionary containing ASV names as keys and taxonomy as values)
        - meta_dict (a dictionary containing all metadata, with sample names as keys)
        - basedir (the base directory to use for saving figures and files to)
    '''
    env, source, inc_time, source_inc, source_inc_dict, source_dict = [], [], [], [], {}, {} #set up some lists
    envs = ['marine', 'freshwater', 'aquatic', 'terrestrial'] #add list of the environments that we want to plot
    samples = list(ft.columns) #get a list of sample names
    for a in range(len(samples)): #for each of these sample names
        env.append(meta_dict[samples[a]][3]) #get the environment this sample was in
        source.append(meta_dict[samples[a]][10]) #get the plastic type/control type of the sample
        source_dict[samples[a]] = meta_dict[samples[a]][10] #add this to a new dictionary containing only this information
        inc_time.append(meta_dict[samples[a]][13]) #get the incubation time (early/late/collection)
        source_inc.append(meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13]) #get the name of the sample type and incubation time together
        source_inc_dict[samples[a]] = meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13] #and also add this to a new dictionary
    comparisons1 = [['aliphatic early', 'aliphatic late'], ['other plastic early', 'other plastic late'], ['biofilm early', 'biofilm late'],
                   ['aliphatic early', 'other plastic early'], ['aliphatic early', 'unknown plastic early'], ['aliphatic early', 'biofilm early'],
                   ['aliphatic late', 'other plastic late'], ['aliphatic late', 'unknown plastic late'], ['aliphatic late', 'biofilm late'],
                   ['other plastic early', 'unknown plastic early'], ['other plastic early', 'biofilm early'],
                   ['other plastic late', 'unknown plastic late'], ['other plastic late', 'biofilm late'],
                   ['unknown plastic early', 'biofilm early'], ['unknown plastic late', 'biofilm late']] #set up a list of all comparisons we will make
    comparisons2 = [['aliphatic', 'other plastic'], ['aliphatic', 'unknown plastic'], ['aliphatic', 'biofilm'], ['aliphatic', 'planktonic'], 
                    ['other plastic', 'unknown plastic'], ['other plastic', 'biofilm'], ['other plastic', 'planktonic'], 
                    ['unknown plastic', 'biofilm'], ['unknown plastic', 'planktonic'], 
                    ['biofilm', 'planktonic']] #and also a second lot of comparisons
    comps = [comparisons1, comparisons2] #add these two sets of comparisons to a list
    asvs = list(ft.index.values) #get a list of all ASVs
    lineage = [] #and an empty list - this will contain the full lineage of each taxon, as required by the R metacoder package for these plots
    for a in range(len(asvs)): #for each asv
        lineage_single = tax_dict[asvs[a]] #get all taxonomic levels for this ASV
        if lineage_single[0] == 'Archaea': #if it is archaea
            ft.drop(asvs[a], axis=0, inplace=True) #remove it from the feature table (there were too few archaea for useful comparisons to be made)
            continue #and continue to the next ASV
        string = '' #add an empty string
        for b in range(len(lineage_single)): #for each level in the taxonomy
            length = 7
            if b > 0 and lineage_single[b] == lineage_single[b-1]: #if this level is the same as the previous level (i.e. it is unclassified at lower taxonomic levels)
                length = b #then this is as many levels as we will want to add
        for b in range(length): #for each of these levels
            if lineage_single[b] != '': #if they're not empty strings
                if string != '': #and we have already add a level to our overall string
                    string += ';' #add a semi-colon to separate levels
                string += lineage_single[b] #and add this level to the string
        lineage.append(string) #and then add it to our overall list
    ft.insert(0, "lineage", lineage, True) #insert this list at the beginning of our feature table
    for a in range(len(envs)): #for each environment
        for z in range(len(comps)): #and each set of comparisons
            comparisons = comps[z] #set the competition we are using
            if comps[z] == comparisons1: #if its the first set of competitions
                treat = source_inc #we're renaming all of our samples using the source and incubation time
                treat_dict = source_inc_dict #and using the dictionary with source and incubation time
            else: #otherwise
                treat = source #we're only renaming them based on the source
                treat_dict = source_dict #and using the source dictionary
            for b in range(len(comparisons)): #for each of our comparisons
                if os.path.exists(basedir+'/figures/metacoder/'+norm+'/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'_labels.pdf'):
                  if os.path.exists(basedir+'/figures/metacoder/'+norm+'/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.pdf'):
                    print('Already have '+basedir+'/figures/metacoder/'+norm+'/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.pdf')
                    continue
                if len(comparisons[b]) < 3: #if we haven't already
                    comparisons[b].append('lineage') #add the 'lineage' to our list
                keeping = [True] #and tell it that we're keeping the first column (with out lineage information)
                for c in range(len(env)): #for each name in our list of environments
                    if envs[a] == env[c] and treat[c] in comparisons[b]: #if it is the same as the environment that we're currently looking at, and the name of the sample matches the current comparison
                        keeping.append(True) #then we'll be keeping this column
                    else: #otherwise
                        keeping.append(False) #we won't be keeping this column
                env_comp_ft = ft.loc[:, keeping] #now make a new feature table with only those columns that correspond to samples in this comparison
                env_comp_ft.rename(columns=treat_dict, inplace=True) #rename all of the samples to just the general sample type
                names = list(set(env_comp_ft.columns)) #now get a list of those names
                if names == ['lineage'] or len(names) == 2: #if we didn't keep any samples, or we only have one sample and the lineage column
                    continue #then carry on to the next comparison as we don't have enough samples for this
                names_dict = {} #make a new name dictionary
                for d in range(len(names)): #for each name in our list
                    if names[d] == comparisons[b][0]: #if it's equal to the first name in our comparison list
                        names_dict[names[d]] = 'Treat1' #then rename it Treat1
                        names[d] = names[d].replace(' early', '') #and in the list (not in the column) remove the 'early' so we can get the right color for this sample type
                        names[d] = names[d].replace(' late', '') #same for late
                        col1 = color_source[names[d]] #get the color
                    elif names[d] == comparisons[b][1]: #otherwise, if it's equal to the second name in our comparison list
                        names_dict[names[d]] = 'Treat2' #rename it Treat2
                        names[d] = names[d].replace(' early', '') #and again remove the 'early'
                        names[d] = names[d].replace(' late', '') #or late
                        col2 = color_source[names[d]] #and get the color
                if comparisons == comparisons1 and b < 3: #if its one of the comparisons between early and late samples in the same sample type
                    col1, col2 = color_source['early'], color_source['late'] #then get the early and late colors
                env_comp_ft.rename(columns=names_dict, inplace=True) #rename the columns with the dictionary we have created
                env_comp_ft.to_csv(basedir+'/metacoder/'+norm+'/metacoder_'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.csv') #save the .csv file with the name of this comparison, incase we need to look in the future
                env_comp_ft.to_csv(basedir+'metacoder.csv') #save the .csv file with the information for running this comparison in metacoder
                colors_metacoder = [col1, col2]
                r.suppressWarnings(r.metacoder(colors_metacoder, norm))
                os.rename(basedir+'metacoder.pdf', basedir+'/figures/metacoder/'+norm+'/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.pdf') #rename the first file, that doesn't have labels on the metacoder plot
                os.rename(basedir+'metacoder_labels.pdf', basedir+'/figures/metacoder/'+norm+'/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'_labels.pdf') #and then rename the second one, that does have labels on the metacoders plot
    return
    
def annotate_heatmap(ax, df, cmap='inferno', yticks=True, xticks=True, rnd=1, annotate=True, italics=False, vmax=False, vmin=False, annotate_only=False, norm=None):
    '''
    Function to plot a heatmap and to optionally annotate each cell of the heatmap
    Takes as input:
        - ax (the axis to plot the heatmap on)
        - df (a dataframe containing the data to plot as floats/integers in the cells. Y labels will be the index values and X labels the column headers)
        - cmap (string of which colormap to use for plotting, default is 'inferno')
        - yticks (optional boolean or string for whether to have ytick labels. Can be True (default), False or 'right' if you want the yticks to show on the right of the plot)
        - xticks (optional boolean or string for whether to have xtick labels. Can be True (default), False or 'top' if you want the xticks to show on the top of the plot)
        - rnd (number of decimal places to round the numbers to when the heatmap is annotated, default is 1)
        - annotate (boolean that decides whether to annotate the heatmap or not. Default is True)
        - italics (boolean of whether to make y values italic or not - i.e. if they are species names)
        - vmax (maximum value to scale color values to in the heatmap, default is False but this should be a number to work)
        - annotate_only (boolean as to whether to only annotate an existing heatmap)
    '''
    #get a heatmap plot using the dataframe df on the axis ax, with the colormap cmap, 
    #with additional options for whether we are to annotate the boxes of the heatmap, 
    #add x or y tick labels and how many decimal places to round each text label to
    plt.sca(ax) #set the axis 
    maxs = [] #get the maximum value of the dataframe
    for a in  df.columns:
        maxs.append(df[a].max())
    ma = max(maxs)
    if annotate: #if we are adding text labels
        for a in range(len(list(df.columns))): #go through each column
            for b in range(len(list(df.index.values))): #and now each row of that column
                num = round(df.iloc[b, a], rnd) #get the value corresponding to the row and columns
                if ma > num:
                    if num < (0.6)*ma: col = 'w' #if the value is below 60% of the maximum value, make the text color white
                    else: col = 'k' #otherwise, make the text color black
                    if cmap == 'PuBu': #if we are using the colormap 'PuBu', do the opposite in terms of text color
                        if num < (0.6)*ma: col = 'k'
                        else: col = 'w'
                else:
                    col = 'k'
                    if num == ma and cmap == 'PuBu': col = 'w'
                if cmap == 'Purples' or cmap == 'Blues':
                    if num > 0.9: col = 'w'
                    else: col = 'k'
                plt.text(a+0.5, b+0.5, str(num), ha='center', va='center', color=col, fontsize=fs_small) #and plot this text label
    if annotate_only:
        return
    if norm == 'clr' or norm == 'log':
      if vmax != False and vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax)) #make the heatmap
      elif vmax != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.LogNorm(vmin=0, vmax=vmax)) #make the heatmap
      elif vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.LogNorm(vmin=vmin)) #make the heatmap
      elif vmin == None:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.LogNorm())
      else:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.LogNorm(vmin=0)) #make the heatmap
    else:
      if vmax != False and vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=vmin, vmax=vmax) #make the heatmap
      elif vmax != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0, vmax=vmax) #make the heatmap
      elif vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=vmin) #make the heatmap
      elif vmin == None:
          plt.pcolor(df, edgecolor='k', cmap=cmap)
      else:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0) #make the heatmap
    rownames = list(df.index.values)
    for a in range(len(rownames)):
        if len(rownames[a]) > 35:
            row = rownames[a].split(' ')
            if len(row) == 1:
                row = rownames[a].split('-')
            string = ''
            for b in range(len(row)):
                if b > 0 and b < 2:
                    string += '\n'+row[b]+' '
                else:
                    string += row[b]+' '
            rownames[a] = string
    if yticks == 'right':
        plt.yticks(np.arange(0.5, len(df.index), 1), rownames, fontsize=fs_small, linespacing=0.8) #add them based in the row values of the dataframe
        ax.yaxis.tick_right()
    elif yticks: #if we are adding y tick labels
        plt.yticks(np.arange(0.5, len(df.index), 1), rownames, fontsize=fs_main, linespacing=0.8) #add them based in the row values of the dataframe
    else:
        plt.yticks([]) #otherwise, add none
    if xticks == 'top':
        colnames = list(df.columns)
        for col in range(len(colnames)):
            ind = col
            col = colnames[col].split(' ')
            string = ''
            for a in range(len(col)):
                if a > 0:
                    string += '\n'+col[a]
                else:
                    string += col[a]+' '
            colnames[ind] = string
        plt.xticks(np.arange(0.5, len(df.columns), 1), colnames, fontsize=fs_main)
        ax.xaxis.tick_top()
    elif xticks: #if we are adding x tick labels (do the same as for the y tick labels)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90, fontsize=fs_main)
    else:
        plt.xticks([])
    return
    
def get_rf_comparison(basedir):
    norm = ['rare', 'rel_abun', 'log', 'clr']
    levels = ['phyla', 'classes', 'orders', 'families', 'genera', 'species']  
    order_cols, order_rows = [], []
    fig = plt.figure(figsize=(13, 10)) 
    for l in range(len(levels)):
        level = levels[l]
        dfs, dfs_tax = [], []
        for n in norm:
            df = pd.read_csv(basedir+'random_forest/'+n+'/'+level+'_overall.csv', header=0, index_col=0)
            try:
                df_tax = df.drop(['Score', 'OOB_score'], axis=0)
            except:
                df_tax = df.drop(['Score'], axis=0)
            df_tax['Mean'] = df_tax.mean(axis=1)
            df_tax = pd.DataFrame(df_tax.loc[:, 'Mean'])
            df = pd.DataFrame(df.loc['Score', :])
            df = df.transpose()
            df.rename(index={'Score':n}, inplace=True)
            dfs.append(df)
            df_tax = df_tax.transpose()
            df_tax.rename(index={'Mean':n}, inplace=True)
            dfs_tax.append(df_tax)
        dfs = pd.concat(dfs)*100
        dfs_tax = pd.concat(dfs_tax).fillna(value=0)
        dfs_tax = dfs_tax.groupby(by=dfs_tax.index.values, axis=0).sum()
        for a in range(2):
            dfs['Mean'] = dfs.mean(axis=1)
            dfs = dfs.sort_values(by='Mean', axis=0, ascending=True)
            if a == 1: 
                dfs = dfs.sort_values(by='Mean', axis=0, ascending=False)
            else:
                df_mean = pd.DataFrame(dfs.loc[:, 'Mean'])
            dfs.drop('Mean', axis=1, inplace=True)
            dfs = dfs.transpose()
        dfs.rename(index=norm_names, columns=rename_plots, inplace=True)
        df_mean.rename(index=norm_names, inplace=True)
        order_rows = ['Rarefied', 'Relative\nabundance', 'Log', 'CLR']
        if order_cols == []:
            order_cols = list(dfs.columns)
        else:
            dfs = dfs.loc[:, order_cols]
        dfs = dfs.loc[order_rows, :]
        df_mean = df_mean.loc[order_rows, :]
        ax = plt.subplot2grid((6,120), (l, 5), colspan=76)
        ax_mean = plt.subplot2grid((6, 120), (l, 0), colspan=4)
        ax_con = plt.subplot2grid((6,6), (l, 4))
        ax_con_tax = plt.subplot2grid((6,6), (l, 5))
        xtcks = False
        if l == 5: xtcks = True
        annotate_heatmap(ax, dfs, xticks=xtcks, yticks=False)
        
        annotate_heatmap(ax_mean, df_mean, xticks=xtcks)
        ax_mean.set_ylabel(level.capitalize(), fontsize=fs_main, fontweight='bold')
        concs = []
        order_rows = ['rare', 'rel_abun', 'log', 'clr']
        dfs_tax = dfs_tax.loc[order_rows, :]
        for a in range(4):
            conc = []
            for b in range(4):
                l1 = dfs.iloc[a, :].values
                l2 = dfs.iloc[b, :].values
                conc.append(concordance_index(l1, l2))
            concs.append(conc)
        concs = pd.DataFrame(concs, index=dfs.index.values, columns=dfs.index.values)
        concs = concs.loc[:, ['CLR', 'Log', 'Relative\nabundance', 'Rarefied']]
        annotate_heatmap(ax_con, concs, cmap='Blues', rnd=2, yticks=False, xticks=xtcks, vmin=0.75)
        concs_tax = []
        for a in range(4):
            conc_tax = []
            for b in range(4):
                l1 = dfs_tax.iloc[a, :].values
                l2 = dfs_tax.iloc[b, :].values
                conc_tax.append(concordance_index(l1, l2))
            concs_tax.append(conc_tax)
        concs_tax = pd.DataFrame(concs_tax, index=dfs_tax.index.values, columns=dfs_tax.index.values)
        concs_tax.rename(columns=norm_names, inplace=True)
        concs_tax = concs_tax.loc[:, ['CLR', 'Log', 'Relative\nabundance', 'Rarefied']]
        annotate_heatmap(ax_con_tax, concs_tax, cmap='Purples', rnd=2, yticks=False, xticks=xtcks, vmin=0.75)
        
        if l == 0:
            ax.set_title('Classification accuracy (%)', fontsize=fs_title, fontweight='bold')
            ax_con.set_title('Concordance in\nclassification accuracy', fontsize=fs_main, fontweight='bold')
            ax_con_tax.set_title('Concordance in\nfeature importance', fontsize=fs_main, fontweight='bold')
    
    plt.savefig(basedir+'/figures/FigS2_RF_compare'+ext, dpi=600, bbox_inches='tight')
    plt.close()
    return 
    
def generate_rf(X_train, y_train, X_test, y_test, rc='cls', est=10000, n_jobs=None): #generate a random forest based on the data being split to train and test
    '''
    Function to generate a random forest model for one dataset, depending on whether we want to use the classification or regression model
    It takes as input:
        - X_train, y_train, X_test, y_test (this should be the output from the sklearn train_test_split function)
        - rc (a string of either 'cls' or 'reg', depending on whether we want to perform classification or regression, respectively. Default='cls' as this will work on numeric or non-numeric data, where regression requires numeric)
        - est (the number of estimators to use. In general, the higher this number the better/more accurate the model can fit. Default=10000)
        - n_jobs (the number of processors to use. Default=None, meaning 1 is used)
    Returns:
        - rf (the random forest model)
        - rf.score (the classification accuracy of this model)
        - rf.feature_importances_ (the importance values for each included feature/taxa)
    '''
    if rc == 'cls': #if we are using a classification (i.e. discrete categories)
        rf = RandomForestClassifier(n_estimators=est, min_samples_leaf=3, n_jobs=n_jobs, random_state=seed, oob_score=False)
    else: #if we are using a regression (i.e. continuous categories)
        rf = RandomForestRegressor(n_estimators=est, min_samples_leaf=3, n_jobs=n_jobs, random_state=seed, oob_score=False)
    rf.fit(X_train, y_train) #fit out data to either the regressor or classifier
    return rf, rf.score(X_test, y_test), rf.feature_importances_ #return the forest (i.e. features), score (how well it classifies the test data) and feature importances (ASV importances)

def get_single_forest(ft, meta_df, sn, level, est=10000, n_jobs=None):
    '''
    Function to calculate random forest models for each metadata category for a feature table at one phylogenetic level
    Takes as input:
        - ft (feature table/dataframe with sample names as columns and ASVs/taxa names as rows)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - sn (name/string to use to save the results of these random forests)
        - level (phylogenetic level that these comparisons are being performed at, where 1 is phylum and 7 is ASV)
        - est (integer, number of estimators to use to calculate the random forests. Default=10000)
        - n_jobs (integer, number of processors to use for calculations. Default is None, meaning that 1 will be used)
    Returns:
        - scores (scores for how well the random forest models performed for each metadata category)
        - dfs (dataframe containing the importance values for each taxon across all metadata categories)
    '''
    meta_df = meta_df.loc[list(ft.index.values), :]
    cls_reg = {'Study':'cls', 'Latitude':'reg', 'Longitude':'reg', 'Environment':'cls', 'WaterOrSediment':'cls',
       'LabOrField':'cls', 'IncubationOrCollection':'cls', 'Source':'cls', 'MaterialType':'cls',
       'PlasticTypeSpecific':'cls', 'PlasticTypeGeneral':'cls', 'DEPTH':'reg', 'IncubationTime':'reg',
       'IncubationGeneral':'cls', 'Temperature':'reg', 'Salinity':'reg', 'Light':'cls', 'Season':'cls',
       'PrimerPair':'cls', 'DNAExtraction':'cls', 'CollectionDate':'cls', 'PlasticOnly':'cls'} #dictionary telling it whether each category should use a classification (cls) or regression (reg) algorithm for the random forests
    ft_in = ft #save the original feature table
    ft = ft.rename_axis('ID').values #rename the feature table axis (i.e. remove the indice and column names)
    max_abs_scaler = preprocessing.MaxAbsScaler() #scale all ASVs to the maximum absolute value
    ft_scale = max_abs_scaler.fit_transform(ft) #apply this to our feature table
    meta_split, reg_cls, col_names = [], [], []
    for col in meta_df.columns: #make lists of all of the sample names based on the sample type within each meta category
        meta_split.append(meta_df[col])
        reg_cls.append(cls_reg[col])
        col_names.append(col)
    scores, importances = [], []
    count = 0
    for a in range(len(meta_split)): #for each meta-category
        meta = meta_split[a].rename_axis('ID').values #split the metadata dataframe to be only the category we are currently looking at
        if len(set(meta)) == 2 and 'unknown' in meta:
            continue
        keeping = []
        for z in range(len(meta)):
            if meta[z] == 'unknown':
                keeping.append(False)
            else:
                keeping.append(True)
        meta_reduced = meta[keeping]
        ft_scale_reduced = ft_scale[keeping]
        X_train, X_test, y_train, y_test = train_test_split(ft_scale_reduced, meta_reduced, test_size=0.2) #split the data to training and test data (with 80% of samples being used for training and 20% for testing)
        RF, RF_score, RF_importances = generate_rf(X_train, y_train, X_test, y_test, rc=reg_cls[a], est=est, n_jobs=n_jobs) #generate the random forest for this meta category
        #print(col_names[a], RF_score) #print the category and random forest score/classification accuracy
        scores.append(RF_score), importances.append(RF_importances) #add the scores and importances for this category to the overall lists
        this_rf = {'ASV':ft_in.columns, 'Importance':RF_importances} #make a dataframe with the information on how important each ASV/taxa is
        df_rf = pd.DataFrame(this_rf, columns = ['ASV', 'Importance'])
        df_rf.sort_values(by=['Importance'], ascending=False, inplace=True) #sort the values by importance
        df_rf_score = df_rf.set_index('ASV')
        df_rf_score = df_rf_score
        df_rf = df_rf[df_rf.max(axis=1) > 0.00]
        df_rf = df_rf[:50] #now only keep the top 50 features
        df_rf = df_rf.set_index('ASV') #if we didn't use get_heatmap_random_forest, then we need to uncomment this line
        df_rf = df_rf_score.rename(columns={'Importance':col_names[a]}) #rename the column to be this meta category
        if count == 0:
            dfs = df_rf
            count += 1
        else:
            dfs = dfs.merge(right=df_rf, on='ASV')
    if scores == []:
        return [], []
    scores = pd.DataFrame([scores], index=['Score'], columns=list(dfs.columns))
    dfs = pd.concat([dfs, scores], sort=True)
    dfs.to_csv(sn+'.csv')
    return scores, dfs

def get_random_forests(ft, tax_dict, meta_df, basedir, est=10000, n_jobs=None, sn=None, norm=None):
    '''
    Function to call a function to get random forests for each phylogenetic level, across all metadata categories
    Takes as input:
        - ft (feature table/dataframe with sample names as columns and ASVs/taxa names as rows)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - basedir (name of the directory to save the files to)
        - est (integer, number of estimators to use to calculate the random forests. Default=10000)
        - n_jobs (integer, number of processors to use for calculations. Default is None, meaning that 1 will be used)
    Returns:
        - nothing, but saves a .csv file for each phylogenetic level, containing information on taxon importance at that level
    '''
    this_sn = str(sn)
    asv = list(ft.index.values)
    for a in [7]:#[1, 2, 3, 4, 5, 6, 7]:
        #if norm != None and a < 5: continue
        this_asv = list(asv)
        new_names, other_levels = {}, {}
        for b in range(len(this_asv)):
            this_tax = tax_dict[this_asv[b]]
            prev = ''
            for c in range(len(this_tax)):
                if this_tax[c] != '':
                    prev = this_tax[c]
                else:
                    this_tax[c] = prev
            if a < 7:
                new_names[this_asv[b]] = this_tax[a]
                this_asv[b] = this_tax[a]
            else:
                new_names[this_asv[b]] = this_asv[b]
            other_levels[this_asv[b]] = this_tax[:a]
        this_ft = ft.rename(index=new_names) #now rename the feature table using the new names that we have just defined
        this_ft = this_ft.groupby(by=this_ft.index, axis=0).sum() #group all of the names that are together, taking a sum for each phyla
        this_ft = this_ft.transpose()
        if sn == None:
            if norm == None:
              lev_sn = basedir+'/random_forest/'+phylo_level_names[a]+'_overall'
            else:
              lev_sn = basedir+'/random_forest/'+norm+'/'+phylo_level_names[a]+'_overall'
        else:
            lev_sn = this_sn+phylo_level_names[a]
        
        if not os.path.exists(lev_sn+'.csv'): #only do this if we don't already have this file
          get_single_forest(this_ft, meta_df, lev_sn, a, est=est, n_jobs=n_jobs) #get all of the individual forests for this phylogenetic level
        else:
          print('Didnt get this forest '+lev_sn+' because we already had a file saved under the save name') 
    return

def get_environment_random_forest(ft, tax_dict, meta_df, meta_dict, basedir, est=10000, n_jobs=1, norm=None):
    '''
    Function to calculate random forest models for specific metadata categories in different environments.
    It takes as input:
        - ft (feature table data frame containing taxa as rows and sample names as columns)
        - tax_dict (dictionary containing the full taxonomy for each ASV)
        - meta_df (dataframe containing all metadata information, with samples as rows and metadata categories as columns)
        - meta_dict (dictionary containing all metadata with sample names as keys)
        - basedir (directory to save .csv output files to)
        - est (number of estimators to use for calculating random forests, default is 10000. Note that lower numbers will lead to lower accuracy)
        - n_jobs (number of processors to use for calculating random forests, default is 1)
    It uses the other general random forest function but with filtered versions of the metadata and feature tables, and could easily be edited to do this separately for lab/field studies and to do more than just 'PlasticTypeGeneral'
    '''
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    labfield, labindex = [['lab', 'field']], 5
    
    all_samples = list(ft.columns)
    all_samples_meta = list(meta_df.index.values)
    for a in range(len(envs)):
        for b in range(len(labfield)):
            meta_names = list(meta_df.columns)
            keeping = []
            for c in range(len(all_samples)):
                if meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] == labfield[b]:
                    keeping.append(True)
                elif meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] in labfield[b]:
                    keeping.append(True)
                else:
                    keeping.append(False)
            single_ft = ft.loc[:, keeping]
            if len(single_ft.index.values) == 0:
                continue
            single_ft = single_ft[single_ft.max(axis=1) > 0]
            keeping_meta = []
            for d in range(len(all_samples_meta)):
                if meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] == labfield[b]:
                    keeping_meta.append(True)
                elif meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] in labfield[b]:
                    keeping_meta.append(True)
                else:
                    keeping_meta.append(False)
            single_meta_df = meta_df.loc[keeping_meta]
            meta_names.remove('PlasticTypeGeneral')
            single_meta_df.drop(meta_names, axis=1, inplace=True)
            if len(single_meta_df.index.values) == 0:
                continue
            if norm != None:
                sn = basedir+'/random_forest/'+norm+'/single_environment/'+envs[a]+'_'
            else:
                sn = basedir+'/random_forest/single_environment/'+envs[a]+'_'
            get_random_forests(single_ft, tax_dict, single_meta_df, basedir, est=est, n_jobs=n_jobs, sn=sn, norm=norm) 
    return
    
def get_ancom(ft_abun):
    '''
    Function to perform ANCOM tests for statistical significance on a given feature table
    Takes as input:
        - ft_abun (a feature table dataframe with sample names as rows and taxon names as columns)
    Returns:
        - significant (a list of taxon names corresponding to all taxa that were determined as being statistically significantly differentially abundant between treatments)
    '''
    other_ft = ft_abun.copy()
    other_ft = other_ft.replace(to_replace=0,value = 0.0001)
    other_ft = other_ft.fillna(value=0.0001)
    other_ft = other_ft.transpose()
    ancom_df, percentile_df = ancom(other_ft, pd.Series(list(other_ft.index.values), index=list(other_ft.index.values)), multiple_comparisons_correction='holm-bonferroni')
    ancom_results, percentile_results, tax_names = ancom_df.values.tolist(), percentile_df.values.tolist(), ancom_df.index.values#, list(percentile_df.columns.values)
    significant, medians, not_sig, not_sig_medians = [], [], [], []
    for d in range(len(ancom_results)):
        if ancom_results[d][1] == True:
            if percentile_results[d][2] == 0.0001 and percentile_results[d][7] == 0.0001:
                not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]])
                not_sig.append(tax_names[d])
                significant.append(tax_names[d])
                medians.append([percentile_results[d][2], percentile_results[d][7]])
            else:
                not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]])
                not_sig.append(tax_names[d])
    return significant

def make_env_rf_plot(ft, tax_dict, basedir, ASV_dict, meta_dict, mx=0.005, norm='rare'):
    '''
    Function to make heatmaps (of abundance and including whether they were significantly differentially abundant according to ANCOM) of all taxa at all phylogenetic levels with a feature importance value of above 0.01 in the random forest models constructed for separate environments
    Takes as input:
        - ft (feature table dataframe containing sample names as columns and ASV names as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - basedir (directory to use to open random forest information as well as to save the resulting figure to)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
    Returns:
        - nothing, but saves the figure environment_random_forest to the figures directory
    '''
    envs, ei = ['marine', 'aquatic', 'freshwater', 'terrestrial'], 3
    envs_rf_imp = []
    for a in range(len(envs)):
        if os.path.exists(basedir+'/random_forest/'+norm+'/single_environment/'+envs[a]+'_'+str(mx)+'.csv'):
            envs_rf_imp.append(pd.read_csv(basedir+'/random_forest/'+norm+'/single_environment/'+envs[a]+'_'+str(mx)+'.csv', header=0, index_col=0))
            continue
        elif os.path.exists(basedir+'/random_forest/'+norm+'/single_environment/'+envs[a]+'.csv'):
            rf_file = pd.read_csv(basedir+'/random_forest/'+norm+'/single_environment/'+envs[a]+'.csv', header=0, index_col=0)
            rf_file = rf_file[rf_file['Importance'] >= mx]
            envs_rf_imp.append(rf_file)
            continue
        samples = list(ft.columns)
        keeping, count, rename = [], 0, {}
        for b in range(len(samples)):
            if meta_dict[samples[b]][ei] == envs[a]:
                keeping.append(True)
                rename[samples[b]] = meta_dict[samples[b]][10]
                count += 1 
            else:
                keeping.append(False)
        env_ft = ft.loc[:, keeping]
        env_ft = env_ft[env_ft.max(axis=1) > 0]
        env_ft.rename(columns=rename, inplace=True)
        scores = []
        asv = list(env_ft.index.values)
        for b in [7, 6, 5, 4, 3, 2, 1]:
            lvl = phylo_level_names[b]
            rf = pd.read_csv(basedir+'random_forest/'+norm+'/single_environment/'+envs[a]+'_'+lvl+'.csv', header=0, index_col=0)
            score = rf.loc[['Score']].values
            scores.append(score[0][0])
            try:
              rf.drop(['Score', 'OOB_score'], axis=0, inplace=True)
            except:
              rf.drop(['Score'], axis=0, inplace=True)
            rf = rf[rf.max(axis=1) > mx]
            lvl_rename, other_levels = {}, {}
            for c in range(len(asv)):
                tax = tax_dict[asv[c]]
                for d in range(len(tax)):
                    if tax[d] == 'Alteromonas macleodii' and d < 6:
                        for e in range(d, len(tax)):
                            tax[e] = 'Marinimicrobia (SAR406 clade)'
                            break
                    if tax[d] == '':
                        if d < 7:
                            if tax[-1] == tax[5]:
                                tax[d] = tax[-1]
                            else:
                                tax[d] = tax[-1]
                        else:
                            tax[d] = tax[-1]
                lvl_rename[asv[c]] = tax[b]
                adding = tax[:b]
                while len(adding) < 7:
                    adding.append(adding[-1])
                if b == 7:
                    other_levels[asv[c]] = adding
                else:
                    other_levels[tax[b]] = adding
            if b < 7:
                env_ft_level = env_ft.rename(index=lvl_rename)
                env_ft_level = env_ft_level.groupby(by=env_ft_level.index, axis=0).sum()
            else:
                env_ft_level = env_ft
            if 'blank' in env_ft_level.columns:
                env_ft_level.drop('blank', axis=1, inplace=True)
            significant = get_ancom(env_ft_level)
            env_ft_level = env_ft_level.groupby(by=env_ft_level.columns, axis=1).mean()
            rows = [['Taxa', 'kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera', 'species']+list(env_ft_level.columns)+['Importance', 'ANCOM']]
            rf_id = list(rf.index.values)
            for e in range(len(rf_id)):
                if b == 7:
                    this_name = ASV_dict[rf_id[e]]+': '+tax_dict[rf_id[e]][-1]
                else:
                    this_name = rename_plots[lvl]+': '+rf_id[e]
                if b < 6 and rf_id[e] == 'Alteromonas macleodii':
                    print(rf_id[e], this_name)
                    rf_id[e] = 'Marinimicrobia (SAR406 clade)'
                    this_name = 'Marinimicrobia (SAR406 clade)'
                try:
                    this_row = [this_name]+other_levels[rf_id[e]]+list(env_ft_level.loc[rf_id[e], :].values)+list(rf.loc[rf_id[e], :].values)
                except:
                    print(rf_id[e], this_name, envs[a], lvl)
                if rf_id[e] in significant:
                    this_row.append('S')
                else:
                    this_row.append('NS')
                if b == 7:
                    rows.append(this_row)
                else:
                    rows.append(this_row)
            if len(this_row) < 2:
                continue
            new_df = pd.DataFrame(rows[1:], columns=rows[0])
            if b == 7:
                all_rows = new_df
            else:
                all_rows = all_rows.append(new_df, ignore_index=True)
        all_rows = all_rows.sort_values(['kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera', 'species'], ascending=False)
        all_rows = all_rows.set_index('Taxa')
        envs_rf_imp.append(all_rows)
        all_rows.to_csv(basedir+'/random_forest/'+norm+'/single_environment/'+envs[a]+'_'+str(mx)+'.csv')
    length, width = [], []
    [(length.append(i.shape[0]), width.append(i.shape[1]-9)) for i in envs_rf_imp]
    ma = max(length)
    fig = plt.figure(figsize=(25, int(ma/5)))
    start = 0
    for a in range(len(envs_rf_imp)):
        env_rf = envs_rf_imp[a]
        env_rf = env_rf.drop(['kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera', 'species'], axis=1)
        treats1, treats2 = list(env_rf.columns)[:-2]+['ANCOM'], list(env_rf.columns)[:-2]+['Importance']
        rf = env_rf.drop(treats1, axis=1)
        sig = env_rf.drop(treats2, axis=1)
        env_rf.drop(['Importance', 'ANCOM'], axis=1, inplace=True)
        ax = plt.subplot2grid((ma, 60), (0, start), colspan=env_rf.shape[1], rowspan=env_rf.shape[0])
        start += env_rf.shape[1]
        ax_imp = plt.subplot2grid((ma, 60), (0, start), rowspan=env_rf.shape[0])
        start += 10
        env_rf.rename(columns=rename_plots, inplace=True)
        if 'Other plastic' in env_rf.columns:
            if 'Unknown plastic' in env_rf.columns:
                env_rf = env_rf.reindex(['Aliphatic', 'Other plastic', 'Unknown plastic', 'Biofilm', 'Planktonic'], axis=1)
            else:
                env_rf = env_rf.reindex(['Aliphatic', 'Other plastic', 'Biofilm', 'Planktonic'], axis=1)
        annotate_heatmap(ax_imp, rf, cmap='inferno', yticks=False, xticks=True, rnd=2, annotate=True)
        env_rf = env_rf.div(env_rf.max(axis=1), axis=0)
        annotate_heatmap(ax, env_rf, cmap='viridis', yticks=False, xticks=True, annotate=False)
        nums = np.arange(0.5, env_rf.shape[0]+0.5, 1)
        tax = list(env_rf.index.values)
        empty = []
        for b in range(len(nums)):
            color = 'k'
            if list(sig.loc[tax[b], :])[0][0] == 'S':
                color = '#900C3F'
            if 'Allorhizobium' in tax[b]:
                    tax[b] = 'Allorhizobium'
            tax[b] = tax[b].replace('Incertae Sedis', '')
            ax_imp.text(1.5, nums[b], tax[b], fontsize=fs_main, color=color, va='center')
            empty.append('')
        plt.sca(ax_imp)
        plt.yticks(nums, empty)
        ax_imp.yaxis.tick_right()
        ax.set_title(envs[a].capitalize(), fontsize=fs_title, fontweight='bold')
    plt.savefig(basedir+'/figures/FigS5_environment_random_forest_'+str(mx)+'_'+norm+'.pdf', dpi=dpi, bbox_inches='tight')
    plt.savefig(basedir+'/figures/FigS5_environment_random_forest_'+str(mx)+'_'+norm+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return
    
def get_rf_comparison_env(basedir):
  norm = ['rare', 'rel_abun', 'log', 'clr']
  levels = ['phyla', 'classes', 'orders', 'families', 'genera', 'species']
  envs = ['marine', 'aquatic', 'freshwater', 'terrestrial']
  
  fig = plt.figure(figsize=(13,10))
  for l in range(len(levels)):
    level = levels[l]
    all_scores = []
    df_norm = [[], [], [], []]
    for n in norm:
      scores = []
      for e in range(len(envs)):
        env = envs[e]
        df = pd.read_csv(basedir+'random_forest/'+n+'/single_environment/'+env+'_'+level+'.csv', header=0, index_col=0)
        score = df.loc['Score', :].values[0]*100
        scores.append(score)
        df.rename(columns={'PlasticTypeGeneral':n}, inplace=True)
        df.drop(['Score'], axis=0, inplace=True)
        df_norm[e].append(df)
      all_scores.append(scores)
    for d in range(len(df_norm)):
      df = df_norm[d]
      df = pd.concat(df).fillna(value=0)
      df = df.groupby(by=df.index.values, axis=0).sum()
      df_norm[d] = df
    score_df = pd.DataFrame(all_scores, columns=envs, index=norm)
    score_df.rename(index=norm_names, columns=rename_plots, inplace=True)
    xtcks = False
    if l == 5: xtcks = True
    df_mean = pd.DataFrame(score_df.mean(axis=1).values, columns=['Mean'], index=norm)
    score_df.rename(index=norm_names, columns={'marine':'Marine', 'aquatic':'Aquatic', 'freshwater':'Freshwater', 'terrestrial':'Terrestrial'}, inplace=True)
    df_mean.rename(index=norm_names, inplace=True)
    ax = plt.subplot2grid((6,25), (l, 1), colspan=4)
    ax_mean = plt.subplot2grid((6, 25), (l, 0), colspan=1)
    annotate_heatmap(ax_mean, df_mean, xticks=xtcks)
    annotate_heatmap(ax, score_df, xticks=xtcks, yticks=False)
    if l == 0:
      ax.set_title('Classification\naccuracy (%)', fontsize=fs_title, fontweight='bold')
    place = [6, 11, 16, 21]
    ax_mean.set_ylabel(level.capitalize(), fontsize=fs_main, fontweight='bold')
    for e in range(len(envs)):
      env = envs[e]
      feat_conc = plt.subplot2grid((6,25), (l,place[e]), colspan=4)
      if l == 0:
        feat_conc.set_title('Concordance in\nfeature importance', fontsize=fs_main, fontweight='bold')
        feat_conc.text(0.5, 1.4, env.capitalize(), fontsize=fs_title, fontweight='bold', transform = feat_conc.transAxes, ha='center', va='center')
      concs = []
      for a in range(4):
          conc = []
          for b in range(4):
              l1 = df_norm[e].iloc[a, :].values
              l2 = df_norm[e].iloc[b, :].values
              conc.append(concordance_index(l1, l2))
          concs.append(conc)
      concs = pd.DataFrame(concs, index=df_norm[e].columns, columns=df_norm[e].columns)
      concs.rename(index=norm_names, columns=norm_names, inplace=True)
      concs = concs.loc[:, ['CLR', 'Log', 'Relative\nabundance', 'Rarefied']]
      annotate_heatmap(feat_conc, concs, cmap='Blues', rnd=2, yticks=False, xticks=xtcks, vmin=0.75)
  plt.savefig(basedir+'/figures/FigS3_RF_env_compare'+ext, dpi=600, bbox_inches='tight')
  return
  
def annotate_heatmap(ax, df, cmap='inferno', yticks=True, xticks=True, rnd=1, annotate=True, italics=False, vmax=False, vmin=False, annotate_only=False, norm=None):
    '''
    Function to plot a heatmap and to optionally annotate each cell of the heatmap
    Takes as input:
        - ax (the axis to plot the heatmap on)
        - df (a dataframe containing the data to plot as floats/integers in the cells. Y labels will be the index values and X labels the column headers)
        - cmap (string of which colormap to use for plotting, default is 'inferno')
        - yticks (optional boolean or string for whether to have ytick labels. Can be True (default), False or 'right' if you want the yticks to show on the right of the plot)
        - xticks (optional boolean or string for whether to have xtick labels. Can be True (default), False or 'top' if you want the xticks to show on the top of the plot)
        - rnd (number of decimal places to round the numbers to when the heatmap is annotated, default is 1)
        - annotate (boolean that decides whether to annotate the heatmap or not. Default is True)
        - italics (boolean of whether to make y values italic or not - i.e. if they are species names)
        - vmax (maximum value to scale color values to in the heatmap, default is False but this should be a number to work)
        - annotate_only (boolean as to whether to only annotate an existing heatmap)
    '''
    #get a heatmap plot using the dataframe df on the axis ax, with the colormap cmap, 
    #with additional options for whether we are to annotate the boxes of the heatmap, 
    #add x or y tick labels and how many decimal places to round each text label to
    plt.sca(ax) #set the axis 
    maxs = [] #get the maximum value of the dataframe
    for a in  df.columns:
        maxs.append(df[a].max())
    ma = max(maxs)
    if annotate: #if we are adding text labels
        for a in range(len(list(df.columns))): #go through each column
            for b in range(len(list(df.index.values))): #and now each row of that column
                num = round(df.iloc[b, a], rnd) #get the value corresponding to the row and columns
                if ma > num:
                    if num < (0.6)*ma: col = 'w' #if the value is below 60% of the maximum value, make the text color white
                    else: col = 'k' #otherwise, make the text color black
                    if cmap == 'PuBu': #if we are using the colormap 'PuBu', do the opposite in terms of text color
                        if num < (0.6)*ma: col = 'k'
                        else: col = 'w'
                else:
                    col = 'k'
                    if num == ma and cmap == 'PuBu': col = 'w'
                if cmap == 'Purples' or cmap == 'Blues':
                    if num > 0.9: col = 'w'
                    else: col = 'k'
                plt.text(a+0.5, b+0.5, str(num), ha='center', va='center', color=col, fontsize=fs_small) #and plot this text label
    if annotate_only:
        return
    if norm == 'clr' or norm == 'log':
      if vmax != False and vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.SymLogNorm(vmin=vmin, vmax=vmax, linthresh=0.03)) #make the heatmap
      elif vmax != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.SymLogNorm(vmin=0, vmax=vmax, linthresh=0.03)) #make the heatmap
      elif vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.SymLogNorm(vmin=vmin, linthresh=0.03)) #make the heatmap
      elif vmin == None:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.SymLogNorm(linthresh=0.03))
      else:
          plt.pcolor(df, edgecolor='k', cmap=cmap, norm=colors.SymLogNorm(vmin=0, linthresh=0.03)) #make the heatmap
    else:
      if vmax != False and vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=vmin, vmax=vmax) #make the heatmap
      elif vmax != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0, vmax=vmax) #make the heatmap
      elif vmin != False:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=vmin) #make the heatmap
      elif vmin == None:
          plt.pcolor(df, edgecolor='k', cmap=cmap)
      else:
          plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0) #make the heatmap
    rownames = list(df.index.values)
    for a in range(len(rownames)):
        if len(rownames[a]) > 35:
            row = rownames[a].split(' ')
            if len(row) == 1:
                row = rownames[a].split('-')
            string = ''
            for b in range(len(row)):
                if b > 0 and b < 2:
                    string += '\n'+row[b]+' '
                else:
                    string += row[b]+' '
            rownames[a] = string
    if yticks == 'right':
        plt.yticks(np.arange(0.5, len(df.index), 1), rownames, fontsize=fs_small, linespacing=0.8) #add them based in the row values of the dataframe
        ax.yaxis.tick_right()
    elif yticks: #if we are adding y tick labels
        plt.yticks(np.arange(0.5, len(df.index), 1), rownames, fontsize=fs_main, linespacing=0.8) #add them based in the row values of the dataframe
    else:
        plt.yticks([]) #otherwise, add none
    if xticks == 'top':
        colnames = list(df.columns)
        for col in range(len(colnames)):
            ind = col
            col = colnames[col].split(' ')
            string = ''
            for a in range(len(col)):
                if a > 0:
                    string += '\n'+col[a]
                else:
                    string += col[a]+' '
            colnames[ind] = string
        plt.xticks(np.arange(0.5, len(df.columns), 1), colnames, fontsize=fs_main)
        ax.xaxis.tick_top()
    elif xticks: #if we are adding x tick labels (do the same as for the y tick labels)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90, fontsize=fs_main)
    else:
        plt.xticks([])
    return

def sort_by_tax(df, level, other_levels, tax_dict=False, names=False):
    '''
    Function to sort a dataframe based on the higher-level taxonomy (i.e. if the dataframe contains Oceanospirillales, Alteromonadales, Vibrionales, Rhodobacterales and Rhizobiales, then these will be grouped alphabetically to Alpha (Rhizobiales & Rhodobacterales) and Gamma (Alteromonadales, Oceanospirillales & Vibrionales) proteobacteria)
    Takes as input:
        - df (dataframe containing taxa as rows and sample names as columns)
        - level (numeric value of the current level, where 1 is phylum and 7 is ASV)
        - other_levels (dictionary containing taxa at this level as the key and the other higher levels for each as a list)
        - tax_dict (the dictionary containing all taxa - this is False as default and only useful is the level is 7)
        - names (this is redundant but left in so as not to mess up other functions. It is False as default, but is over-written within the function)
    Returns:
        - df_tax (an ordered dataframe containing taxa as rows and sample names as columns)
        - names (a dictionary containing the current taxa as keys and the lowest classification also if the level we are working at is ASV)
    '''
    #function to sort the dataframe df (with index values being the taxon/ASV names) by the higher phylogeny
    current_tax = df.index.values #get the current values
    df_tax = df
    sort_levels, names = [], {}
    for a in range(level): #for each phylogenetic level that is higher
        new_col = []
        for b in range(len(current_tax)):
            if level < 6 and current_tax[b] == 'Alteromonas macleodii':
                current_tax[b] = 'Marinimicrobia (SAR406 clade)'
            if current_tax[b][0] == '_':
                current_tax[b] = current_tax[b].replace('_', '')
            new_col.append([other_levels[current_tax[b]][a], current_tax[b]]) #append the appropriate value for the current higher level as well as the current taxon name
            if tax_dict != False: #if tax_dict isn't False (i.e. if we are looking at the level of ASV - we use tax_dict as here we already worked out the lowest classification for each ASV)
                names[current_tax[b]] = tax_dict[current_tax[b]][-1] #get the lowest classification for each ASV
        new_col = pd.DataFrame(new_col, columns=[str(a), 'ASV']) #turn the classifications at the current level into a dataframe
        new_col = new_col.set_index(['ASV']) #set the ASV (note: this does not necessarily mean this is an ASV, rather it is the taxon names that we are using for this level) column as the index column
        df_tax = new_col.merge(right=df_tax, on='ASV') #and add it to the precious dataframe
        sort_levels.append(str(a)) #add the name of the column that we made with the classifications
    df_tax = df_tax.sort_values(sort_levels, ascending=False) #and now sort the dataframe based on these classifications (this will go in order, so it will keep all proteobacteria together, with all e.g. alpha-, gamma-proteobacteria)
    df_tax.drop(sort_levels, axis=1, inplace=True) #now remove these columns from the sorted dataframe
    return df_tax, names #return them along with the dataframe

def plot_colorbar(ax, cmap, df, name, fs=fs_main, orientation='horizontal', norm='rel_abun'):
    '''
    Function to add a colorbar to an axis
    Takes as input:
        - ax (the axis to plot the colorbar on - note that this will take up the whole axis, rather than just adding a padded version on the side)
        - cmap (the colormap to plot on the bar, as a string not a matplotlib colormap object)
        - df (the data that the colorbar corresponds to - this can be a list of values or a dataframe)
        - name (a string of the name that should me added as a label)
        - fs (the fontsize to use for the axis tick labels - default is the main size set at the top of the script)
        - orientation (string of whether to make this a horizontal or vertical colorbar, default is 'horizontal')
    '''
    #function to plot a colorbar on the axis ax using the colormap cmap
    #it also takes the dataframe of values that it should be based on, a name and the fontsize (using fs_main as default)
    mins, maxs = [], []
    if isinstance(df, list):
        mins.append(df[0])
        maxs.append(df[1])
    else:
        for a in df.columns: #get the minimum and maximum values for each column
            mins.append(df[a].min())
            maxs.append(df[a].max())
    mi = min(mins) #now get the minimum and maximum values overall
    ma = max(maxs)
    """
    if ma < 0:
        mi = 0
        ma = 0.1
    """
    cmap = mpl.cm.get_cmap(cmap, 256) #tell it which colormap to use
    if norm == 'clr':
      norm = mpl.colors.SymLogNorm(vmin=mi, vmax=ma, linthresh=0.03) #the minimum and maximum values to normalize to
    else:
      norm = mpl.colors.Normalize(vmin=mi, vmax=ma) #the minimum and maximum values to normalize to
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation=orientation) #get the colorbar for the heatmap
    plot_names = []
    for a in range(len(name)):
        if name[a] in rename_plots:
            plot_names.append(rename_plots[name[a]])
        else:
            plot_names.append(name[a])
    cb1.set_label(name, fontsize=fs) #set the label names
    cb1.ax.tick_params(labelsize=fs_small) #and the size of the ticks
    return
    
def get_heatmap_random_forest(df_rf, ft_in, colnames, basedir, level, other_levels, meta_name, predict, tax_dict, ASV_dict, other_folder=False, norm='rare'):
    '''
    This function will get an individual heatmap plot for one category of one random forest, by sorting out all of the data and ordering taxa etc
    It takes as input:
        - df_rf (dataframe containing list of ASVs in one column with a list of their importances in a second column)
        - ft_in (the original feature table containing a list of ASVs as rows and samples as columns)
        - colnames (the groupings within the metadata category)
        - basedir (the base directory to save everything to)
        - level (taxonomic level that we are plotting)
        - meta_name (the metadata category that we are plotting)
        - predict (the accuracy with which this category performed)  
    '''
    #ft = ft_in.transpose() #transpose the feature table
    #df_rf = df_rf.set_index(['ASV']) #in our list of
    ft = ft_in.reset_index() #reset the indices in our feature table
    ft = ft.rename(columns={'index':'ASV'}) #rename the index column
    ft = ft.set_index(['ASV']) #set the index column as ASV again
    numeric, non_numeric = [], []
    for a in colnames: #for each grouping name
        try: #if it is numeric
            colnames[a] = float(colnames[a]) #get it as a float
            colnames[a] = round(colnames[a]) #and then round it to the nearest whole number (the random forest regression was still based on the full values, but we do this to simplify plotting)
            numeric.append(colnames[a])
        except:
            colnames[a] = colnames[a] #otherwise, just leave it as it is (but the try: function doesn' work without an except: option)
            non_numeric.append(colnames[a])
    new_col_order = sorted(list(set(numeric))) + sorted(list(set(non_numeric))) #add together the numeric and non-numeric lists of column names, removing duplicates and sorting them in ascending order
    ft = ft.rename(columns=colnames) #rename the columns
    ft = ft.groupby(by=ft.columns, axis=1).mean() #if we have the same names now, group them together (taking a mean for each abundance value)
    ft = ft.reindex(new_col_order, axis=1) #re-index based on the order of our unique sorted columns
    df_rf.sort_values(by=['Importance'], ascending=False, inplace=True)
    df_rf = df_rf[df_rf.max(axis=1) > 0.00]
    if df_rf.shape[0] > 50:
      df_rf = df_rf[:50]
    df_rf_tax, n = sort_by_tax(df_rf, level, other_levels) #sort the taxa/ASVs in the random forest data frame by their higher phylogenetic grouping
    merge = df_rf_tax.merge(right=ft, on='ASV') #merge the importance dataframe with the feature table (this means we only keep those ASVs that are important)
    merge.drop('Importance', axis=1, inplace=True) #now remove the importance column again
    if norm not in ['log', 'clr']:
      merge = merge.div(merge.max(axis=1), axis=0) #normalise within each ASV (to make the colors visible for all even if there are large differences in abundance)
    #set up the figure size based on the number of metadata groupings we have, but checking that it won't be too small if we only have a couple of categories
    width = len(merge.columns)+3
    if width < 15:
        width = 15
    start_col2 = int(width-(width/5)*2)
    width_col = int(width/5)
    plt.figure(figsize=(max([width, 25])/3, 20))
    ax1 = plt.subplot2grid((merge.shape[0], width+5), (2,0), colspan=2, rowspan=merge.shape[0]-2)
    ax2 = plt.subplot2grid((merge.shape[0], width+5), (2,2), colspan=width, rowspan=merge.shape[0]-2)
    ax3 = plt.subplot2grid((merge.shape[0], width+5), (2,width+3), colspan=1, rowspan=int((merge.shape[0]-2)/6))
    ax4 = plt.subplot2grid((merge.shape[0], width+5), (int((merge.shape[0]-2)/3)+2, width+3), colspan=1, rowspan=int((merge.shape[0]-2)/3))
    #get the heatmaps for the taxon importances (df_rf_tax) and taxon abundances (merge) separately and with different color schemes
    if level == 7:
        asv_plot = {}
        prev_asv = list(df_rf_tax.index.values)
        for asv in prev_asv:
            if asv in ASV_dict:
                asv_plot[asv] = ASV_dict[asv]+': '+tax_dict[asv][-1]
            else:
                asv_plot[asv] = tax_dict[asv][-1]
        df_rf_tax_plot = df_rf_tax.rename(columns=rename_plots, index=asv_plot)
    else:
        df_rf_tax_plot = df_rf_tax.rename(columns=rename_plots, index=rename_plots)
    merge_plot = merge.rename(columns=rename_plots, index=rename_plots)
    annotate_heatmap(ax1, df_rf_tax_plot, cmap='inferno', rnd=3)
    if norm in ['clr', 'log']:
      annotate_heatmap(ax2, merge_plot, cmap='viridis', yticks=False, rnd=1, annotate=False, vmin=None, norm=norm)
    else:
      annotate_heatmap(ax2, merge_plot, cmap='viridis', yticks=False, rnd=1, annotate=False, norm=norm)
    #now get the appropriate color bars for each dataframe
    plot_colorbar(ax3, 'inferno', df_rf_tax_plot, orientation='vertical', name='Taxon importance')
    plot_colorbar(ax4, 'viridis', merge_plot, orientation='vertical', name=norm_names[norm]+' abundance', norm=norm)
    if meta_name in rename_plots:
        meta_name_plot = rename_plots[meta_name]
    else:
        meta_name_plot = meta_name
    plt.subplots_adjust(hspace=0.5)
    if other_folder == False:
        ax2.text(0.5, 1.025, meta_name_plot+'\n\n', va='center', ha='center', fontsize=fs_title, fontweight='bold', transform=ax2.transAxes)
        ax2.text(0.5, 1.025, '\nInformative '+phylo_level_names[level]+'\nClassification accuracy = '+str(predict)+'%', va='center', ha='center', fontsize=fs_main, fontweight='bold', transform=ax2.transAxes)
        plt.savefig(basedir+'/figures/random_forest/single_category/'+norm+'_'+phylo_level_names[level]+'_'+meta_name+'_single_forest'+ext, dpi=dpi, bbox_inches='tight')
    else: 
        of = other_folder.split('_')
        if of[0] in rename_plots:
            of[0] = rename_plots[of[0]]
        if len(of) > 1:
            if of[1] in rename_plots:
                of[1] = rename_plots[of[1]]
            of = of[0]+' '+of[1]
        else:
            of = of[0]
        ax2.set_title(of+'\n'+meta_name_plot+'\nInformative '+phylo_level_names[level]+'\nClassification accuracy = '+str(predict)+'%')
        plt.savefig(basedir+'/figures/random_forest/single_environment/'+norm+'_'+other_folder+'_'+phylo_level_names[level]+'_'+meta_name+'_single_forest'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return df_rf_tax

def get_random_forest_plots(ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=False, skip_individual=False, norm='rare'):
    '''
    Function to get all random forest plots, i.e. at all Phylogenetic levels and including the overall heatmap showing how accurate all random forests were for all metadata categories at all phylogenetic levels
    Takes as input:
        - ft (feature table as a dataframe with samples as columns and taxa as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to use to open random forest information as well as to save figures to)
        - other_folder (optional input string (file path) to tell the function to save to a folder that is not the default directory for random forests, default is False)
        - skip_individual (optional boolean, default is False. If True no individual heatmaps (i.e. for each metadata category and each phylogenetic level) will be made)
    '''
    asv = list(ft.index.values)
    count = 0
    for a in [1, 2, 3, 4, 5, 6, 7]:
        this_asv = list(asv)
        new_names, other_levels = {}, {}
        for b in range(len(this_asv)):
            this_tax = tax_dict[this_asv[b]]
            prev = ''
            for c in range(len(this_tax)):
                if this_tax[c] != '':
                    prev = this_tax[c]
                else:
                    this_tax[c] = prev
            if a < 7:
                new_names[this_asv[b]] = this_tax[a]
                this_asv[b] = this_tax[a]
            else:
                new_names[this_asv[b]] = this_asv[b]
            other_levels[this_asv[b]] = this_tax[:a]
        this_ft = ft.rename(index=new_names) #now rename the feature table using the new names that we have just defined
        this_ft = this_ft.groupby(by=this_ft.index, axis=0).sum() #group all of the names that are together, taking a sum for each phyla
        if other_folder == False:
            df_rf = pd.read_csv(basedir+'/random_forest/'+norm+'/'+phylo_level_names[a]+'_overall.csv', header=0, index_col=0)
        else:
            df_rf = pd.read_csv(basedir+'/random_forest/'+norm+'/single_environment/'+other_folder+'_'+phylo_level_names[a]+'.csv', header=0, index_col=0)
        meta_names = list(df_rf.columns)
        try:
          scores, oob_scores = df_rf.loc[['Score'], :], df_rf.loc[['OOB_score'], :]
          df_rf.drop(['Score', 'OOB_score'], axis=0, inplace=True)
        except:
          scores = df_rf.loc[['Score'], :]
          df_rf.drop(['Score'], axis=0, inplace=True)
        df_rf = df_rf.reset_index()
        df_rf = df_rf.rename(columns={'index':'ASV'})
        df_rf.index = df_rf['ASV']
        df_rf.drop(['ASV'], axis=1, inplace=True)
        df_rf_full = df_rf.copy()
        if a < 7:
            df_rf.rename(index={'Alteromonas macleodii':'Marinimicrobia (SAR406 clade)'}, inplace=True)
            df_rf_sorted, names = sort_by_tax(df_rf, a, other_levels, tax_dict=False, names=False)
        else:
            df_rf_sorted, names = sort_by_tax(df_rf, a, other_levels, tax_dict=tax_dict, names=True)
        df_rf_sorted = df_rf_sorted.groupby(by=df_rf_sorted.index.values, axis=0).mean()
        df_rf_sorted["Mean"] = df_rf_sorted.mean(axis=1)
        df_rf_sorted.sort_values(by=['Mean'], ascending=False, inplace=True)
        df_rf_sorted = df_rf_sorted[df_rf_sorted.max(axis=1) > 0.00]
        if df_rf_sorted.shape[0] > 50:
          df_rf_sorted = df_rf_sorted[:50]
        df_mean = df_rf_sorted.drop(meta_names, axis=1, inplace=False)
        df_rf_sorted.drop(['Mean'], axis=1, inplace=True)
        if other_folder == False:
            fig = plt.figure(figsize=(len(df_rf_sorted.columns)/3, 0.22*(df_rf_sorted.shape[0]+5)))
            
            ax1 = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (8,2), colspan=len(df_rf_sorted.columns), rowspan=2)
            ax2 = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (10,0), colspan=2, rowspan=df_rf_sorted.shape[0]-2)
            ax3 = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (10,2), colspan=len(df_rf_sorted.columns), rowspan=df_rf_sorted.shape[0]-2)
            ax4 = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (0,0), colspan=int(len(df_rf_sorted.columns)/5))
            
            n1 = int(len(df_rf_sorted.columns)-len(df_rf_sorted.columns)/5)+3
            ax5 = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (0,n1), colspan=int(len(df_rf_sorted.columns)/5))
            
            n2 = int(n1-len(df_rf_sorted.columns)/5)-1
            ax6 = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (0,n2), colspan=int(len(df_rf_sorted.columns)/5))
            
            ax_dendro = plt.subplot2grid((df_rf_sorted.shape[0]+5, len(df_rf_sorted.columns)+3), (4,2), rowspan=4, colspan=len(df_rf_sorted.columns), frameon=False)
            plt.sca(ax_dendro)
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            plt.tick_params(axis='y', which='both', left=False, right=False, labelbottom=False)
            
            Z = hierarchy.linkage(df_rf_sorted.transpose(), 'average', metric='braycurtis')
            mpl.rcParams['lines.linewidth'] = 2
            hierarchy.set_link_color_palette(['k'])
            dn = hierarchy.dendrogram(Z, above_threshold_color='k', orientation='top')
            
            y_labels, locs, xlocs, labels = list(ax_dendro.get_xticklabels()), list(ax_dendro.get_xticks()), list(ax2.get_yticks()), []
            for y in y_labels:
              labels.append(y.get_text())
            grouping = list(df_rf_sorted.columns)
            plot_labels = []
            for l in range(len(labels)):
              plot_labels.append(grouping[int(labels[l])-1])
            plt.xticks([])
            plt.yticks([])
            df_rf_sorted = df_rf_sorted.reindex(columns=plot_labels)
            scores = scores.reindex(columns=plot_labels)
            
        else:
            fig = plt.figure(figsize=(5, 0.22*(df_rf_sorted.shape[0]+5)))
            ax1 = plt.subplot2grid((df_rf_sorted.shape[0]+5, 9), (7,1), colspan=5, rowspan=1)
            ax2 = plt.subplot2grid((df_rf_sorted.shape[0]+5, 9), (9,0), colspan=1, rowspan=df_rf_sorted.shape[0]-3)
            ax3 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (9,1), colspan=5, rowspan=df_rf_sorted.shape[0]-3)
            ax4 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (0,6), colspan=3)
            ax5 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (3,6), colspan=3)
            ax6 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (6,6), colspan=3)
            
        #plot colorbars for each of the three sets of values (the sums, the importances within categories and the scores, respectively)
        scores = scores*100
        scores = scores.astype('int32')
        plot_colorbar(ax4, 'summer', df_mean, name='Mean taxon importance', fs=fs_small)
        plot_colorbar(ax5, 'inferno', df_rf_sorted, name='Taxon\nimportance', fs=fs_small)
        plot_colorbar(ax6, 'PuBu', scores, name='Classification\naccuracy (%)', fs=fs_small)
        if other_folder == False:
            #ax1.set_title('Top informative                              \n'+phylo_level_names[a]+'                              ', fontsize=fs_main)
            do_nothing = True
        else:
            of = other_folder.split('_')
            if of[0] in rename_plots:
                of[0] = rename_plots[of[0]]
            if len(of) > 1:
                if of[1] in rename_plots:
                    of[1] = rename_plots[of[1]]
                of = of[0]+' '+of[1]
            else:
                of = of[0]
            ax1.set_title('Top informative '+phylo_level_names[a]+'\n'+of, fontsize=fs_main) 
        #now get the heatmaps for each dataframe (same as for the colorbars)
        df_rf_sorted_names = df_rf_sorted.rename(columns=rename_plots, index=rename_plots)
        df_rf_sorted_names = df_rf_sorted_names.div(df_rf_sorted_names.max(axis=1), axis=0)
        
        if a == 7:
            mean_names = list(df_rf_sorted_names.index.values)
            new_names = {}
            for b in range(len(mean_names)):
                if mean_names[b] in tax_dict and mean_names[b] in ASV_dict:
                    new_names[mean_names[b]] = ASV_dict[mean_names[b]]+': '+tax_dict[mean_names[b]][-1]
            df_rf_sorted_names.rename(index=new_names, inplace=True)   
        annotate_heatmap(ax1, scores, cmap='PuBu', yticks=False, xticks=False, rnd=0, annotate=True)
        annotate_heatmap(ax2, df_mean, cmap='summer', yticks=False, xticks=True, rnd=3, annotate=True)
        annotate_heatmap(ax3, df_rf_sorted_names, cmap='inferno', yticks='right', xticks=True, rnd=3, annotate=False, italics=True)
        if other_folder == False:
            plt.savefig(basedir+'/figures/random_forest/'+phylo_level_names[a]+'_overall_forest_'+norm+ext, dpi=dpi, bbox_inches='tight') #save the figure
        #else:
        #    plt.savefig(basedir+'/figures/random_forest/single_environment/'+other_folder+'_'+phylo_level_names[a]+ext, dpi=dpi, bbox_inches='tight')
        plt.close()
        meta_order = open_txt(basedir+'metadata.txt')[0]
        del meta_order[0]
        for n in range(len(meta_names)):
            if skip_individual:
                continue
            mn = list(meta_names)
            mn.remove(meta_names[n])
            new_names = {}
            colnames = list(this_ft.columns)
            for b in range(len(colnames)):
                i = -10
                for c in range(len(meta_order)):
                  if meta_order[c] == meta_names[n]:
                    i = c
                new_names[colnames[b]] = meta_dict[colnames[b]][i]
            score = scores.loc['Score', meta_names[n]]
            this_df = df_rf_full.drop(mn, axis=1, inplace=False)
            try:
                this_df.drop(['Score', 'OOB_score'], axis=0, inplace=True)
            except:
                try:
                  this_df.drop(['Score'], axis=0, inplace=True)
                except:
                  do_nothing = True
            this_df = this_df.rename(columns={meta_names[n]:'Importance'})
            get_heatmap_random_forest(this_df, this_ft, new_names, basedir, a, other_levels, meta_names[n], score, tax_dict, ASV_dict, other_folder, norm=norm)
        scores = scores.rename(index={'Score':phylo_level_names[a]})
        if count == 0:
            all_scores = scores
            count += 1
        else:
            all_scores = pd.concat([scores, all_scores], sort=True)
    fig = plt.figure(figsize=(5, len(list(all_scores.columns))/5))
    ax1 = plt.subplot(111)
    all_scores = all_scores.transpose()
    all_scores = all_scores.rename(columns=rename_plots, index=rename_plots)
    all_scores = all_scores[all_scores.columns[::-1]]
    all_scores['Mean'] = all_scores.mean(axis=1)
    all_scores = all_scores.sort_values(by=['Mean'], ascending=True)
    all_scores.drop('Mean', axis=1, inplace=True)
    annotate_heatmap(ax1, all_scores, cmap='inferno', yticks=True, xticks=True, rnd=0, annotate=True)
    if other_folder == False:
        plt.savefig(basedir+'/figures/random_forest/overall_forest_scores_'+norm+ext, dpi=dpi, bbox_inches='tight')
    else:
        of = other_folder.split('_')
        if of[0] in rename_plots:
            of[0] = rename_plots[of[0]]
        if len(of) > 1:
            if of[1] in rename_plots:
                of[1] = rename_plots[of[1]]
            if of[1] != 'Laboratory' and of[0] != 'Freshwater':
                of = of[0]+' '+of[1]
            else:
                of = of[0]+'\n'+of[1]
        else:
            of = of[0]
        plt.ylabel(of, fontsize=fs_main, fontweight='bold')
        plt.savefig(basedir+'/figures/random_forest/single_environment/'+norm+'_'+other_folder+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return all_scores
    
def get_environment_random_forest_plots(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir, norm='rare', skip_individual=True):
    '''
    Function to make figures for all random forests calculated on PlasticTypeGeneral in different environments
    Takes as input:
        - ft (feature table as a dataframe with samples as columns and taxa as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to use to open random forest information as well as to save figures to)
    Uses the get_random_forest_plots function from above but using reduced feature and metadata tables 
    This could easily be edited to take lab and field separately as well as together
    '''
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    labfield, labindex = [['lab', 'field']], 5
    all_samples = list(ft.columns)
    all_samples_meta = list(meta_df.index.values)
    for a in range(len(envs)):
        for b in range(len(labfield)):
            if envs[a] == 'freshwater' and labfield[b] == 'lab':
                continue
            meta_names = list(meta_df.columns)
            keeping = []
            for c in range(len(all_samples)):
                if meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] == labfield[b]:
                    keeping.append(True)
                elif meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] in labfield[b]:
                    keeping.append(True)
                else:
                    keeping.append(False)
            single_ft = ft.loc[:, keeping]
            if len(single_ft.index.values) == 0:
                continue
            single_ft = single_ft[single_ft.max(axis=1) > 0]
            keeping_meta = []
            for d in range(len(all_samples_meta)):
                if meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] == labfield[b]:
                    keeping_meta.append(True)
                elif meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] in labfield[b]:
                    keeping_meta.append(True)
                else:
                    keeping_meta.append(False)
            single_meta_df = meta_df.loc[keeping_meta]
            single_meta_df.drop(meta_names, axis=1, inplace=True)
            if len(single_meta_df.index.values) == 0:
                continue
            name = envs[a]
            scores = get_random_forest_plots(single_ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=name, skip_individual=skip_individual, norm=norm)
            scores.rename(index={'Plastic type (general)':envs[a].capitalize()}, inplace=True)
            if a == 0:
                all_scores = scores
            else:
                all_scores = pd.concat([all_scores, scores], sort=False)
    return all_scores

def get_overall_random_forest_plot(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir, norm='rare'):
    '''
    Function to make figures for all random forests calculated on PlasticTypeGeneral in different environments and different metadata categories in all environments
    Takes as input:
        - ft (feature table as a dataframe with samples as columns and taxa as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to use to open random forest information as well as to save figures to)
    '''
    metadata_scores = get_random_forest_plots(ft, tax_dict, ASV_dict, meta_dict, basedir, skip_individual=True, norm=norm)
    env_scores = get_environment_random_forest_plots(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir, norm=norm)
    env_scores['Mean'] = env_scores.mean(axis=1)
    env_scores = env_scores.sort_values(by=['Mean'], ascending=True)
    env_scores.drop('Mean', axis=1, inplace=True)
    
    fig = plt.figure(figsize=(5, 8))
    ax1 = plt.subplot2grid((26, 1), (0, 0), rowspan=20)
    ax2 = plt.subplot2grid((26, 1), (22, 0), rowspan=4)
    
    annotate_heatmap(ax1, metadata_scores, cmap='inferno', yticks=True, xticks=False, rnd=0, annotate=True)
    annotate_heatmap(ax2, env_scores, cmap='inferno', yticks=True, xticks=True, rnd=0, annotate=True)
    ax1.set_title('A', loc='left', fontsize=fs_title, fontweight='bold')
    ax2.set_title('B', loc='left', fontsize=fs_title, fontweight='bold')
    ax1.set_title('Metadata categories', fontsize=fs_title, fontweight='bold')
    ax2.set_title('Plastic type (general)', fontsize=fs_title, fontweight='bold')
    plt.savefig(basedir+'/figures/Fig5_random_forest_'+norm+ext, dpi=dpi, bbox_inches='tight')
    return
    
def calculate_richness(ft_df_rare, ft_df_nr, ft_df_raw, meta_df, basedir):
    ft_df_raw.drop(['taxonomy'], axis=1, inplace=True)
    ft_df_raw = pd.DataFrame(ft_df_raw.sum(axis=0)).transpose()
    studies = sorted(list(set(meta_df.loc[:, 'Study'])))
    rich_rare = list(ft_df_rare[ft_df_rare > 0].count())
    rich_rare = pd.DataFrame([rich_rare], columns=ft_df_rare.columns, index=['Richness'])
    rich_nr = list(ft_df_nr[ft_df_nr > 0].count())
    rich_nr = pd.DataFrame([rich_nr], columns=ft_df_nr.columns, index=['Richness'])
    rename_samples = {}
    for sam in meta_df.index.values:
        rename_samples[sam] = meta_df.loc[sam, 'Study']
    rich_rare = rich_rare.rename(columns=rename_samples)
    rich_nr = rich_nr.rename(columns=rename_samples)
    ft_df_raw = ft_df_raw.rename(columns=rename_samples)
    
    plt.figure(figsize=(15,8))
    ax_reads = plt.subplot2grid((3,3), (0,0), colspan=2)
    ax1 = plt.subplot2grid((3,3), (1,0), colspan=2)
    ax2 = plt.subplot2grid((3,3), (2,0), colspan=2)
    ax_corr = plt.subplot2grid((3,3), (0,2), rowspan=3, colspan=1)
    
    ax_reads.set_title('A', fontweight='bold', loc='left')
    ax1.set_title('B', fontweight='bold', loc='left')
    ax2.set_title('C', fontweight='bold', loc='left')
    ax_corr.set_title('D', fontweight='bold', loc='left')
    
    count = 1
    counts = []
    plot_names = []
    for study in studies:
        if study not in rich_rare.columns:
            continue
        ax_reads.boxplot(ft_df_raw.loc[:, study], positions=[count])
        ax2.boxplot(rich_rare.loc[:, study], positions=[count])
        ax1.boxplot(rich_nr.loc[:, study], positions=[count])
        plot_names.append(rename_plots[study])
        counts.append(count)
        count += 1
        l1, l2 = list(ft_df_raw.loc[:, study].values), list(rich_nr.loc[:, study].values)
        for x in range(len(l1)):
            ax_corr.scatter(l1[x], l2[x], color='k', s=5)
    plt.sca(ax1)
    plt.xticks([])
    plt.sca(ax_reads)
    plt.xticks([])
    plt.semilogy()
    plt.sca(ax2)
    plt.xticks(counts, plot_names, rotation=90)
    plt.sca(ax_corr)
    plt.ylabel('Richness')
    plt.xlabel('Reads per sample')
    plt.semilogx()
    ax2.set_ylabel('Rarefied\nRichness'), ax1.set_ylabel('Not rarefied\nRichness'), ax_reads.set_ylabel('Not rarefied\nReads per sample')
    plt.tight_layout()
    plt.xlim(left=2000)
    plt.savefig(basedir+'/figures/FigS1_richness'+ext, bbox_inches='tight', dpi=600)
    return
    
def plot_individual_study(study_name, study_df, study_meta_df, study_w_uf, study_uw_uf, meta_dict, tax_dict, ASV_dict, basedir, est=10000, n_jobs=1):
    '''
    Get a figure for an individual study, including weighted and unweighted unifrac nMDS plots 
    and heatmaps to show the classification accuracy of random forest models as well as the top ASVs
    Input:
        - study_name (the name of the study that we are plotting, a string)
        - study_df (a feature table dataframe containing only samples for this study, with samples as columns and ASVs as rows)
        - study_meta_df (a dataframe containing all metadata for the samples in this study, with samples as rows and categories as columns. It should only contain the categories (columns) that you want to carry out random forest models for)
        - study_w_uf (dataframe containing the weighted unifrac distance matrix with only samples for this study)
        - study_uw_uf (dataframe containing the unweighted unifrac distance matrix with only samples for this study)
        - meta_dict (a dictionary with sample names as keys and all metadata information returned)
        - tax_dict (a dictionary containing ASVs as keys that returns all other taxonomic levels)
        - basedir (path to the base directory for analyses as a string)
        - est (number of estimators to use for random forests; default is 10,000)
        - n_jobs (number of threads to use for calculating random forests; default is 1)
    '''
    this_fol = basedir+'/separate_studies/'+study_name+'/' #define the folder name for this study
    if os.path.exists(this_fol): #if the folder already exists, do nothing
        do_nothing = True
    else: #otherwise, make the folder
        os.system("mkdir "+this_fol)
    #set up the figure and subplots, using 
    fig = plt.figure(figsize=(20,10))
    fig.suptitle(rename_plots[study_name], fontsize=fs_title+4, fontweight='bold')
    ax_nmds_w = plt.subplot2grid((2,3), (0,0))
    ax_nmds_uw = plt.subplot2grid((2,3), (1,0))
    ax_nmds_w.set_title('Weighted Unifrac', fontsize=fs_title, fontweight='bold')
    ax_nmds_uw.set_title('Unweighted Unifrac', fontsize=fs_title, fontweight='bold')
    uf, ax_uf = [study_w_uf, study_uw_uf], [ax_nmds_w, ax_nmds_uw]
    filter_on, filter_index, second_filter, second_filter_ind = 'none', 'none', '', ''
    source, source_index, color_src = ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic', 'blank'], 10, ['#5F9EA0', '#8B008B', '#3593FC', '#FCA94A', 'yellow', '#C3C3C3']
    scores, dfs = get_single_forest(study_df.transpose(), study_meta_df, this_fol+'/random_forest', 7, est=est, n_jobs=n_jobs)
    for a in range(len(uf)):
        names = list(uf[a].columns) #get a list of column names
        dist_matr = uf[a].rename_axis('ID').values
        save_name = this_fol+'nmds.df'+str(a)
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, source, source_index, ax_uf[a], 'best', color_src, names, meta_dict, save_name, second_filter, second_filter_ind, '', n_jobs=1)
        plt.sca(ax_uf[a])
        plt.legend(handles=handles, loc='best', fontsize=fs_main)
        plt.xlabel('nMDS1', fontsize=fs_small), plt.ylabel('nMDS2', fontsize=fs_small)
        plt.yticks(fontsize=fs_small), plt.xticks(fontsize=fs_small)
    if len(list(study_meta_df.columns)) == 0 or (isinstance(scores, list) and scores == []) or np.isnan(scores.loc['Score', :].max()) or scores.loc['Score', :].max() <= 0:
        cs = len(study_df.columns)
        if cs <= 8:
            cs = 16
        elif cs > 30:
            cs = 34
        ax1 = plt.subplot2grid((20,66), (0,30), rowspan=20, colspan=cs)
        ax_colbar1 = plt.subplot2grid((3,66), (0,30+cs+1))
        asv_tax = {}
        asv = list(study_df.index.values)
        for a in range(len(asv)):
            asv_tax[asv[a]] = ASV_dict[asv[a]]+': '+r'$'+tax_dict[asv[a]][-1]+'$'
        study_df.rename(index=asv_tax, inplace=True)
        study_df = study_df*100
        study_df["Mean"] = study_df.mean(axis=1)
        study_df.sort_values(by=['Mean'], ascending=False, inplace=True)
        study_df = study_df[:30]
        study_df = study_df.drop(['Mean'], axis=1, inplace=False)
        study_df = study_df.reindex(sorted(study_df.columns), axis=1)
        if cs > 20:
            rename = {}
            nums = []
            for a in range(len(study_df.columns)):
                rename[study_df.columns[a]] = a+1
                if a+1 % 10 == 0:
                    nums.append(a+1)
            study_df.rename(columns=rename, inplace=True)
            annotate_heatmap(ax1, study_df, cmap='viridis', yticks=True, xticks=False, rnd=3, annotate=False, italics=True)
            plt.sca(ax1)
            plt.xticks(nums)
        else:
            annotate_heatmap(ax1, study_df, cmap='viridis', yticks=True, xticks=True, rnd=3, annotate=True, italics=True)
        plot_colorbar(ax_colbar1, 'viridis', [0, max(study_df.max())], 'ASV relative abundance (%)', fs=fs_main, orientation='vertical')
        string = r'$\bf{Relative}$ $\bf{abundance}$ $\bf{of}$ $\bf{top}$ $\bf{30}$ $\bf{ASVs}$'
        ax1.set_title(string, fontsize=fs_title)
        plt.savefig(this_fol+'overall_plot'+ext, dpi=dpi, bbox_inches='tight')
        plt.subplots_adjust(hspace=0.2, wspace=0.2)
        plt.close()
        return
    ax_rf_scores = plt.subplot2grid((20,33), (0,13), colspan=19)
    ax_mean_importance = plt.subplot2grid((20,33), (2,13), rowspan=18)
    ax_rf_1_importance = plt.subplot2grid((20,33), (2,14), rowspan=18)
    ax_rf_1 = plt.subplot2grid((20,33), (2,15), rowspan=18, colspan=5)
    ax_rf_2_importance = plt.subplot2grid((20,33), (2,20), rowspan=18)
    ax_rf_2 = plt.subplot2grid((20,33), (2,21), rowspan=18, colspan=5)
    ax_rf_3_importance = plt.subplot2grid((20,33), (2,26), rowspan=18)
    ax_rf_3 = plt.subplot2grid((20,33), (2,27), rowspan=18, colspan=5)
    ax_colbar1 = plt.subplot2grid((3,66), (0,65))
    ax_colbar2 = plt.subplot2grid((3,66), (1,65))
    ax_colbar3 = plt.subplot2grid((3,66), (2,65))
    ax1, ax2, ax3 = [ax_rf_1_importance, ax_rf_1], [ax_rf_2_importance, ax_rf_2], [ax_rf_3_importance, ax_rf_3]
    axes = [ax1, ax2, ax3]
    try:
      scores.drop('OOB_score', axis=0, inplace=True)
      dfs.drop(['Score', 'OOB_score'], axis=0, inplace=True)
    except:
      dfs.drop(['Score'], axis=0, inplace=True)
    scores = scores*100
    scores_rename = scores.rename(columns=rename_plots)
    ax_rf_scores.set_title('Random forest classification accuracy (%)', fontsize=fs_title, fontweight='bold')
    annotate_heatmap(ax_rf_scores, scores_rename, cmap='PuBu', yticks=False, xticks='top', rnd=1, annotate=True, italics=False)
    meta_cats = list(scores.columns)
    scores_list = []
    for cat in meta_cats:
        scores_list.append(scores.loc['Score', cat])
    plot_colorbar(ax_colbar1, 'PuBu', [0, max(scores_list)], 'Classification accuracy (%)', fs=fs_main, orientation='vertical')
    scores_list = [x for _,x in sorted(zip(scores_list,meta_cats))]
    scores_list.reverse()
    scores_list, drop = scores_list[:3], scores_list[3:]
    count, count1 = 0, 0
    for a in range(len(scores_list)):
        if 'incubation' in scores_list[a] or 'Incubation' in scores_list[a]:
            if count > 0:
                rem = str(scores_list[a])
                scores_list[a] = drop[0]
                drop.remove(drop[0])
                drop.append(rem)
            else:
                count += 1
        if 'type' in scores_list[a] or 'Type' in scores_list[a]:
            if count1 > 0:
                rem = str(scores_list[a])
                scores_list[a] = drop[0]
                drop.remove(drop[0])
                drop.append(rem)
            else:
                count1 += 1
    dfs.drop(drop, axis=1, inplace=True)
    dfs["Mean"] = dfs.mean(axis=1)
    dfs.sort_values(by=['Mean'], ascending=False, inplace=True)
    dfs = dfs[:30]
    df_mean = dfs.drop(scores_list, axis=1, inplace=False)
    dfs.drop(['Mean'], axis=1, inplace=True)
    asv_tax = {}
    asv = list(dfs.index.values)
    for a in range(len(asv)):
        asv_tax[asv[a]] = ASV_dict[asv[a]]+': '+r'$'+tax_dict[asv[a]][-1]+'$'
    df_mean = df_mean.reset_index()
    df_mean.rename(columns={'index':'ASV'}, inplace=True)
    df_mean = df_mean.set_index('ASV')
    study_df = study_df.merge(right=df_mean, on='ASV')
    study_df.sort_values(by=['Mean'], ascending=False, inplace=True)
    study_df.drop('Mean', axis=1, inplace=True)
    df_mean.rename(columns=rename_plots, index=asv_tax, inplace=True)
    annotate_heatmap(ax_mean_importance, df_mean, cmap='inferno', yticks=True, xticks=True, rnd=3, annotate=True, italics=True)
    plot_colorbar(ax_colbar2, 'inferno', [0, df_mean.iloc[0,0]], 'ASV importance', fs=fs_main, orientation='vertical')
    mas = []
    this_dfs = []
    for a in range(len(scores_list)):
        ax = axes[a]
        this_df = dfs.copy()
        mn = list(this_df.columns)
        mn.remove(scores_list[a])
        this_df.drop(mn, axis=1, inplace=True)
        this_df.rename(columns={scores_list[a]:'Importance'}, inplace=True)
        annotate_heatmap(ax[0], this_df, cmap='inferno', yticks=False, xticks=True, rnd=3, annotate=True)
        rename_cols = {}
        samples = list(study_df.columns)
        for s in range(len(samples)):
            group = study_meta_df.loc[samples[s], scores_list[a]]
            rename_cols[samples[s]] = group
        this_study_df = study_df.rename(columns=rename_cols, inplace=False)
        this_study_df.rename(columns=rename_plots, inplace=True)
        this_study_df = this_study_df.groupby(by=this_study_df.columns, axis=1).mean()
        this_study_df = this_study_df*100
        mas.append(max(this_study_df.max()))
        if scores_list[a] in rename_plots:
            scores_list[a] = rename_plots[scores_list[a]]
        if a == 1:
            string = r'$\bf{Top}$ $\bf{3}$ $\bf{random}$ $\bf{forest}$ $\bf{models}$ $\bf{and}$ $\bf{important}$ $\bf{ASVs}$'
            ax[1].set_title(string+'\n'+scores_list[a], fontsize=fs_title)
        else:
            ax[1].set_title(scores_list[a], fontsize=fs_title)
        this_dfs.append(this_study_df)
    annotate_heatmap(ax_rf_1, this_dfs[0], cmap='viridis', yticks=False, xticks=True, rnd=1, annotate=True, vmax=max(mas))
    if len(this_dfs) > 1:
        annotate_heatmap(ax_rf_2, this_dfs[1], cmap='viridis', yticks=False, xticks=True, rnd=1, annotate=True, vmax=max(mas))
        if len(this_dfs) > 2:
            annotate_heatmap(ax_rf_3, this_dfs[2], cmap='viridis', yticks=False, xticks=True, rnd=1, annotate=True, vmax=max(mas))
        else:
            ax_rf_3_importance = plt.subplot2grid((20,33), (2,26), rowspan=18, frameon=False)
            ax_rf_3 = plt.subplot2grid((20,33), (2,27), rowspan=18, colspan=5, frameon=False)
            plt.sca(ax_rf_3_importance)
            plt.yticks([])
            plt.xticks([])
            plt.sca(ax_rf_3)
            plt.yticks([])
            plt.xticks([])
    else:
        ax_rf_2_importance = plt.subplot2grid((20,33), (2,20), rowspan=18, frameon=False)
        ax_rf_2 = plt.subplot2grid((20,33), (2,21), rowspan=18, colspan=5, frameon=False)
        plt.sca(ax_rf_2_importance)
        plt.yticks([])
        plt.xticks([])
        plt.sca(ax_rf_2)
        plt.yticks([])
        plt.xticks([])
        ax_rf_3_importance = plt.subplot2grid((20,33), (2,26), rowspan=18, frameon=False)
        ax_rf_3 = plt.subplot2grid((20,33), (2,27), rowspan=18, colspan=5, frameon=False)
        plt.sca(ax_rf_3_importance)
        plt.yticks([])
        plt.xticks([])
        plt.sca(ax_rf_3)
        plt.yticks([])
        plt.xticks([])
    plot_colorbar(ax_colbar3, 'viridis', [0, max(mas)], 'ASV relative abundance (%)', fs=fs_main, orientation='vertical')
    plt.subplots_adjust(hspace=0.4, wspace=0.5)
    plt.savefig(this_fol+'overall_plot'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def plots_per_study(ft_df, meta_df, meta_dict, w_uf, uw_uf, tax_dict, ASV_dict, basedir, est=10000, n_jobs=1):
    '''
    Get figures for each study, showing weighted and unweighted unifrac distances between samples, random forest models based on ASVs and
    plots of the top 3 random forest models and 30 most informative ASVs
    Input:
        - ft_df (a dataframe of a feature table with all samples across all studies, with ASVs as rows and samples as columns)
        - meta_df (a dataframe with samples as rows and metadata categories as columns)
        - meta_dict (a dictionary with sample names as keys and returning all sample metadata)
        - w_uf (the file name of a weighted unifrac distance matrix for all samples)
        - uw_uf (the file name of an unweighted unifrac distance matrix for all samples)
        - tax_dict (a dictionary containing ASVs as keys that returns all other taxonomic levels)
        - basedir (path to the base directory for analyses as a string)
        - est (number of estimators to use for random forests; default is 10,000)
        - n_jobs (number of threads to use for calculating random forests; default is 1)
    '''
    studies = sorted(list(set(meta_df['Study']))) #get the unique study names from the metadata file
    samples = list(ft_df.columns) #get a list of all sample names
    samples_meta = list(meta_df.index.values) #get a list of all samples names in the metadata file
    #w_uf = pd.read_csv(w_uf, header=0, index_col=0) #read in the weighted unifrac distance matrix
    #uw_uf = pd.read_csv(uw_uf, header=0, index_col=0) #read in the unweighted unifrac distance matrix
    cols_uf = list(w_uf.columns) #get a list of sample names from the distance matrix
    ft_df = ft_df.reset_index() #rest the index column
    ft_df.rename(columns={'index':'ASV'}, inplace=True) #rename the previous index column
    ft_df = ft_df.set_index('ASV') #now reset the original index column, now named 'ASV'
    adding = False
    for a in range(len(studies)): #for each study
        if os.path.exists(basedir+'/separate_studies/'+studies[a]+'/overall_plot'+ext):
          print('Already had the plot for this study: '+studies[a]+"\nDelete the plot (and the 'random_forest.csv' file) if you'd like to rerun this\n")
          continue
        keeping, keeping_meta, keeping_uf = [], [], [] #set up lists
        #check each of the feature table, metadata dataframe and unifrac sample names to see whether they are in this study. 
        #if they are, add True to the appropriate list, if not, add False
        for b in range(len(samples)): 
            if meta_dict[samples[b]][0] == studies[a]:
                keeping.append(True)
            else:
                keeping.append(False)
        for c in range(len(samples_meta)):
            if meta_dict[samples_meta[c]][0] == studies[a]:
                keeping_meta.append(True)
            else:
                keeping_meta.append(False)
        for d in range(len(cols_uf)):
            if meta_dict[cols_uf[d]][0] == studies[a]:
                keeping_uf.append(True)
            else:
                keeping_uf.append(False)
        study_df = ft_df.loc[:, keeping] #now get only those that we are keeping as a new dataframe for the feature table
        study_df = study_df[study_df.max(axis=1) > 0] #and only keep ASVs that are present in this study
        study_w_uf = w_uf.loc[keeping_uf, keeping_uf] #now get only those we are keeping as a new dataframe for the weighted unifrac
        study_uw_uf = uw_uf.loc[keeping_uf, keeping_uf] #now get only those we are keeping as a new dataframe for the unweighted unifrac
        study_meta_df = meta_df.loc[keeping_meta, :] #now get only those we are keeping as a new dataframe for the metadata dataframe
        meta_cols = list(study_meta_df.columns) #get a list of the metadata categories
        not_adding = []
        for col in meta_cols: #for each of those columns
            vals = list(study_meta_df[col]) #get a list of unique values
            if len(list(set(vals))) < 2: #if there are fewer than two unique values in the category, then don't use this for random forests
                not_adding.append(col) #so add it to not adding
        study_meta_df.drop(not_adding, axis=1, inplace=True) #and drop those columns that we don't want from the dataframe
        plot_individual_study(studies[a], study_df, study_meta_df, study_w_uf, study_uw_uf, meta_dict, tax_dict, ASV_dict, basedir, est=est, n_jobs=n_jobs) #now get the plot for this study
    return
    
def group_ft_level(ft, level, tax_dict, basedir, rename=False, saving=False): #This script was written specifically for grouping before ancom, but the rename option is there to make it compatible with other functions that might want to regroup but not rename to an ASV (this renaming is only to have a representative sequence for each group)
    '''
    Function to group a feature table to a particular phylogenetic level
    Takes as input:
        - ft (dataframe with sample names as columns and ASVs as rows)
        - level (integer indicating the phylogenetic level to group to, where 1 is phylum and 7 is ASV)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - basedir (name of the directory to save the figures to)
        - rename (boolean indicating whether to rename the grouped ASVs or not, default is False)
        - saving (boolean indicating whether to save the resulting taxonomy file (this is used to rename to the level that we have used for plots made in R), default is False)
    '''
    asv = list(ft.index.values) #get a list of ASVs
    asv_tax_dict, asv_tax = {}, [] #set up dictionary and list
    for a in range(len(asv)): #for each ASV
        asv_tax_dict[asv[a]] = tax_dict[asv[a]][level] #get the taxonomy at the level that we want, and add this to a dictionary
        asv_tax.append(tax_dict[asv[a]][level]) #also add it to a list
    unique_tax = list(set(asv_tax)) #get a list of the unique taxa at this level
    representative_asv, rep_asv_dict = [], {} #set up a dictionary and list
    for a in range(len(unique_tax)): #for each of the unique taxa
        for b in range(len(asv_tax)): #for each ASV
            if unique_tax[a] == asv_tax[b]: #if tthe ASV is the unique taxa that we are looking at
                representative_asv.append([asv[b], unique_tax[a]]) #add this ASV (i.e. the first) as the representative ASV for this taxa
                rep_asv_dict[unique_tax[a]] = asv[b] #and add it to the dictionary
                break #now break this loop, as we only need to do this once
    ft = ft.rename(index=asv_tax_dict) #now rename the feature table at the level that we want
    ft = ft.groupby(by=ft.index, axis=0).sum() #and group all ASVs that now have the same name together (taking a sum of the abundance for each)
    if rename: #if we're renaming the taxa
        ft = ft.rename(index=rep_asv_dict) #rename the feature table again, with the representative ASV for each taxa
        if saving: #if we are saving this file
            with open(basedir+'taxonomy_name_only.csv', 'w') as f: #open the csv file and add the rows to it
                writer = csv.writer(f)
                writer.writerow(['OTUID', 'Species name'])
                for asv in representative_asv:
                    writer.writerow(asv)
    return ft

def tree_heatmap(ft, meta_dict, basedir, tax_dict, level=7):
    '''
    Function to perform ANCOM tests for statistical significance at a given phylogenetic level between treatments at different time points, and plot the results of this as a tree and heatmap using an R script
    Takes as input:
        - ft (dataframe containing sample names as columns and ASVs as rows)
        - meta_dict (dictionary containing metadata for all samples, with sample names as keys)
        - basedir (name of the directory to save the figures to)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - level (phylogenetic level at which to perform the comparison, where 1 is phylum and 7 is ASV)
    Returns:
        - plots for each environment at the given phylogenetic level. Figure names will take the form environment_ancom_significance_level.pdf
    '''
    lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'specific_species']
    ft = ft*100 #multiply the relative abundance for each value in the feature table by 100 
    if level < 7: #if this isn't at the ASV level, then group by taxa
        ft = group_ft_level(ft, level, tax_dict, basedir+'/ancom/', rename=True, saving=True)
    samples = list(ft.columns) #get a list of sample names
    env, source, inc_time, source_inc, source_inc_dict = [], [], [], [], {} #set up lists for environment, incubation time, source and a combination of these
    envs = ['marine', 'freshwater', 'aquatic', 'terrestrial']
    #os.rename(basedir+'agglom/reduced_tree.tree', 'reduced_tree.tree') #rename the tree to save it to the directory that contains the script (this means that we don't need to put anything additional into the R script, rather than needing to change the file names etc.)
    for a in range(len(samples)): #for each sample
        env.append(meta_dict[samples[a]][3]) #get the environment that it came from, using the metadata dictionary
        source.append(meta_dict[samples[a]][10]) #get the sample type/source for this sample
        inc_time.append(meta_dict[samples[a]][13]) #get the incubation time
        source_inc.append(meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13]) #add a combination of the source and incubation time to another list
        source_inc_dict[samples[a]] = meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13] #and also to a dictionary
    #set up lists of the comparisons to be carried out, as well as the names to use
    comparisons = [['aliphatic early', 'aliphatic late'], ['other plastic early', 'other plastic late'], ['biofilm early', 'biofilm late'],
                   ['aliphatic early', 'other plastic early'], ['aliphatic early', 'unknown plastic early'], ['aliphatic early', 'biofilm early'],
                   ['aliphatic late', 'other plastic late'], ['aliphatic late', 'unknown plastic late'], ['aliphatic late', 'biofilm late'],
                   ['other plastic early', 'unknown plastic early'], ['other plastic early', 'biofilm early'],
                   ['other plastic late', 'unknown plastic late'], ['other plastic late', 'biofilm late'],
                   ['unknown plastic early', 'biofilm early'], ['unknown plastic late', 'biofilm late']]
    comp_names = ['Aliphatic: early vs late', 'Other plastic: early vs late', 'Biofilm: early vs late',
                  'Early: aliphatic vs other plastic', 'Early: aliphatic vs unknown plastic', 'Early: aliphatic vs biofilm', 
                  'Late: aliphatic vs other plastic', 'Late: aliphatic vs unknown plastic', 'Late: aliphatic vs biofilm', 
                  'Early: other plastic vs unknown plastic', 'Early: other plastic vs biofilm', 
                  'Late: other plastic vs unknown plastic', 'Late: other plastic vs biofilm',
                  'Early: unknown plastic vs biofilm', 'Late: unknown plastic vs biofilm']
    for a in range(len(envs)): #for each environment
        if os.path.exists(basedir+'/figures/ancom/'+envs[a]+'_ancom_significance_'+lvls[level]+'.pdf'):
          print('Already had this one: '+envs[a]+' '+lvls[level])
          continue
        this_env, keeping = [], []
        for b in range(len(env)): #for each environment in the list of environments that the samples belong to
            if env[b] == envs[a]: #if it is the same as the environment that we are looking at
                keeping.append(True) #add True to the keeping list
            else:
                keeping.append(False) #otherwise, add False to the keeping list
        env_ft = ft.loc[:, keeping] #get a new dataframe with only those samples that are in envs[a]
        env_ft.rename(columns=source_inc_dict, inplace=True) #rename the remaining samples using the dictionary of sample type and incubation time
        all_comps, colors_list = [], [[], []] #set up lists
        for c in range(len(comparisons)): #for each comparison in the list of comparisons
            try: #try means that the script won't break if this doesn't work - this means that if we don't have the right samples for the comparison in this environment, the script will just continue to the next comparison
                comp_env_ft = env_ft.loc[:, comparisons[c]] #only keep the samples with names matching this comparison
                if len(comp_env_ft.loc[:, comparisons[c][0]].columns) < 2 or len(comp_env_ft.loc[:, comparisons[c][1]].columns) < 2: #if we don't have enough samples remaining
                    all_comps.append(False) #add False to all_comps
                    continue #and continue to the next comparison
                else:
                    all_comps.append(True) #otherwise, add True to all_comps
            except: #if this broke the script
                all_comps.append(False) #add False to all_comps
                continue #and continue to the next comparison
            comp_env_ft = comp_env_ft.transpose() #transpose the feature table
            comp_env_ft = comp_env_ft[comp_env_ft.columns[comp_env_ft.max() > 0]] #only keep taxa if they are present in at least one sample
            comp_env_ft = comp_env_ft.replace(to_replace=0,value = 0.0001) #if we have any zeroes, replace these with 
            comp_env_ft = comp_env_ft.fillna(value=0.0001) #if we have any na values, fill these with 0.0001
            #perform the ANCOM comparison with this transposed feature table, a list of the sample names and holm-bonferroni false discovery rate correction:
            ancom_df, percentile_df = ancom(comp_env_ft, pd.Series(list(comp_env_ft.index.values), index=list(comp_env_ft.index.values)), multiple_comparisons_correction='holm-bonferroni')
            ancom_results, percentile_results, tax_names = ancom_df.values.tolist(), percentile_df.values.tolist(), ancom_df.index.values #get lists of the parts of the ancom data that we are interested in
            significant, medians, not_sig, not_sig_medians = [], [], [], [] #set up some lists
            for d in range(len(ancom_results)): #for each of the ancom results (i.e. each taxa)
                if ancom_results[d][1] == True: #if the ancom result is significant
                    if percentile_results[d][2] == 0.0001 and percentile_results[d][7] == 0.0001: #if the percentile results aren't both 0.0001 (i.e. added in artificially)
                        not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]]) #add the medians to the list of non-significant results
                        not_sig.append(tax_names[d]) #add the name of this taxa to the list of non-significant results
                    else:
                        significant.append(tax_names[d]) #otherwise, add the name of this taxa to the list of significant results
                        medians.append([percentile_results[d][2], percentile_results[d][7]]) #add the medians to the list of significant results
                else: #if the ancom result is not significant
                    not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]]) #add the medians to the list of non-significant results
                    not_sig.append(tax_names[d]) #add the name of this taxa to the list of non-significant results
            differences = []
            for e in range(len(medians)): #for each of the significant sets of medians
                if 0 in medians[e] or 1 in medians[e]: #if either value is 0 or 1
                    medians[e][0] += 0.0001 #add 0.0001 to both numbers
                    medians[e][1] += 0.0001
                if medians[e][0] == medians[e][1]: #if the medians are equal to one another
                    diff = 0 #make the difference be 0
                else: #otherwise
                    diff = math.log2(medians[e][1])/math.log2(medians[e][0]) #calculate the log2 fold change
                    diff = math.pow(2, diff) #and then convert the number back to not log2
                    if diff < 1 and diff > 0: #if the difference is between 0 and 1, then make this a negative number and divide by 1
                        diff = -(1/diff)
                differences.append(diff) #add the difference to the list
            if len(significant) < 1: #if we only have one significant result
                all_comps[c] = False #don't bother plotting this competition
                continue #and continue
            if isinstance(this_env, list): #if this_env is a list (i.e. we haven't made  a dataframe with the results yet)
                this_env = pd.DataFrame(data=differences, index=significant, columns=[comp_names[c]]) #get a dataframe with the results of this competition
            else: #otherwise
                new_df = pd.DataFrame(data=differences, index=significant, columns=[comp_names[c]]) #make a new dataframe with the results of this competition
                this_env = pd.concat([this_env, new_df], axis=1, sort=False) #and add it to the dataframe with the results of the other competitions
            str1, str2 = comparisons[c][0].split(" "), comparisons[c][1].split(" ") #split the comparisons based on spaces, so we can look at the first and second parts of each
            if str1[-1] == 'early' and str2[-1] == 'late': #if the comparison was with the same sample type, between early and late samples
                col1, col2 = color_source['early'], color_source['late'] #set the color to be that of early/late in the dictionary
            else: #otherwise
                if len(str1) > 2: #if the length of the comparison is more than two (the sample name was longer/too long after splitting based on spaces)
                    str1[0] = str1[0]+' '+str1[1] #then add the names back together
                    str1[1] = str1[2]
                    str1 = str1[:-1]
                if len(str2) > 2: #otherwise if this is true for the second part, do as above
                    str2[0] = str2[0]+' '+str2[1]
                    str2[1] = str2[2]
                    str2 = str2[:-1]
                col1 = color_source[str1[0]] #get the color for the first sample type
                col2 = color_source[str2[0]] #and the second sample type 
            colors_list[0].append(str(col1)+", "+str(col2)), colors_list[1].append(comp_names[c]) #add a string to the color list
        if sum(all_comps) > 0: #if we had results for some competitions
            cols_df = pd.DataFrame(data=[colors_list[0]], index=['Colors'], columns=colors_list[1]) #get the colors as a dataframe
            this_env = this_env.append(cols_df) #add this dataframe to the one of all comparisons for this environment
            this_env = this_env.fillna(value=0) #and fill in any na values
            this_env.to_csv(basedir+'/ancom/'+envs[a]+'_'+lvls[level]+'_ancom_significant.csv') #save this dataframe with the differences for each taxa, for each comparison
            this_env.to_csv(basedir+'ancom_significant.csv') #and save it as a generic file, that the R script can get to
            ft_fn = basedir+'/ancom/'+envs[a]+'_'+lvls[level]+'_ancom_significant.csv'
            tree_fn = basedir+'agglom/reduced_tree.tree'
            tax_fn = basedir+"taxonomy_name_only.csv"
            save_name = basedir+'/figures/ancom/'+envs[a]+'_ancom_significance_'+lvls[level]+'.pdf'
            r.suppressWarnings(r.ancom_tree())
            os.rename(basedir+'ancom_heatmap.pdf', basedir+'/figures/ancom/'+envs[a]+'_ancom_significance_'+lvls[level]+'.pdf') #rename the output of the R script
    try:
      os.remove('ancom_significant.csv') #remove this from the folder with the scripts (each environment will be saved in the ancom folder anyway)
    except:
      do_nothing = True
    #os.remove('taxonomy_name_only.csv') #remove this, too
    #os.rename('reduced_tree.tree', basedir+'agglom/reduced_tree.tree') #and put the reduced tree back where it came from
    return
```
Note: basedir should be the directory containing the files listed above and you should also change n_jobs to be the number of processors that you want to use (this will affect the speed with which many functions run) and est to the number of estimators that you want to use for the random forest sections. 1 will run very quickly but will not be robust - the analyses in the paper used 10,000.

R:
```{R, R_input_functions, results='hide', fig.keep='all', message=FALSE, warning=FALSE}
#"agglom/reduced_tree.tree"
r_unifrac <- function(tree_file) {
  #Import all of the data (these are large files so this may take a while)
  asv_table <- read.csv(paste(py$basedir, "grouped_samples_for_unifrac.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, tree_file, sep=""))
  
  #Convert these to phyloseq objects
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  w_uf <- UniFrac(physeq, weighted=TRUE, normalized=FALSE, fast=TRUE)
  
  w_uf_df <- as.data.frame(as.matrix(w_uf))
  write.csv(w_uf_df, paste(py$basedir, "/weighted_unifrac_grouped_samples.csv", sep=""))
}

metacoder <- function(colors_metacoder, norm) {

color1 = colors_metacoder[1]
color2 = colors_metacoder[2]

asvs <- read.csv(paste(py$basedir, "metacoder.csv", sep=""), header=TRUE, sep = ",")

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

if (norm == 'clr') {
obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund", cols = cn, groups = treats, combinations=list(c('Treat1', 'Treat2')), func=function(abund_1, abund_2) {
  log_ratio <- median(abund_1) - median(abund_2)
  if (is.nan(log_ratio)) {
    log_ratio <- 0
  }
  list(log2_median_ratio = log_ratio,
       median_diff = median(abund_1) - median(abund_2),
       mean_diff = mean(abund_1) - mean(abund_2),
       wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
})
} else if (norm == 'log') {
  obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund", cols = cn, groups = treats, combinations=list(c('Treat1', 'Treat2')), func=function(abund_1, abund_2) {
  log_ratio <- median(abund_1) - median(abund_2)
  if (is.nan(log_ratio)) {
    log_ratio <- 0
  }
  list(log2_median_ratio = log_ratio,
       median_diff = median(abund_1) - median(abund_2),
       mean_diff = mean(abund_1) - mean(abund_2),
       wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
})
} else {
  obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund", cols = cn, groups = treats, combinations=list(c('Treat1', 'Treat2')))
}
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method="fdr")
range(obj$data$diff_table$wilcox_p_value, finite=TRUE)
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0 
heat_tree(obj, node_label=taxon_names, node_size=n_obs, node_color=log2_median_ratio, node_color_interval=c(-3,3), node_color_range=c(color2, middle, color1), layout="davidson-harel", initial_layout="reingold-tilford", make_node_legend=FALSE, make_edge_legend=FALSE, output_file=paste(py$basedir, "metacoder.pdf", sep=""), node_label_max=0) 
heat_tree(obj, node_label=taxon_names, node_size=n_obs, node_color=log2_median_ratio, node_color_interval=c(-3,3), node_color_range=c(color2, middle, color1), layout="davidson-harel", initial_layout="reingold-tilford", make_node_legend=FALSE, make_edge_legend=FALSE, output_file=paste(py$basedir, "metacoder_labels.pdf", sep=""))
dev.off()
}

ancom_tree <- function() {
  library("microbiome")
  library("phyloseq")
  ft_fn = paste(py$basedir, 'ancom_significant.csv', sep='')
  tree_fn = paste(py$basedir, 'agglom/reduced_tree_rare.tree', sep='')
  tax_fn = paste(py$basedir, 'ancom/taxonomy_name_only.csv', sep='')
  save_name = paste(py$basedir, 'ancom_heatmap.pdf', sep='')
  asv_table <- read.csv(ft_fn, sep=",", row.names=1, check.names=FALSE)
  asv_table = as.matrix(asv_table)
  colors = asv_table['Colors', ]
  asv_table <- asv_table[!rownames(asv_table) %in% c('Colors'), ]
  class(asv_table) <- "numeric"
  
  phy_tree <- read_tree(tree_fn)
  taxonomy <- read.csv(tax_fn, sep=",", row.names=1)
  taxonomy = as.matrix(taxonomy)
  asv_df = as.data.frame(asv_table)
  tax_df = as.data.frame(taxonomy, stringsAsFactors=F)
  
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  tree = phy_tree(physeq)
  
  detach("package:microbiome", unload=TRUE, force=TRUE)
  detach("package:phyloseq", unload=TRUE, force=TRUE)
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
  p
  ggsave(save_name)
}
```

#### (II) Make empty folders, if they don't already exist:
```{python, make_folders, results='hide', fig.keep='all', message=FALSE}
folder_names = ["agglom", "picrust", "figures", "ancom", "separate_studies", "figures/ancom", "figures/metacoder", "random_forest", "random_forest/rare", "random_forest/log", "random_forest/clr", "random_forest/rel_abun", "random_forest/rare/single_environment", "random_forest/log/single_environment", "random_forest/clr/single_environment", "random_forest/rel_abun/single_environment", "figures/random_forest", "figures/random_forest/single_environment", "figures/random_forest/single_category"]
for fn in folder_names:
    if not os.path.exists(basedir+"/"+fn):
       os.system("mkdir "+basedir+"/"+fn)
```

#### (*III) Reformat QIIME2 output files for R:
```{python, reformat_output_Q2, results='hide', fig.keep='all', message=FALSE}
if os.path.exists(basedir+'tax_dict_rare.dictionary'):
    tax_dict_rare = open_pickle(basedir+'tax_dict_rare.dictionary')
else:
    ft_rare, tax_dict_rare = format_R(ft_tax_rare, basedir, 'rare')

if os.path.exists(basedir+'tax_dict_nr.dictionary'):
    tax_dict_nr = open_pickle(basedir+'tax_dict_nr.dictionary')
else:
    ft_nr, tax_dict_nr = format_R(ft_tax_nr_filt, basedir, 'nr')
```

#### (IV) Perform agglomeration and calculate distances:

Do agglomeration:
```{R, agglom, results='hide', fig.keep='all', cache=TRUE, message=FALSE}
#first the rarefied data
if (!file.exists(paste(py$basedir, "agglom/otu_table_rare.csv", sep=""))) {
  if (!file.exists(paste(py$basedir, "agglom/reduced_tree_rare.tree", sep=""))) {
    asv_table <- read.csv(paste(py$basedir, "feature_table_rare.csv", sep=""), sep=",", row.names=1)
    asv_table = as.matrix(asv_table)
    phy_tree <- read_tree(paste(py$basedir, "qiime_output/tree_rare.nwk", sep=""))
    
    #Convert these to phyloseq objects
    ASV = otu_table(asv_table, taxa_are_rows = TRUE)
    physeq = phyloseq(ASV,phy_tree)
    
    #Remove ASVs below 1% relative abundance (without doing this, we have too many taxa to agglomerate and everything will most likely just crash)
    rel_abun <- transform_sample_counts(physeq, function(x) x/sum(x))
    remove_low_abun <- filter_taxa(rel_abun, function(x) max(x) > 0.01, TRUE)
    
    #See how many ASVs we had before and after removing the low abundance
    physeq
    remove_low_abun
    
    #Now do the agglomeration - 0.2 will have the least (~2000 from 7401 for the original data), while 0.05 will keep the most tips (~5000). You can find more information on doing this here: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/tip_glom
    #0.2 is the default value, but for the analyses shown in the paper we have compromised on 0.1, which gives 4469 taxa remaining (tips)
    agglom <- tip_glom(remove_low_abun, h=0.1)
    
    #And look at how many taxa we now have:
    agglom
    
    tree = phy_tree(agglom)
    write_phyloseq(agglom, type="all", path=paste(py$basedir, "agglom/", sep=""))
    file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_rare.csv", sep=""))
    ape::write.tree(tree, paste(py$basedir, "agglom/reduced_tree_rare.tree", sep=""))
  }
}

if (!file.exists(paste(py$basedir, "agglom/otu_table_nr.csv", sep=""))) {
  if (!file.exists(paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))) {
    #Add the rest in here at some point!
      asv_table <- read.csv(paste(py$basedir, "feature_table_nr.csv", sep=""), sep=",", row.names=1)
      asv_table = as.matrix(asv_table)
      phy_tree <- read_tree(paste(py$basedir, "qiime_output/tree_nr.nwk", sep=""))
          
      #Convert these to phyloseq objects
      ASV = otu_table(asv_table, taxa_are_rows = TRUE)
      
      #See how many ASVs we have
      ASV
      #24,873 for this one for me
      
      #For this there are too many taxa in the original object, so we split this up into three to make it manageable (not sure if this is a general R issue or an issue with running out of memory) Keeping the numbers to under ~15,000 seems to fix this, so you may need to adjust this depending on how many you have originally
      ASV_1 = ASV[1:8000]
      ASV_2 = ASV[8000:16000]
      ASV_3 = ASV[16000:24873]
      physeq_1 = phyloseq(ASV_1,phy_tree)
      physeq_2 = phyloseq(ASV_2,phy_tree)
      physeq_3 = phyloseq(ASV_3,phy_tree)
      
      #Agglomerate each separately
      agglom_1 <- tip_glom(physeq_1, h=0.1)
      agglom_2 <- tip_glom(physeq_2, h=0.1)
      agglom_3 <- tip_glom(physeq_3, h=0.1)
      
      #these are now each approximately 5,000, so we can add all of them together
      agglom_1
      agglom_2
      agglom_3
      
      #first save them
      write_phyloseq(agglom_1, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_nr_1.csv", sep=""))
      write_phyloseq(agglom_2, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_nr_2.csv", sep=""))
      write_phyloseq(agglom_3, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_nr_3.csv", sep=""))
      
      #and then make them into otu tables again
      asv_1 <- as.matrix(otu_table(agglom_1))
      asv_2 <- as.matrix(otu_table(agglom_2))
      asv_3 <- as.matrix(otu_table(agglom_3))
      
      #add them together and convert them to a phyloseq object
      asv <- rbind(asv_1, asv_2, asv_3)
      asv <- otu_table(asv, taxa_are_rows = TRUE)
      physeq <- phyloseq(asv, phy_tree)
      
      #and agglomerate
      agglom <- tip_glom(physeq, h=0.1)
      
      #Check how many ASVs we have now (12635)
      agglom
      
      #save this one
      write_phyloseq(agglom, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_nr.csv", sep=""))
      
      #and write the tree 
      tree = phy_tree(agglom)
      ape::write.tree(tree, paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))
      
      #Do the transformations:
      #relative abundance transformation
      physeq_ra <- transform_sample_counts(agglom, function(x) x/sum(x))
      write_phyloseq(physeq_ra, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_relabun.csv", sep=""))
              
      #log transformation (with pseudocount)
      physeq_log <- transform_sample_counts(agglom, function (x) log10(x+1))
      write_phyloseq(physeq_log, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_log.csv", sep=""))
              
      #CLR transformation (with pseudocount)
      physeq_clr <- microbiome::transform(agglom, "clr")
      write_phyloseq(physeq_clr, type="all", path=paste(py$basedir, "agglom/", sep=""))
      file.rename(paste(py$basedir, "agglom/otu_table.csv", sep=""), paste(py$basedir, "agglom/otu_table_clr.csv", sep=""))
  }
}
```

Make a taxonomy matrix for R:
```{python, tax_matrix, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
ft_nr = open_and_sort(basedir+'agglom/otu_table_nr.csv')
tax_matrix = pd.read_csv(basedir+'tax_dict_nr.csv', header=0, index_col=0)
tax_matrix = tax_matrix.loc[ft_nr.index.values, :]
tax_matrix.drop(['Species name'], axis=1, inplace=True)
tax_matrix = tax_matrix.reset_index()
```

Calculate distances:
```{R, unifrac_PHILR, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
#Get the unifrac distances for the rarefied data
if (!file.exists(paste(py$basedir, "agglom/weighted_unifrac_rare.csv", sep=""))) {
  #read in files again
  asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_rare.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_rare.tree", sep=""))
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  
  #calculate distances
  w_uf <- UniFrac(agglom, weighted=TRUE, normalized=FALSE, fast=TRUE)
  w_uf_df <- as.data.frame(as.matrix(w_uf))
  #And finally we can write the data to file so that we can continue our analysis in Python
  write.csv(w_uf_df, paste(py$basedir, "agglom/weighted_unifrac_rare.csv", sep=""))
}

if (!file.exists(paste(py$basedir, "agglom/unweighted_unifrac_rare.csv", sep=""))) {
  #read in files again
  asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_rare.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_rare.tree", sep=""))
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  
  #calculate distances
  uw_uf <- UniFrac(agglom, weighted=FALSE, normalized=FALSE, fast=TRUE)
  uw_uf_df <- as.data.frame(as.matrix(uw_uf))
  #And finally we can write the data to file so that we can continue our analysis in Python
  write.csv(uw_uf_df, paste(py$basedir, "agglom/unweighted_unifrac_rare.csv", sep=""))
}

#Get the unifrac distances for the relative abundance data
if (!file.exists(paste(py$basedir, "agglom/weighted_unifrac_rel_abun.csv", sep=""))) {
  #read in files again
  asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_relabun.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  
  #calculate distances
  w_uf <- UniFrac(physeq, weighted=TRUE, normalized=FALSE, fast=TRUE)
  w_uf_df <- as.data.frame(as.matrix(w_uf))
  #And finally we can write the data to file so that we can continue our analysis in Python
  write.csv(w_uf_df, paste(py$basedir, "agglom/weighted_unifrac_rel_abun.csv", sep=""))
}

if (!file.exists(paste(py$basedir, "agglom/unweighted_unifrac_rel_abun.csv", sep=""))) {
  #read in files again
  asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_relabun.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  
  #calculate distances
  uw_uf <- UniFrac(physeq, weighted=FALSE, normalized=FALSE, fast=TRUE)
  uw_uf_df <- as.data.frame(as.matrix(uw_uf))
  #And finally we can write the data to file so that we can continue our analysis in Python
  write.csv(uw_uf_df, paste(py$basedir, "agglom/unweighted_unifrac_rel_abun.csv", sep=""))
}

#Get the unifrac distances for the log data
if (!file.exists(paste(py$basedir, "agglom/weighted_unifrac_log.csv", sep=""))) {
  #read in files again
  asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_log.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  
  #calculate distances
  w_uf <- UniFrac(physeq, weighted=TRUE, normalized=FALSE, fast=TRUE)
  w_uf_df <- as.data.frame(as.matrix(w_uf))
  #And finally we can write the data to file so that we can continue our analysis in Python
  write.csv(w_uf_df, paste(py$basedir, "agglom/weighted_unifrac_log.csv", sep=""))
}

if (!file.exists(paste(py$basedir, "agglom/unweighted_unifrac_log.csv", sep=""))) {
  #read in files again
  asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_log.csv", sep=""), sep=",", row.names=1)
  asv_table = as.matrix(asv_table)
  phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))
  ASV = otu_table(asv_table, taxa_are_rows = TRUE)
  physeq = phyloseq(ASV,phy_tree)
  
  #calculate distances
  uw_uf <- UniFrac(physeq, weighted=FALSE, normalized=FALSE, fast=TRUE)
  uw_uf_df <- as.data.frame(as.matrix(uw_uf))
  #And finally we can write the data to file so that we can continue our analysis in Python
  write.csv(uw_uf_df, paste(py$basedir, "agglom/unweighted_unifrac_log.csv", sep=""))
}

#Get the PHILR distances for the not rarefied data
if (!file.exists(paste(py$basedir, "agglom/philr_distance.csv", sep=""))) {
    asv_table <- read.csv(paste(py$basedir, "agglom/otu_table_clr.csv", sep=""), sep=",", row.names=1)
    asv_table = as.matrix(asv_table)
    taxmat <- py$tax_matrix
    taxmat2 <- taxmat[,-1]
    rownames(taxmat2) <- taxmat[,1]
    sampledata <- py$meta_df
    sampledata <- sample_data(sampledata)
    colnames(taxmat2) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    phy_tree <- read_tree(paste(py$basedir, "agglom/reduced_tree_nr.tree", sep=""))
    
    #Convert these to phyloseq objects
    ASV = otu_table(asv_table, taxa_are_rows = TRUE)
    TAX = tax_table(taxmat2)
    taxa_names(TAX) <- taxmat[,1]
    physeq = phyloseq(ASV,phy_tree,TAX, sampledata)
    is.rooted(phy_tree(physeq))
    is.binary.tree(phy_tree(physeq))
    phy_tree(physeq) <- makeNodeLabel(phy_tree(physeq), method="number", prefix='n')
    name.balance(phy_tree(physeq), tax_table(physeq), 'n1')
    
    #Add pseudocount 
    physeq <- transform_sample_counts(physeq, function(x) x+1)
    
    #now the philr part
    otu.table <- t(otu_table(physeq))
    tree <- phy_tree(physeq)
    metadata <- sample_data(physeq)
    tax <- tax_table(physeq)
    physeq.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
    
    #now calculate the distance
    physeq.dist <- dist(physeq.philr, method="euclidean")
    physeq.dist.mat <- as.matrix(physeq.dist)
    write.csv(physeq.dist.mat, paste(py$basedir, "agglom/philr_distance.csv", sep=""))
}
```

#### (*V) Read these files into Python again and open the information that we need:

```{python, read_agglom_uf_files, results='hide', fig.keep='all', message=FALSE, warning=FALSE}
ft_df_rare = open_and_sort(basedir+'/agglom/otu_table_rare.csv')
ft_df_log = open_and_sort(basedir+'/agglom/otu_table_log.csv')
ft_df_clr = open_and_sort(basedir+'/agglom/otu_table_clr.csv')
ft_df_rel_abun = open_and_sort(basedir+'/agglom/otu_table_relabun.csv')
ft_df_nr = open_and_sort(basedir+'/agglom/otu_table_nr.csv')
tree_agglom_nr = basedir+'/agglom/reduced_tree_nr.tree'
tree_agglom_rare = basedir+'/agglom/reduced_tree_rare.tree'
tax_dict_rare = open_pickle(basedir+'tax_dict_rare.dictionary')
tax_dict_nr = open_pickle(basedir+'tax_dict_nr.dictionary')
tax_matrix_nr = pd.read_csv(basedir+'tax_dict_nr.csv')
tax_matrix_rare = pd.read_csv(basedir+'tax_dict_rare.csv')

ft_full_nr = basedir+'/feature_table_nr.csv'
ft_full_rare = basedir+'/feature_table_rare.csv'
meta, meta_names, meta_dict = get_meta(meta_fn)
meta_df = get_meta_df(meta, meta_names, list(ft_df_rare.columns))

#Read in distance matrices
w_uf_rare, uw_uf_rare = open_and_sort(basedir+'/agglom/weighted_unifrac_rare.csv'), open_and_sort(basedir+ '/agglom/unweighted_unifrac_rare.csv') #file names for unifrac distances based on agglomerated, rarefied data
w_uf_rel_abun, uw_uf_rel_abun = open_and_sort(basedir+'/agglom/weighted_unifrac_rel_abun.csv'), open_and_sort(basedir+ '/agglom/unweighted_unifrac_rel_abun.csv') #file names for unifrac distances based on agglomerated, relative abundance data
w_uf_log, uw_uf_log = open_and_sort(basedir+'/agglom/weighted_unifrac_log.csv'), open_and_sort(basedir+ '/agglom/unweighted_unifrac_log.csv') #file names for unifrac distances based on agglomerated, log transformed data
philr_clr = open_and_sort(basedir+'/agglom/philr_distance.csv')
eucl_dist = get_eucl_dist(ft_df_clr.transpose())
eucl_dist_df = pd.DataFrame(eucl_dist, index=ft_df_clr.columns, columns=ft_df_clr.columns)
eucl_dist_df.to_csv(basedir+'agglom/euclidean_distance.csv')
eucl_clr = open_and_sort(basedir+'/agglom/euclidean_distance.csv')

w_uf_rare.to_csv(basedir+'agglom/sorted_weighted_unifrac_rare.csv'), uw_uf_rare.to_csv(basedir+'agglom/sorted_unweighted_unifrac_rare.csv')
w_uf_log.to_csv(basedir+'agglom/sorted_weighted_unifrac_log.csv'), uw_uf_log.to_csv(basedir+'agglom/sorted_unweighted_unifrac_log.csv')
w_uf_rel_abun.to_csv(basedir+'agglom/sorted_weighted_unifrac_rel_abun.csv'), uw_uf_rel_abun.to_csv(basedir+'agglom/sorted_unweighted_unifrac_rel_abun.csv')
philr_clr.to_csv(basedir+'agglom/sorted_philr_distance.csv')
eucl_clr.to_csv(basedir+'agglom/sorted_euclidean_distance.csv')

if os.path.exists(basedir+'ASV_dict_nr.dictionary'):
    ASV_dict_nr = open_pickle(basedir+'ASV_dict_nr.dictionary')
else:
    ASV_dict_nr = get_ASV_dict(ft_df_nr, seqs, basedir)
    write_pickle(basedir+'ASV_dict_nr.dictionary', ASV_dict_nr)
if os.path.exists(basedir+'ASV_dict_rare.dictionary'):
    ASV_dict_rare = open_pickle(basedir+'ASV_dict_rare.dictionary')
else:
    ASV_dict_rare = get_ASV_dict(ft_df_rare, seqs_rare, basedir)
    write_pickle(basedir+'ASV_dict_rare.dictionary', ASV_dict_rare)
```
Note here that R seems to sometimes save things with a different encoding than is expected by pandas. If you get errors related to NaNs here, just try opening the .csv files in excel and resaving them as .csv files (not actually changing anything, the default for excel is to save in utf-8, which is what pandas wants) should fix this. It otherwise reads in the dataframes as all NaNs sometimes. 

#### Now we are ready to do all of the statistical analysis and make all of the plots.

#### (VI) Get basic study map and metrics (Figure 1):
```{python, map_metrics, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
if not os.path.exists(basedir+'/figures/Fig1_map_metrics'+ext):
  study_map_metrics(study_dates, study_locs, basedir, map_img, meta_df) 
```
This should save Figure 1 as 'figures/Fig1_map_metrics.png' (if you have changed the extension then this will be different).

And also the richness/number of reads for all studies (Figure S1):
```{python, richness, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
if os.path.exists(basedir+'feature_table_raw.csv'):
    ft_df_raw = open_and_sort(basedir+'feature_table_raw.csv')
else:
    ft_df_raw, tax_dict_raw = format_R(ft_tax_nr_filt, basedir, 'raw')

if not os.path.exists(basedir+'figures/FigS1_richness'+ext):
  calculate_richness(ft_df_rare, ft_df_nr, ft_df_raw, meta_df, basedir)
```

This should save Figure S1 as 'figures/FigS1_richness.png' (if you have changed the extension then this will be different).

#### (VII) Get the NMDS plots (Figure 2):
First get the groupings:
```{python, stats_groups, results='hide', fig.keep='all', message=FALSE}
def get_stats_groups(uf, meta_dict):
  group_env, group_study = [], []
  for a in uf.index.values:
    meta = meta_dict[a]
    group_env.append(meta[3]), group_study.append(meta[0])
  return group_env, group_study

if os.path.exists(basedir+'anosim_permanova_all.list'):
  stats = open_pickle(basedir+'anosim_permanova_all.list')
  run_tests = False
else:
  group_env, group_study = get_stats_groups(w_uf_rare, meta_dict)
  run_tests = True
```

Carry out the stats tests:
```{R, perm_anosim, results='hide', fig.keep='all', message=FALSE}
if (!file.exists(paste(py$basedir, "anosim_permanova_all.list", sep=""))) {
  #rarefied data
  as_uw_study <- anosim(py$uw_uf_rare, py$group_study, permutations=999)
  as_uw_env <- anosim(py$uw_uf_rare, py$group_env, permutations=999)
  as_w_study <- anosim(py$w_uf_rare, py$group_study, permutations=999)
  as_w_env <- anosim(py$w_uf_rare, py$group_env, permutations=999)

  perm_uw_study <- adonis(t(py$uw_uf_rare) ~ py$group_study, data=py$uw_uf_rare)
  perm_uw_env <- adonis(t(py$uw_uf_rare) ~ py$group_env, data=py$uw_uf_rare)
  perm_w_study <- adonis(t(py$w_uf_rare) ~ py$group_study, data=py$w_uf_rare)
  perm_w_env <- adonis(t(py$w_uf_rare) ~ py$group_env, data=py$w_uf_rare)

  rare_stats <- c(c(as_w_env$statistic, as_w_env$signif, perm_w_env$aov.tab$F.Model[1], perm_w_env$aov.tab$`Pr(>F)`[1]), c(as_w_study$statistic, as_w_study$signif, perm_w_study$aov.tab$F.Model[1], perm_w_study$aov.tab$`Pr(>F)`[1]), c(as_uw_env$statistic, as_uw_env$signif, perm_uw_env$aov.tab$F.Model[1], perm_uw_env$aov.tab$`Pr(>F)`[1]), c(as_uw_study$statistic, as_uw_study$signif, perm_uw_study$aov.tab$F.Model[1], perm_uw_study$aov.tab$`Pr(>F)`[1]))
  
  #relative abundance data
  as_uw_study <- anosim(py$uw_uf_rel_abun, py$group_study, permutations=999)
  as_uw_env <- anosim(py$uw_uf_rel_abun, py$group_env, permutations=999)
  as_w_study <- anosim(py$w_uf_rel_abun, py$group_study, permutations=999)
  as_w_env <- anosim(py$w_uf_rel_abun, py$group_env, permutations=999)
  
  perm_uw_study <- adonis(t(py$uw_uf_rel_abun) ~ py$group_study, data=py$uw_uf_rel_abun)
  perm_uw_env <- adonis(t(py$uw_uf_rel_abun) ~ py$group_env, data=py$uw_uf_rel_abun)
  perm_w_study <- adonis(t(py$w_uf_rel_abun) ~ py$group_study, data=py$w_uf_rel_abun)
  perm_w_env <- adonis(t(py$w_uf_rel_abun) ~ py$group_env, data=py$w_uf_rel_abun)
  
  rel_abun_stats <- c(c(as_w_env$statistic, as_w_env$signif, perm_w_env$aov.tab$F.Model[1], perm_w_env$aov.tab$`Pr(>F)`[1]), c(as_w_study$statistic, as_w_study$signif, perm_w_study$aov.tab$F.Model[1], perm_w_study$aov.tab$`Pr(>F)`[1]), c(as_uw_env$statistic, as_uw_env$signif, perm_uw_env$aov.tab$F.Model[1], perm_uw_env$aov.tab$`Pr(>F)`[1]), c(as_uw_study$statistic, as_uw_study$signif, perm_uw_study$aov.tab$F.Model[1], perm_uw_study$aov.tab$`Pr(>F)`[1]))
  
  #Log-transformed data
  as_uw_study <- anosim(py$uw_uf_log, py$group_study, permutations=999)
  as_uw_env <- anosim(py$uw_uf_log, py$group_env, permutations=999)
  as_w_study <- anosim(py$w_uf_log, py$group_study, permutations=999)
  as_w_env <- anosim(py$w_uf_log, py$group_env, permutations=999)
  
  perm_uw_study <- adonis(t(py$uw_uf_log) ~ py$group_study, data=py$uw_uf_log)
  perm_uw_env <- adonis(t(py$uw_uf_log) ~ py$group_env, data=py$uw_uf_log)
  perm_w_study <- adonis(t(py$w_uf_log) ~ py$group_study, data=py$w_uf_log)
  perm_w_env <- adonis(t(py$w_uf_log) ~ py$group_env, data=py$w_uf_log)
  
  log_stats <- c(c(as_w_env$statistic, as_w_env$signif, perm_w_env$aov.tab$F.Model[1], perm_w_env$aov.tab$`Pr(>F)`[1]), c(as_w_study$statistic, as_w_study$signif, perm_w_study$aov.tab$F.Model[1], perm_w_study$aov.tab$`Pr(>F)`[1]), c(as_uw_env$statistic, as_uw_env$signif, perm_uw_env$aov.tab$F.Model[1], perm_uw_env$aov.tab$`Pr(>F)`[1]), c(as_uw_study$statistic, as_uw_study$signif, perm_uw_study$aov.tab$F.Model[1], perm_uw_study$aov.tab$`Pr(>F)`[1]))

  #CLR-transformed data
  as_eucl_study <- anosim(py$eucl_clr, py$group_study, permutations=999)
  as_eucl_env <- anosim(py$eucl_clr, py$group_env, permutations=999)
  as_philr_study <- anosim(py$philr_clr, py$group_study, permutations=999)
  as_philr_env <- anosim(py$philr_clr, py$group_env, permutations=999)
  
  perm_eucl_study <- adonis(t(py$eucl_clr) ~ py$group_study, data=py$eucl_clr)
  perm_eucl_env <- adonis(t(py$eucl_clr) ~ py$group_env, data=py$eucl_clr)
  perm_philr_study <- adonis(t(py$philr_clr) ~ py$group_study, data=py$philr_clr)
  perm_philr_env <- adonis(t(py$philr_clr) ~ py$group_env, data=py$philr_clr)
  
  clr_stats <- c(c(as_eucl_env$statistic, as_eucl_env$signif, perm_eucl_env$aov.tab$F.Model[1], perm_eucl_env$aov.tab$`Pr(>F)`[1]), c(as_eucl_study$statistic, as_eucl_study$signif, perm_eucl_study$aov.tab$F.Model[1], perm_eucl_study$aov.tab$`Pr(>F)`[1]), c(as_philr_env$statistic, as_philr_env$signif, perm_philr_env$aov.tab$F.Model[1], perm_philr_env$aov.tab$`Pr(>F)`[1]), c(as_philr_study$statistic, as_philr_study$signif, perm_philr_study$aov.tab$F.Model[1], perm_philr_study$aov.tab$`Pr(>F)`[1]))
}
```

And then make the NMDS plots:
```{python, uf_nmds, results='hide', fig.keep='all', message=FALSE, cache=TRUE, warning=FALSE}
if run_tests:
  stats_rare = r.rare_stats
  stats_rel_abun = r.rel_abun_stats
  stats_log = r.log_stats
  stats_clr = r.clr_stats
  stats = [stats_rare, stats_rel_abun, stats_log, stats_clr]
  write_pickle(basedir+'anosim_permanova_all.list', stats)

count = 0
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  if not os.path.exists(basedir+'/figures/Fig2_nmds_overall_'+norm+ext):
    this_stats = [[stats[count][:4], stats[count][4:8]], [stats[count][8:12], stats[count][12:]]]
    if norm != 'clr':
      nmds_plot_study_env('agglom/sorted_weighted_unifrac_'+norm+'.csv', 'agglom/sorted_unweighted_unifrac_'+norm+'.csv', meta_dict, basedir, stats=this_stats, norm_name=norm)
    else:
      nmds_plot_study_env('agglom/sorted_euclidean_distance.csv', 'agglom/sorted_philr_distance.csv', meta_dict, basedir, stats=this_stats, norm_name=norm)
  count += 1
```

This should save Figure 2 as 'figures/Fig2_nmds_overall_rare.png' (if you have changed the extension then this will be different).

It will also save the NMDS plots that use the other normalisation methods, but these are only shown in [Supplementary Section 2](https://doi.org/10.6084/m9.figshare.12915317).

#### (VIII) Get the heatmaps using weighted and unweighted unifrac distances that shows the average distance within and between studies (Figure 3):
```{python, uf_heatmaps, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
dist_mats = [[w_uf_rare, uw_uf_rare], [w_uf_rel_abun, uw_uf_rel_abun], [w_uf_log, uw_uf_log], [eucl_clr, philr_clr]]
count = 0
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  if not os.path.exists(basedir+'/figures/Fig3_heatmap_combined_'+norm+ext):
    similarity_heatmap_combined(dist_mats[count][0], dist_mats[count][1], basedir, norm_name=norm) 
  count += 1
```

This should save Figure 3 as 'figures/Fig3_unifrac_heatmap_combined_rare.png' (if you have changed the extension then this will be different).

It will also plot the heatmaps for the other normalisation methods, but these are only shown in [Supplementary Section 2](https://doi.org/10.6084/m9.figshare.12915317). 

#### (IX) Get the figure that summarises sample types, groupings and the number of taxa shared within environments and sample types (Figure 4):
```{python, dendro_heatmap, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
count = 0
ft_full_nr = basedir+'/feature_table_nr.csv'
ft_full_rare = basedir+'/feature_table_rare.csv'
full_fts = [ft_full_rare, ft_full_nr, ft_full_nr, ft_full_nr]
norm_fts = [ft_df_rare, ft_df_rel_abun, ft_df_log, ft_df_clr]
tax_dicts = [tax_dict_rare, tax_dict_nr, tax_dict_nr, tax_dict_nr]
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  if not os.path.exists(basedir+'/figures/Fig4_dendro_venn_'+norm+ext):
    bar_dendro_venn(norm_fts[count], full_fts[count], meta_dict, basedir, tax_dicts[count], str_norm=norm)
  count += 1
```

This should save Figure 4 as 'figures/Fig4_dendro_venn_rare.png' (if you have changed the extension then this will be different).

It will also make these plots for the other normalisation methods, but these are only shown in [Supplementary Section 2](https://doi.org/10.6084/m9.figshare.12915317).

#### (X) Get the overall random forests:
This part should skip for you if you already have the files (and will tell you that it didn't do it because it already had a file with that name saved) - if you want to re-run it then delete the contents of the random_forests folder.
*Note that if you re-run this section then you might get slightly different results than the paper due to different random selections of 20%/80% of the data for testing/training*
```{python, get_rf_models_overall, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
get_random_forests(ft_df_rare, tax_dict_rare, meta_df, basedir, est=est, n_jobs=n_jobs, norm='rare')
get_random_forests(ft_df_log, tax_dict_nr, meta_df, basedir, est=est, n_jobs=n_jobs, norm='log')
get_random_forests(ft_df_clr, tax_dict_nr, meta_df, basedir, est=est, n_jobs=n_jobs, norm='clr')
get_random_forests(ft_df_rel_abun, tax_dict_nr, meta_df, basedir, est=est, n_jobs=n_jobs, norm='rel_abun')
```

#### (XI) Get the random forests for each environment for general plastic type:
```{python, get_rf_models_env, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
#add norm name to these
get_environment_random_forest(ft_df_rare, tax_dict_rare, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='rare') 
get_environment_random_forest(ft_df_log, tax_dict_nr, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='log') 
get_environment_random_forest(ft_df_clr, tax_dict_nr, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='clr') 
get_environment_random_forest(ft_df_rel_abun, tax_dict_nr, meta_df, meta_dict, basedir, est=est, n_jobs=n_jobs, norm='rel_abun') 
```

#### (XII) Get the random forest plots comparing the different normalisation methods for the overall random forests:
```{python, get_rf_plots_comparison, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
if not os.path.exists(basedir+'figures/FigS2_RF_compare'+ext):
  get_rf_comparison(basedir)
```

This should save Figure S2 as 'figures/FigS2_RF_comp.png' (if you have changed the extension then this will be different).

#### (XIII) Get the random forest plots comparing the different normalisation methods for the random forests for each environment on plastic type:
```{python, get_rf_plots_comparison_env, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
if not os.path.exists(basedir+'/figures/FigS3_RF_env_compare'+ext):
  get_rf_comparison_env(basedir)
```

This should save Figure S3 as 'figures/FigS3_RF_env_compare.png' (if you have changed the extension then this will be different).

#### (XIV) Get the random forest figures: 
```{python, get_rf_plots_overall, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
fts = [ft_df_rare, ft_df_rel_abun, ft_df_log, ft_df_clr]
tax_dicts = [tax_dict_rare, tax_dict_nr, tax_dict_nr, tax_dict_nr]
ASV_dicts = [ASV_dict_rare, ASV_dict_nr, ASV_dict_nr, ASV_dict_nr]
count = 0
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  get_overall_random_forest_plot(fts[count], meta_df, tax_dicts[count], ASV_dicts[count], meta_dict, basedir, norm=norm)
  count += 1
```

This should save Figure 5 as 'figures/Fig5_random_forest_rare.png' (if you have changed the extension then this will be different).

It will also save the same figure for the other normalisation methods, but these are only shown in [Supplementary Section 2](https://doi.org/10.6084/m9.figshare.12915317).

This bit doesn't have any checks to see if it has been run already, and will take 15 minutes or so to run, to change eval to =TRUE if you want to run it:
Run  the bottom lines to get the overall plots for each taxonomic level (for each normalisation method) and the top to get all of the plots for each metadata category.
```{python, get_rf_plots_separate, cache=TRUE}
get_random_forest_plots(ft_df_rare, tax_dict_rare, ASV_dict_rare, meta_dict, basedir, norm='rare')
get_random_forest_plots(ft_df_clr, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='clr')
get_random_forest_plots(ft_df_log, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='log')
get_random_forest_plots(ft_df_rel_abun, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='rel_abun')

get_random_forest_plots(ft_df_rare, tax_dict_rare, ASV_dict_rare, meta_dict, basedir, skip_individual=True, norm='rare')
get_random_forest_plots(ft_df_clr, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, skip_individual=True, norm='clr')
get_random_forest_plots(ft_df_log, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, skip_individual=True, norm='log')
get_random_forest_plots(ft_df_rel_abun, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, skip_individual=True, norm='rel_abun')

os.system('cp '+basedir+'/figures/random_forest/ASVs_overall_forest_rare.png '+basedir+'/figures/FigS4_rare_ASV_accuracy.png')
```

This is a lot of figures, so I haven't plotted them all here, but they are in Supplementary Sections [2](https://doi.org/10.6084/m9.figshare.12915317) (log, relative abundance and CLR-normalised) and [3](https://doi.org/10.6084/m9.figshare.12233759) (rarefied data).
The figure that we show in the supplementary (Figure S4) should also get copied/renamed to the main figures folder.

#### (XV) Get the environment random forest figures for general plastic type (including carrying out ANCOM tests for statistical significance): 
```{python, get_rf_plots_env, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
get_environment_random_forest_plots(ft_df_rare, meta_df, tax_dict_rare, ASV_dict_rare, meta_dict, basedir, norm='rare', skip_individual=False)
get_environment_random_forest_plots(ft_df_log, meta_df, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='log', skip_individual=False)
get_environment_random_forest_plots(ft_df_rel_abun, meta_df, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='rel_abun', skip_individual=False)
get_environment_random_forest_plots(ft_df_clr, meta_df, tax_dict_nr, ASV_dict_nr, meta_dict, basedir, norm='clr', skip_individual=False)
if not os.path.exists(basedir+'figures/FigS5_environment_random_forest_0.005_rare.pdf'):
  if not os.path.exists('figures/FigS5_environment_random_forest_0.005_rare'+ext):
    make_env_rf_plot(ft_df_rare, tax_dict_rare, basedir, ASV_dict_rare, meta_dict, norm='rare')
if not os.path.exists(basedir+'figures/FigS5_environment_random_forest_0.01_rare.pdf'):
  if not os.path.exists('figures/FigS5_environment_random_forest_0.01_rare'+ext):
    make_env_rf_plot(ft_df_rare, tax_dict_rare, basedir, ASV_dict_rare, meta_dict, mx=0.01, norm='rare')
```

These are shown in Supplementary Sections [2](https://doi.org/10.6084/m9.figshare.12915317) (log, relative abundance and CLR-normalised) and [3](https://doi.org/10.6084/m9.figshare.12233759) (rarefied data).

This is shown as Figure S5 and in [Supplementary Section 3](https://doi.org/10.6084/m9.figshare.12233759).

#### (XVI) Get the metacoder plots for early and late incubation times:
```{python, metacoder_plot, results='hide', fig.keep='all', message=FALSE, cache=TRUE, warning=FALSE}
count = 0
ft_dfs = [ft_df_rare, ft_df_rel_abun, ft_df_log, ft_df_clr]
tax_dicts = [tax_dict_rare, tax_dict_nr, tax_dict_nr, tax_dict_nr]
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  metacoder_py(ft_dfs[count], tax_dicts[count], meta_dict, basedir, norm=norm)
  ft_dfs[count].drop(['lineage'], axis=1, inplace=True)
  count += 1
```

This should save lots of figures in figures/metacoder/, which we will convert from PDF to whichever extension we have been using here:
```{python, convert_pdf_metacoder, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
for norm in ['rare', 'rel_abun', 'log', 'clr']:
  path_to_pdf = basedir+'figures/metacoder/'+norm+'/'
  to_convert = os.listdir(path_to_pdf)
  for fi in to_convert:
    if '.pdf' in fi:
      images = convert_from_path(path_to_pdf+fi, dpi=600)
      ext_caps = 'PNG'
      if ext == '.jpg':
        ext_caps = 'JPEG'
      for image in images:
        image.save(path_to_pdf+fi.replace('.pdf', ext), ext_caps)
```

Make a colorbar:
```{python, metacoder_colorbar, results='hide', fig.keep='all', message=FALSE, eval=FALSE}
make_colorbar_fig(basedir)
```

I have put these plots together outside of this script (just in Powerpoint) and have manually added some labels to the metacoder plot, but this is Figure 6 from the manuscript.

All of these plots are shown in [Supplementary Section 4](https://doi.org/10.6084/m9.figshare.12233762).

#### (XVII) Get Ancom trees and heatmaps for early and late incubation times: 

```{python, ancom_tree_heatmap, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
[tree_heatmap(ft_df_rare, meta_dict, basedir, tax_dict_rare, level=al) for al in [1, 2, 3, 4, 5, 6]]
```

This should save lots of figures in figures/ancom/, which we will convert from PDF to whichever extension we have been using here:
```{python, convert_pdf_png, results='hide', fig.keep='all', message=FALSE, cache=TRUE}
path_to_pdf = basedir+'figures/ancom/'
to_convert = os.listdir(path_to_pdf)
for fi in to_convert:
  if '.pdf' in fi:
    images = convert_from_path(path_to_pdf+fi, dpi=600)
    ext_caps = 'PNG'
    if ext == '.jpg':
      ext_caps = 'JPEG'
    for image in images:
      image.save(path_to_pdf+fi.replace('.pdf', ext), ext_caps)
```

These are all shown in [Supplementary Section 5](https://doi.org/10.6084/m9.figshare.12233765). 

#### (XVIII) Get the separate plots for each study: 
```{python, separate_studies, results='hide', fig.keep='all', message=FALSE, warning=FALSE, cache=TRUE}
if fs_main == 10:
  fs_main -= 2
  fs_small -= 2
  fs_title -= 2
plots_per_study(ft_df_rare, meta_df, meta_dict, w_uf_rare, uw_uf_rare, tax_dict_rare, ASV_dict_rare, basedir, est=est, n_jobs=n_jobs) #get the plots summarising the results per study
```

These are all shown in [Supplementary Section 1](https://doi.org/10.6084/m9.figshare.12233753).