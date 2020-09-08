# QIIME2 analysis

Running through the code here will allow you to go from raw reads (*i.e.* from a MiSeq sequencing run) to the files found [on FigShare here](https://doi.org/10.6084/m9.figshare.12227522).

If you want to re-run all of the analysis for all studies, then start at the beginning and download [the raw reads files from here](https://doi.org/10.6084/m9.figshare.12931372). If you just want to add an additional study to the studies already included, then you can just download [the merged QIIME2 objects from here](https://doi.org/10.6084/m9.figshare.12217682.v3). You can then skip ahead to section 2.

## Setup

I used 12 threads on an Ubuntu server for most of my analyses. I think some parts will probably struggle with so many samples if you try to do this locally, but I haven't tried this so can't say for sure. If you have fewer threads available then you should change the '12' in parts to be whatever you have available.
This follows the Microbiome helper tutorial [here](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2020.2)).
You can view all of the summary files (when you use the summarise commands) [at the QIIME2 view website here](https://view.qiime2.org/).

Activate the QIIME2 environment (if you do not already have this installed then follow the instructions [here](https://docs.qiime2.org/2020.6/install/):

```bash
conda activate qiime2-2019.10
```

## 1. Method for individual studies

Note that this uses [Deblur](https://github.com/biocore/deblur). [DADA2](https://benjjneb.github.io/dada2/) could also be used, but given that we didn't know whether all samples from each study came from the same sequencing run, we chose the per sample denoising approach of Deblur. 

This assumes that your .fastq.gz files are in a directory called ```raw_data```. If you are just downloaded the raw files from FigShare then you can skip ahead to (III) (but note that rather than just ```reads.qza``` these files are in the format of ```StudyName_reads.qza```).

### (I) Run quality checks on the data:
```bash
mkdir fastqc_out
fastqc -t 12 raw_data/*.fastq.gz -o fastqc_out
multiqc fastqc_out/
```
You should now look at the ```multiqc_report.html``` to ensure that everything is high enough quality to proceed.

### (II) Import your data to the QIIME2 format:
```bash
qiime tools import \
            --type SampleData[PairedEndSequencesWithQuality] \
            --input-path raw_data/ \
            --output-path reads.qza \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
```

### (III) Trim primers (if present) with cutadapt. 
The primer sequences shown here are for 341F and 802R - these will need changing if you have used different primers:
```bash
qiime cutadapt trim-paired \
            --i-demultiplexed-sequences reads.qza \
            --p-cores 12 \
            --p-front-f CCTACGGGNGGCWGCAG \
            --p-front-r GACTACHVGGGTATCTAATCC \
            --p-discard-untrimmed \
            --p-no-indels \
            --o-trimmed-sequences reads_trimmed.qza
```

### (IV) Summarize the trimmed files:
```bash
qiime demux summarize \
            --i-data reads_trimmed.qza \
            --o-visualization reads_trimmed_summary.qzv
```

### (V) Join paired ends 
If the reads were already trimmed then just use reads.qza as the input here:
```bash
qiime vsearch join-pairs \
            --i-demultiplexed-seqs reads_trimmed.qza \
            --o-joined-sequences reads_joined.qza
```

### (VI) Summarize the joined pairs 
If too many reads were removed then you may need to play around with some of the other options at [here](https://docs.qiime2.org/2020.2/plugins/available/vsearch/join-pairs/):
```bash
qiime demux summarize \
            --i-data reads_joined.qza \
            --o-visualization reads_joined_summary.qzv
```

### (VII) Filter out low quality reads:
```bash
qiime quality-filter q-score-joined \
            --i-demux reads_joined.qza \
            --o-filter-stats filt_stats.qza \
            --o-filtered-sequences reads_joined_filtered.qza
```

### (VIII) Summarize these reads (and look at where to trim) 
You should look at the positions where the quality starts to drop below 30 and use these as trim lengths:
```bash
qiime demux summarize \
            --i-data reads_joined_filtered.qza \
            --o-visualization reads_joined_filtered_summary.qzv
```

### (IX) Run deblur 
You can remove the --p-left-trim-len if you don't need to remove any from this end:
```bash
qiime deblur denoise-16S \
            --i-demultiplexed-seqs reads_joined_filtered.qza \
            --p-trim-length 402 \
            --p-left-trim-len 0 \
            --p-sample-stats \
            --p-jobs-to-start 12 \
            --p-min-reads 1 \
            --output-dir deblur_output_quality
```

### (X) Summarize the feature table to see how many reads we now have:
```bash
qiime feature-table summarize \
            --i-table deblur_output_quality/table.qza  \
            --o-visualization deblur_table_summary.qzv
```

## 2. Merge studies

Once you have performed all of the steps in ```1.``` for all studies that you are including, you can merge the studies into a representative sequences and a table object. 

### (I) Rename the current representative sequences and merged tables:
(Do this if you are just adding a study to the current studies used)
```bash
mv merged_representative_sequences.qza previous_merged_representative_sequences.qza
mv merged_table.qza previous_merged_table.qza
```

### (II) Combine feature tables:
You will need to replace 'your_folder_name' with the folder that contains your tables to be added (when I did this for all studies, I just added additional --i-tables table_name.qza lines):
```bash
qiime feature-table merge \
            --i-tables your_folder_name/deblur_output_quality/table.qza \
            --i-tables previous_merged_table.qza \
            --o-merged-table merged_table.qza
```

### (III) Combine the sequence objects:
You will again need to replace 'your_folder_name' with the folder that contains your sequences to be added (when I did this for all studies, I just added additional --i-data representative_sequences_name.qza lines):
```bash
qiime feature-table merge-seqs \
            --i-data your_folder_name/deblur_output/representative_sequences.qza \
            --i-data previous_merged_representative_sequences.qza \
            --o-merged-data merged_representative_sequences.qza
```

## 3. Combined processing

Now that all of the samples that we are looking at are combined into the merged sequences and table files, we can classify and analyze them.

### (I) Summarize the combined feature tables (this is to check that everything looks OK after the merges, and can be skipped if not necessary):
```bash
qiime feature-table summarize \
            --i-table merged_table.qza  \
            --o-visualization merged_table_summary.qzv
```

### (II) Classify the features (this part will probably take the longest - it may take at least a day or so and is the part that may not be possible on a local computer):
```bash
qiime feature-classifier classify-sklearn \
            --i-reads merged_representative_sequences.qza \
            --i-classifier ref_alignments/classifier_silva_132_99_16S.qza \
            --p-n-jobs 12 \
            --output-dir taxa
```
As these sequences come from different 16S regions, I downloaded the full length 16S classifier from [here](https://docs.qiime2.org/2020.6/data-resources/). There is now an updated SILVA version, but I used the Silva 132 classifier (this can only improve upon classification accuracy, so I recommend using the latest one).

### (III) Export this file to look at the classifications:
```bash
qiime tools export \
            --input-path taxa/classification.qza \
            --output-path taxa
```

### (IV) Filter low abundance features:
```bash
qiime feature-table filter-features \
            --i-table merged_table.qza \
            --p-min-frequency 10 \
            --p-min-samples 1 \
            --o-filtered-table merged_table_filtered.qza
```

### (V) Filter potential contaminants and those not classified at the kingdom level:
```bash
qiime taxa filter-table \
            --i-table merged_table_filtered.qza \
            --i-taxonomy taxa/classification.qza \
            --p-include D_1__ \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table merged_table_filtered_contamination.qza
```

### (VI) Summarize the filtered table:
```bash
qiime feature-table summarize \
            --i-table merged_table_filtered_contamination.qza \
            --o-visualization merged_table_filtered_contamination_summary.qzv
```
Now find out how many features you have as well as the maximum sample depth (this is the "Maximum Frequency" in the "Frequency per sample" section).

### (VII) Obtain rarefaction curves for samples:
```bash
qiime diversity alpha-rarefaction \
            --i-table merged_table_filtered_contamination.qza \
            --p-max-depth 995391 \
            --p-steps 20 \
            --p-metrics 'observed_otus' \
            --o-visualization merged_rarefaction_curves.qzv
```

### (VIII) Filter samples that have below 2000 reads:
```bash
qiime feature-table filter-samples \
            --i-table merged_table_filtered_contamination.qza \
            --p-min-frequency 2000 \
            --o-filtered-table  merged_table_final.qza
```

### (IX) Rarefy remaining samples to 2000:
```bash
qiime feature-table rarefy \
            --i-table merged_table_final.qza \
            --p-sampling-depth 2000 \
            --o-rarefied-table merged_table_final_rarefied.qza
```

### (X) Filter the sequences to contain only those that are in the rarefied feature table:
```bash
qiime feature-table filter-seqs \
            --i-data merged_representative_sequences.qza \
            --i-table merged_table_final_rarefied.qza \
            --o-filtered-data  representative_sequences_final_rarefied.qza
```

### (XI) Export feature table and sequences for rarefied data:
```bash
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

### (XII) Obtain phylogenetic tree using SEPP fragment insertion and the silva reference database for rarefied data:
```bash
qiime fragment-insertion sepp \
            --i-representative-sequences representative_sequences_final_rarefied.qza \
            --i-reference-database ref_alignments/sepp-refs-silva-128.qza \
            --o-tree insertion_tree_rarefied.qza \
            --o-placements insertion_placements_rarefied.qza \
            --p-threads 12
```
You can download the reference file [here](https://docs.qiime2.org/2020.6/data-resources/). At the time of writing, this still used Silva 128, but I would recommend using an updated version if there is one.

### (XIII) Export the resulting insertion tree:
```bash
qiime tools export \
            --input-path insertion_tree_rarefied.qza \
            --output-path exports
```

### (XIV) The files inside the exports folder should then be copied to the folder that the subsequent analyses will be carried out in, e.g.:
```bash
for i in exports/* ; cp $i paper_data/qiime_output/; done
```

And rename the files:
```bash
mv paper_data/qiime_output/feature-table_w_tax.txt paper_data/qiime_output/feature-table_w_tax_rare.txt
mv paper_data/qiime_output/dna-sequences.fasta paper_data/qiime_output/dna-sequences_rare.fasta
mv paper_data/qiime_output/tree.nwk paper_data/qiime_output/tree_rare.nwk
```

### Now do the same for the non-rarefied data

### (XV) Move the data to a new folder and summarize the number of features etc:
```bash
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
```bash
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

### (XVI) Export feature table and sequences for not rarefied data (before and after filtering, because we still want to know how many sequences in each sample to calculate relative abundance):
```bash
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

### (XVI) Now put these sequences into the tree:
```bash
qiime fragment-insertion sepp \
            --i-representative-sequences representative_sequences_final_filtered.qza \
            --i-reference-database /home/robyn/other/ref_alignments/sepp-refs-silva-128.qza \
            --o-tree insertion_tree_not_norm.qza \
            --o-placements insertion_placements_not_norm.qza \
            --p-threads 24
```

### (XVII) Export the tree:
```bash
qiime tools export \
           --input-path insertion_tree_not_norm.qza \
           --output-path exports
```

### (XVIII) The files inside the exports folder should then be copied to the folder that the subsequent analyses will be carried out in, e.g.:
```bash
for i in exports/* ; cp $i paper_data/qiime_output/; done
```

And rename the files:
```bash
mv paper_data/qiime_output/feature-table_w_tax.txt paper_data/qiime_output/feature-table_w_tax_not_rare.txt
mv paper_data/qiime_output/feature-table_w_tax_filtered.txt paper_data/qiime_output/feature-table_w_tax_not_rare_filtered.txt
mv paper_data/qiime_output/dna-sequences.fasta paper_data/qiime_output/dna-sequences_not_rare.fasta
mv paper_data/qiime_output/tree.nwk paper_data/qiime_output/tree_not_rare.nwk
```