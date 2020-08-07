# Details on re-running the Plastisphere meta-analysis
---
### Wright, Langille and Walker (submitted) Food or just a free ride? A meta-analysis reveals the global diversity of the Plastisphere

**I am currently in the process of creating an R notebook that includes all of the analysis presented in the paper, but for now, everything is here:**

**There are several Jupyter notebooks that are intended for re-running the analyses presented in the Plastisphere meta-analysis. All of these can be re-run, or you can choose to re-run any part of this. The basic steps used here are:**
1. **QIIME2 analysis separate per study (includes running quality checks on data, importing reads to QIIME2 format, trimming primers, joining paired end sequences, filtering out low quality reads and running deblur)**<br/>
    - You should only do this if you wish to add your own study to the 35 studies that we have used, and you can do this using the notebook at metaanalysis_run_QIIME2/separate_studies.ipynb<br/>
    - If you do this, then you should then merge the results from your study with those in [this file](https://doi.org/10.6084/m9.figshare.12217682) in the QIIME_input folder. You can do this using the notebook at metaanalysis_run_QIIME2/merge_studies.ipynb
    
2. **QIIME2 analysis on all merged studies. You can either do this with the new files created in step 1 (if you are adding your own study to these results) or with the [these files](https://doi.org/10.6084/m9.figshare.12217682) in the QIIME_input folder**<br/>
    - You can do this using the notebook at metaanalysis_run_QIIME2/QIIME2_analysis.ipynb <br/>
    
3. **ONLY if you want to add additional genes to the PICRUSt2 reference genomes - instructions for adding to this database can be found in the notebook PICRUST2_run_HMM/add_gene_to_picrust.ipynb**

After this stage, we have a .fasta file with all representative sequences, a feature table (.txt) containing the abundances for all ASVs and an insertion tree. <br/><br/>
If you do not want to re-run the sequence processing or add in additional studies, these can be [downloaded from here](https://doi.org/10.6084/m9.figshare.12227303). The files at this link contain all files produced by running this script, as well as an almost empty directory that has only the three files ('dna-sequences.fasta', 'feature-table_w_tax.txt' and 'tree.nwk') and the picrust folder containing 'ko_all.txt' and 'kegg_list.csv'. If you carried out step 3 to add genes to the PICRUSt2 reference, then replace 'ko_all.txt' with the file that you created.

**If you have added a study then you should also add details of this study to three files that are in the 'python_analysis_20-04-14' folder:**<br/>
    - Study_dates.csv<br/>
    - Study_location.csv<br/>
    - metadata.txt<br/>

4. **Carry out all subsequent analyses (basic summaries, random forests, ANCOM, metacoder, PICRUSt2, etc.)**
    - You can do this using the notebook python_analysis_20-04-14/Plastisphere_metaanalysis.ipynb and other scripts contained in the 'python_analysis_20-04-14' folder
