# Food or just a free ride? A meta-analysis reveals the global diversity of the Plastisphere


This repository contains instructions for reproducing the meta-analysis of Plastisphere studies performed in:
**Food or just a free ride? A meta-analysis reveals the global diversity of the Plastisphere** (currently under revision in the ISME journal)</br>
You can do this with only the studies that are included in the meta-analysis, or you can include new studies too (in order to evaluate how they fit with previous studies as well as whether they change any of the conclusions on the composition of the Plastisphere).

Please contact [Robyn Wright](mailto:robyn.wright@dal.ca) with any questions.

All of the files that you need to reproduce the analyses presented in the meta-analysis are in the [metaanalysis_files folder](https://github.com/R-Wright-1/Plastisphere-MetaAnalysis/tree/master/metaanalysis_files). Some of the larger files have been zipped to allow them to be uploaded here, so you should unzip them before trying to run anything.
To unzip all files (in command line):
```
for i in paper_data/*.tar.bz2 ; do tar -jxvf $i ; done
for i in paper_data/agglom/*.tar.bz2 ; do tar -jxvf $i ; done
for i in paper_data/qiime_output/*.tar.bz2 ; do tar -jxvf $i ; done
```
You can also find them (unzipped) in [this Figshare file](https://doi.org/10.6084/m9.figshare.12227303). Note that the figures can only be found in the Figshare file as they were too large to upload to Github.

If you just want to read through the analyses, then you can find the .html file either [in this Figshare file](https://doi.org/10.6084/m9.figshare.12923855) or you can unzip the .html file that is in the metaanalysis_files folder. 

If you want to repeat the analysis (with or without the inclusion of additional studies) then you can follow the instructions in the[Rmarkdown](https://rmarkdown.rstudio.com/) and html files. If you are unfamiliar with Rmarkdown files, they allow you to run 'chunks' of code written in different languages within them, and this has therefore allowed all of the code used for the analysis to be included in a single document. (They also allow you to just run everything all at once, generating a .html output file, which theoretically should work as long as you have the files and packages that it is expecting, but we all know that that's not how science *actually* works). 

Older analysis scripts (in the [Older folder](https://github.com/R-Wright-1/Plastisphere-MetaAnalysis/tree/master/Older)) used Jupyter notebooks or Python scripts to call R scripts. You can find these in the folder. Feel free to ask me questions if you'd like to use them, but I suggest that using the Rmarkdown file will be much easier. 
