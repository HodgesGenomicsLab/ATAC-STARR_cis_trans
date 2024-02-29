[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10728132.svg)](https://doi.org/10.5281/zenodo.10728132)
# Divergence in both cis and trans drives human gene regulatory evolution  
This is the code repository for the analyses performed in Hansen, T., Fong, S., et al. bioRxiv 2022. Any questions or comments, please email Emily Hodges at emily.hodges(at)vanderbilt.edu. 

Please also visit our ATAC-STARR-seq github and read our ATAC-STARR-seq method publication. 
- https://github.com/HodgesGenomicsLab/ATAC-STARR-seq   
- Hansen TJ, Hodges E. Genome Res. 2022 (https://doi.org/10.1101/gr.276766.122)

## FASTQ Proccessing
Files were processed in the same manner as described in the github (https://github.com/HodgesGenomicsLab/ATAC-STARR-seq). For replicate 3 in all conditions, we used the _fastq_processing.py_ script from the ATAC-STARR-seq github.  For repicates 1 & 2, we used custom bash scripts because they were processed before the _fastq_processing.py_ script had been written. Both scripts follow the same process outlined in the manuscript. The merge_runs.slrm script is how we merged runs 1 and 2 for LGR2 and GLR1 samples. 

## Region Calling
This section contains code for calling all of our region sets, including the cis & trans regions. 

- call_activity.ipynb details how we called accessible peaks and active regions for each condition. This also includes important supplementary analyses.
- HH-vs-MM.ipynb details how we called conserved active and species-specific active regions.
- cis-or-trans.ipynb details how we called cis effects and trans effects as well as the cis only, trans only, and cis & trans groups.
- shuffle_analysis.ipynb details the observed vs. expected overlap analyses.
- activity_bigwigs.ipynb details how we created the ATAC-STARR-seq bigwigs.

## Functional Characterization of Cis and Trans Effects
This section contains the code used for the "functional" analyses described in the manuscript. 

## TF Footprinting
This section contains the code used for TF footprinting-related analyses. 

## Evolutionary Characterization of Cis and Trans Effects
This section contains the code used for the "evolutionary" analyses described in the manuscript. This folder is subdiveded into code used for conservation, acceleration, and transposable element analyses. 

## Human Variant Enrichment Analysis
This section contains the code used for the human variant enrichment analyses described in the manuscript. This folder is subdiveded into code used for eQTL and UKBB GWAS trait analyses. 

## Gene Expression Analysis
This section contains the code used for the gene expression analyses performed in the manuscript. 

## Miscellaneous 
This section contains miscellaneous code used in the study. 
