# Neurospora_crassa_transcriptomics
Characterizing the gene-environment interaction underlying natural morphological variation in Neurospora crassa conidiophores using high-throughput phenomics and transcriptomics

### Quantification:
The quantificaiton scripts are written for the alignment of short read illumina RNA sequencing data. The data was quantified using two standard pipelines. The counts and differentially expressed transcripts were cross validated to ensure the reproducibility of results. The genome annotation was built on GCF_000182925.2_NC12_genomic.fna.gz 

CLI:
```
nohup ./star_quantification.sh &
```


### Analysis:
A simple differential expression analysis was conducted using the quantifications.tab files. The transcript counts were analyzed in R with canaonical transcriptomics softwares from bioconductor. 

### Genetic model: 
Table of MATLAB scripts to calculate varied inheritance models for the discrete complex trait of conidiophore type by the Method of Maximum likelihood using iteratively reweighted least squares.  The first 8 columns taken from Table 6 is some of the output of the MATLAB script.  The last column is the name of the corresponding MATLAB script.  Each script is self contained and can simply be run within MATLAB by hitting run once the script is clicked on and opened.



