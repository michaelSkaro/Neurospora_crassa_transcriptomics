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

| Model                | Χ^2    | df | p       | X^2 H0 − X^2 HA          | rho for HA vs. H0 | Alternative HA                                                                                       | MATLAB scripts                                     |
|----------------------|--------|----|---------|--------------------------|-------------------|------------------------------------------------------------------------------------------------------|----------------------------------------------------|
| Full epistatic       | 0      | 0  | -       | -                        | -                 | -                                                                                                    | estimator_inh_model_full_espistatic_loglinear_V8.m |
| αβ=0                 | 13.37  | 1  | <0.0001 | 13.37 – 0.00 = 13.37     | <0.0001           | HA = full epistatic                                                                                  | estimator_inh_model_ab_loglinear_V8.m              |
| αβ=α=0               | 16.27  | 2  | <0.0001 | 16.27 – 13.37 = 9.12     | 0.002             | HA = αβ=0                                                                                            | estimator_inh_model_ab_a_loglinear_V8.m            |
| βγ=0                 | 43.42  | 1  | <0.0001 | 43.42 – 0.00 = 43.42     | <0.0001           | H0 = full epistatic                                                                                  | Estimator_inh_model_bc_loglinear_V8.m              |
| βγ=β=0               | 43.21  | 2  | <0.0001 | 43.21 – 43.42            | 0.64              | HA = βγ=0                                                                                            | estimator_inh_model_bc_b_loglinear_V8.m            |
| αβ=β=0               | 29.88  | 2  | <0.0001 | 29.88 – 13.37 = 16.51    | <0.0001           | HA = αβ=0                                                                                            | estimator_inh_model_ab_a_loglinear_V8.m            |
| αγ=0                 | 11.11  | 1  | 0.0009  | 11.11 – 0.00 = 11.11     | 0.0009            | HA = full epistatic                                                                                  | estimator_inh_model_ac_loglinear_V8.m              |
| αγ=γ=0               | 126.37 | 2  | <0.0001 | 126.37 – 11.11 = 115.26  | <0.0001           | HA = αγ=0                                                                                            | estimator_inh_model_ac_c_loglinear_V8.m            |
| αβ=βγ=αγ=0, additive | 44.61  | 3  | <0.0001 | 44.61 – 13.37 = 31.24    | <0.0001           | HA = αβ=0                                                                                            | estimator_inh_model_additive_fixed_sizes_V8.m      |
| environmental        | 84.27  | 6  | <0.0001 | 84.27 – 44.61 = 39.66    | <0.0001           | HA = additive, H2 = (84.27 – 44.61)/84.27 = 0.47, H0 = environmental model, H1 = full additive model | Estimator_inh_model_environmental_fixed_sizes_V8.m |


TO RUN:
```
$ alias matlab='/Applications/{YOUR MATLAB}/bin/matlab -nodesktop -nosplash $*'

$ matlab
```

MATLAB: 

```Matlab

>>> estimator_inh_model_full_espistatic_loglinear_V8

```





