#!/bin/bash
# author: Michael Skaro
# date: 2/22/2021
# Purpose: Salmon quant using indicies, check for discrepencies with star

for i in `cat /home/mskaro1/storage/neurospora_crassa/fastq/ID`;
do
        echo $i\ ;
#       echo $i\_R1_trimmed_001.fastq.gz ;
#       echo $i\_R2_trimmed_001.fastq.gz ;

        salmon quant -i /home/mskaro1/storage/neurospora_crassa/Annotation/salmon_index --libType A \
          -1 $i\_L001_R1_001.fastq.gz \
          -2 $i\_L001_R2_001.fastq.gz \
          -p 8 --validateMappings \
          -o quant/$i;


done
