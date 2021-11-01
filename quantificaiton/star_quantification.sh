#!/bin/bash
# author: Michael Skaro
# date: 2/22/2021
# Purpose: star quant using indicies

for i in `cat /home/mskaro1/storage/neurospora_crassa/fastq/ID`;
do
        echo $i\ ;
        echo  $i\_R1_001.fastq.gz \ ;
        echo  $i\_R2_001.fastq.gz \ ;

        STAR --genomeDir /home/mskaro1/storage/neurospora_crassa/Annotation/ \
        --readFilesIn $i\_L001_R1_001.fastq.gz $i\_L001_R2_001.fastq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix $i \
        --readFilesCommand zcat \
        --runThreadN 16 \
        --twopassMode Basic


done
