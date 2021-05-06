#!/bin/bash

DIRECTORIES="$1"
LOCATION=$(pwd)
SAMPLE_MANIFEST="${LOCATION}/manifest.csv"


for directory in $DIRECTORIES/*/;
do
    cd $directory
    samplename=$(basename $directory)
    mkdir fastqc_result_raw
    fastqc -o fastqc_result_raw *.gz
    
    mkdir raw
    mv *.gz raw/

    fastp -i raw/*R1* -I raw/*R2*\
    -o "${samplename}_trimmed_R1.fq.gz"\
    -O "${samplename}_trimmed_R2.fq.gz"
    
    dir_forward=$(realpath *_trimmed_R1.fq.gz)
    dir_reverse=$(realpath *_trimmed_R2.fq.gz)

    printf "${samplename},${dir_forward},forward\n" >> $SAMPLE_MANIFEST
    printf "${samplename},${dir_reverse},reverse\n" >> $SAMPLE_MANIFEST

    mkdir fastqc_result_trimmed
    fastqc -o fastqc_result_trimmed *_trimmed_R*

done
