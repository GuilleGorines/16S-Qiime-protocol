#!/bin/bash

# 1: directory containing the raw data (after Sample_catalog.py) (typically RAW/)
DIRECTORY_input="$1"
DIRECTORY=$(realpath $DIRECTORY_input)

# 2: directory on which the full analysis is performed (typically .) 
MOTHER_DIR_input="$2"
MOTHER_DIR=$(realpath $MOTHER_DIR_input)

# create trimmed folder if it doesnt exist
if [ ! -d ${MOTHER_DIR}/TRIMMED ];
then
    mkdir ${MOTHER_DIR}/TRIMMED
fi

# create analysis folder if it doesnt exist
if [ ! -d ${MOTHER_DIR}/ANALYSIS ];
then
    mkdir ${MOTHER_DIR}/ANALYSIS
fi

printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > ${MOTHER_DIR}/ANALYSIS/manifest.tsv

# create results folder if it doesnt exist
if [ ! -d ${MOTHER_DIR}/RESULTS ];
then
    mkdir ${MOTHER_DIR}/RESULTS
fi

for directory in $DIRECTORY/*/
do
    # samplename = dirname
    samplename=$(basename $directory)
    
    # first quality result
    cd $directory
    mkdir ${MOTHER_DIR}/RESULTS/${samplename}_quality
    fastqc -o ${MOTHER_DIR}/RESULTS/${samplename}_quality *.gz
    
    mkdir ${MOTHER_DIR}/RESULTS/${samplename}_trimming
    # perform trimming
    fastp -i *R1* -I *R2* \
    --trim_poly_x \
    --cut_front --cut_tail \
    --cut_mean_quality 20 \
    --length_required 250 \
    --qualified_quality_phred 30 \
    --cut_window_size 4 \
    --html  ${MOTHER_DIR}/RESULTS/${samplename}_trimming/${samplename}_trim_report.html \
    --json  ${MOTHER_DIR}/RESULTS/${samplename}_trimming/${samplename}_trim_report.json \
    -o ${samplename}_trimmed_R1.fq.gz \
    -O ${samplename}_trimmed_R2.fq.gz  

    # second quality analysis
    fastqc -o ${MOTHER_DIR}/RESULTS/${samplename}_quality *_trimmed_*
    mv ${samplename}_trimmed_R1.fq.gz $MOTHER_DIR/TRIMMED
    mv ${samplename}_trimmed_R2.fq.gz $MOTHER_DIR/TRIMMED

    mv ${samplename}_trim.html $MOTHER_DIR/RESULTS 
     
    # add the data to the sample manifest
    printf "${samplename}\t${MOTHER_DIR}/TRIMMED/${samplename}_trimmed_R1.fq.gz\t${MOTHER_DIR}/TRIMMED/${samplename}_trimmed_R2.fq.gz\n" >> ${MOTHER_DIR}/ANALYSIS/manifest.tsv
done

# Ready to launch multiqc on RESULTS