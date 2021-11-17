#!/bin/bash

# 1: directory containing the raw data (after sample_catalog.py) (typically 00-reads/)
DIRECTORY_input="$1"
DIRECTORY=$(realpath $DIRECTORY_input)

# create trimmed folder if it doesnt exist
if [ ! -d 01-quality_control ];
then
    mkdir 01-quality_control
    echo "Created 01-quality_control directory."
else
    echo "Directory 01-quality_control already there."
fi

# create results folder if it doesnt exist

if [ ! -d ../RESULTS ];
then
    mkdir ../RESULTS
    mkdir ../RESULTS/Quality_control
    echo "Created RESULTS and RESULTS/Quality_control directories."
else
    if [ ! -d ../RESULTS/Quality_control ]
    then
        mkdir ../RESULTS/Quality_control
        echo "Created RESULTS/Quality_control directory."
    else
        echo "Directories RESULTS and RESULTS/Quality_control already there"
    fi
fi

printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv

absolute_filepath=$(pwd)

for directory in $DIRECTORY/*/
do
    echo "Starting ${directory} quality control."
    
    # samplename = dirname
    samplename=$(basename $directory)
    
    # make the directories for the quality results: pre-trimming, fastp, post-trimming
    mkdir -p ../RESULTS/Quality_control/${samplename}_quality/{fastqc_pre-trimming_reports,fastp_reports,fastqc_post-trimming}


    fastqc -o ../RESULTS/Quality_control/${samplename}_quality/fastqc_pre-trimming_reports $directory/*.gz
    
    # make the directory for trimmed sequences
    mkdir 01-quality_control/${samplename}

    fastp -i $directory/*R1* -I $directory/*R2* \
    --trim_poly_x \
    --cut_front --cut_tail \
    --cut_mean_quality 20 \
    --length_required 250 \
    --qualified_quality_phred 30 \
    --cut_window_size 4 \
    --html ../RESULTS/Quality_control/${samplename}_quality/fastp_reports/${samplename}_trim_report_fastp.html \
    --json ../RESULTS/Quality_control/${samplename}_quality/fastp_reports/${samplename}_trim_report_fastp.json \
    -o 01-quality_control/${samplename}/${samplename}_trimmed_R1.fq.gz \
    -O 01-quality_control/${samplename}/${samplename}_trimmed_R2.fq.gz

    # second quality analysis
    fastqc -o ../RESULTS/${samplename}_quality/fastqc_post-trimming 01-quality_control/${samplename}/*_trimmed_*
    
    # add the data to the sample manifest
    printf "${samplename}\t${absolute_filepath}/01-quality_control/${samplename}/${samplename}_trimmed_R1.fq.gz\t${absolute_filepath}/01-quality_control/${samplename}/${samplename}_trimmed_R2.fq.gz\n" >> manifest.tsv
    echo "${directory} quality control finished."

done

echo "Launching MultiQC"

mkdir ../RESULTS/multiqc_results
multiqc ../RESULTS/*/*/* -o ../RESULTS/multiqc_results

echo "MultiQC ended"
# Ready to launch multiqc on RESULTS (use the correct environment for it)
# multiqc RESULTS/*
