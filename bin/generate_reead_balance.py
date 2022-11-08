'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 0
CREATED: 13-1-2022
REVISED: 07-11-2022
DESCRIPTION: 
    Ad-hoc script to get the read balance in a 16S protocol
INPUT:
    1. multiqc_data.json generated by multiqc data
OUTPUT:
    Balance of reads obtained in the whole service 
'''

import json
import sys

multiqc_json = "multiqc_data/multiqc_data.json"
stats_file = "stats.tsv"

def clean_name(dirty_name):
    clean_name = dirty_name.replace("_R1","").replace("_R2","").replace("_adapter_removed","").replace("_filtered","").replace("_trimmed","").replace("_fastp", "")
    return clean_name

with open(multiqc_json,"r") as infile:
    mqc_json = json.load(infile)

read_data_dict = {}

# PART 1: parse multiqc data json file

# Find all samplenames from the beginning
# Initialize all values so they are easier to work with
for program, sample_list in mqc_json["report_data_sources"].items():
    for item in sample_list["all_sections"].keys():

        digested_sample_name = clean_name(item)

        if digested_sample_name not in read_data_dict.keys():
            read_data_dict[digested_sample_name] = {}

            read_data_dict[digested_sample_name]["Sample name"] = digested_sample_name
            read_data_dict[digested_sample_name]["RAW reads"] = ["-"]
            read_data_dict[digested_sample_name]["Reads AFTER Cutadapt"] = ["-"]
            read_data_dict[digested_sample_name]["Total reads AFTER trimming"] = ["-"]
            read_data_dict[digested_sample_name]["R1 mean length AFTER trimming"] = ["-"]
            read_data_dict[digested_sample_name]["R2 mean length AFTER trimming"] = ["-"]
            read_data_dict[digested_sample_name]["R1 mean length BEFORE trimming"] = ["-"]
            read_data_dict[digested_sample_name]["R2 mean length BEFORE trimming"] = ["-"]

            read_data_dict[digested_sample_name]["Reads after Qiime2 filtering"] = ["-"]
            read_data_dict[digested_sample_name]["Reads after Qiime2 denoising"] = ["-"]
            read_data_dict[digested_sample_name]["Reads merged by Qiime2"] = ["-"]
            read_data_dict[digested_sample_name]["Non-chimeric merged reads"] = ["-"]


# Find the dicts with data inside "report_saved_raw_data"
for program, data in mqc_json["report_saved_raw_data"].items():

    if program == "multiqc_cutadapt":
        for sample_name, cutadapt_dict in data.items():
            digested_sample_name = clean_name(sample_name)
            read_data_dict[digested_sample_name]["RAW reads"] = int(cutadapt_dict["pairs_processed"])
            read_data_dict[digested_sample_name]["Reads AFTER Cutadapt"] = int(cutadapt_dict["pairs_written"])

    elif program == "multiqc_fastp":
        for sample_name, fastp_dict in data.items():
            digested_sample_name = clean_name(sample_name)
            read_data_dict[digested_sample_name]["Total reads AFTER trimming"] = int(fastp_dict["summary"]["after_filtering"]["total_reads"])/2
            read_data_dict[digested_sample_name]["R1 mean length BEFORE trimming"] = int(fastp_dict["summary"]["before_filtering"]["read1_mean_length"])
            read_data_dict[digested_sample_name]["R2 mean length BEFORE trimming"] = int(fastp_dict["summary"]["before_filtering"]["read2_mean_length"])
            read_data_dict[digested_sample_name]["R1 mean length AFTER trimming"] = int(fastp_dict["summary"]["after_filtering"]["read1_mean_length"])
            read_data_dict[digested_sample_name]["R2 mean length AFTER trimming"] = int(fastp_dict["summary"]["after_filtering"]["read2_mean_length"])        

# PART 2: qiime2 stats file
with open(stats_file,"r") as statfile:
    statfile = statfile.readlines()
    statfile = [item.split("\t") for item in statfile[2:]]

# item 0 will always (theoretically at least) match the digested_sample_name
for item in statfile:
    read_data_dict[item[0]]["Reads after Qiime2 filtering"] = int(item[2])
    read_data_dict[item[0]]["Reads after Qiime2 denoising"] = int(item[4])
    read_data_dict[item[0]]["Reads merged by Qiime2"] = int(item[5])
    read_data_dict[item[0]]["Non-chimeric merged reads"] = int(item[7])

# 0: Sample id
# 1: Raw reads
# 2: Reads after Cutadapt
# 3: Reads after Trimming
# 4: R1 Raw mean length
# 5: R2 Raw mean length
# 6: R1 Trimmed mean length
# 7: R2 Trimmed mean length
# 8: Reads after Qiime2 filtering
# 9: Reads after Qiime2 denoising
# 10: Reads merged by Qiime2
# 11: Non-chimeric merged reads

with open("Read_balance.tsv","w") as outfile:
    headers_list = ["Sample id", "RAW Reads",
                    "Reads after Cutadapt","% of reads after Cutadapt",
                    "Reads after Trimming","% of reads after Trimming",
                    "R1 reads RAW mean length","R2 reads RAW mean length",
                    "R1 reads Trimmed mean length","R2 reads Trimmed mean length",
                    "Reads after Qiime2 filtering","% of reads after Qiime2 filtering",
                    "Reads after Qiime2 denoising","% of reads after Qiime2 denoising",
                    "Reads merged by Qiime2","% of merged reads after Qiime2",
                    "Non-chimeric merged reads","% of non-chimeric merged reads"]
    outfile.write("\t".join(headers_list))
    outfile.write("\n")

    for line in read_data_dict.values():

        samplename = line["Sample name"]
        raw = line["RAW reads"]

        after_cutadapt = line["Reads AFTER Cutadapt"]
        after_cutadapt_percentage = after_cutadapt*100/raw

        # trimmomatic counts single reads instead of pairs
        after_trimmomatic = line["Total reads AFTER trimming"]
        after_trimmomatic_percentage = after_trimmomatic*100/raw

        R1_raw = line["R1 mean length BEFORE trimming"]
        R2_raw = line["R2 mean length BEFORE trimming"]
        R1_trim = line["R1 mean length AFTER trimming"]
        R2_trim = line["R2 mean length AFTER trimming"]

        qiime_filtering = line["Reads after Qiime2 filtering"]
        qiime_filtering_percentage = qiime_filtering*100/raw

        qiime_denoising = line["Reads after Qiime2 denoising"]
        qiime_denoising_percentage = qiime_denoising*100/raw

        qiime_merged = line["Reads merged by Qiime2"]
        qiime_merged_percentage = qiime_merged*100/raw

        qiime_non_chimeric = line["Non-chimeric merged reads"]
        qiime_non_chimeric_percentage = qiime_non_chimeric*100/raw

        data_list = [samplename, raw, after_cutadapt, after_cutadapt_percentage,
                     after_trimmomatic, after_trimmomatic_percentage, R1_raw,
                     R2_raw, R1_trim, R2_trim, qiime_filtering, qiime_filtering_percentage,
                     qiime_denoising, qiime_denoising_percentage, qiime_merged,
                     qiime_merged_percentage, qiime_non_chimeric, 
                     qiime_non_chimeric_percentage]

        data_list = map(lambda x: str(x).replace(".",","),data_list)


        outfile.write("\t".join(data_list))