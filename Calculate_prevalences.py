'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 0
CREATED: 15-3-2022
REVISED: 15-3-2022
DESCRIPTION: 
    Ad-hoc script to get the list of present organsims by group in the metadata
INPUT:
    1. The raw counts (lvl 7)
    2. the metadata file
OUTPUT:
    -Json with all presence by level and category
    -Directory for each category
        -tsv file for each level (1-7) in that category 
'''
import sys
import os

import pandas as pd
from qiime2 import Artifact

qza_in = sys.argv[1]
metadata_file = sys.argv[2]
lvl = sys.argv[3]

def normalize_dataframe(dataframe, criteria=0):
    """
    Change the dataframe to an absence-presence matrix
    based on a criteria (by now, a number)
    """
    
    row_number, col_number = dataframe.shape
    
    for row in range(0, row_number):
        for col in range(0, col_number):
            if dataframe.iloc[row, col] >= criteria:
                dataframe.iloc[row, col] = 1
            else:
                dataframe.iloc[row, col] = 0
                
    return dataframe

def create_category_dict(metadata):
    """
    Create, from the metadata dataframe, a dict with
    key: category; val: values in that category
    if only one category, it wont be taken into account
    """
    valid_categories = dict()
    category_names_list = list(metadata.columns)

    # get all different possibilities for each metadata column
    for col_index in range(metadata.shape[1]):
        
        # list from a set to avoid repeating
        groups = (list(set(metadata[category_names_list[col_index]])))
        
        # if more than 1 different category, add it to the dict
        if len(groups) > 1:
            category_name = category_names_list[col_index]
            valid_categories[category_name] = [item for item in groups]

    return valid_categories, category_names_list

# Create the output directory
os.mkdir(f"lvl_{lvl}/prevalence")

# Import the counts
qza = Artifact.load(qza_in)
counts = qza.view(pd.DataFrame)

# Import the metadata file
metadata = pd.read_csv(
    metadata_file,
    sep='\t',
    header=0,
    index_col=0
    )

# Generate the full dataframe
# Concat metadata and counts
full_df = pd.concat([metadata, counts], axis=1)

# Get the valid categories
valid_categories, category_names_list = create_category_dict(metadata)

# Generate a sub-df for each category
for category, values in valid_categories.items():
    for value in values:
        # Drop metadata columns
        sub_df = full_df[full_df[category] == value].drop(category_names_list, axis=1)
        # Normalize (0: absence, 1: presence)
        leveled_df = normalize_dataframe(sub_df, criteria=1)    
        leveled_df["Total"] = leveled_df.sum()

        # normalize it (I had to do it with a for loop, sorry viewer)
        row_number, col_number = leveled_df.shape
        
        # last row, row_number-1, is the "All" row
        for column in range(0, col_number):
            # first, get the relative abundance of each taxon on each group 
            leveled_df.iloc[row_number-1, column] = leveled_df.iloc[row_number-1, column]*100/(row_number-1) 
            

            # file will contain row: category name, columns: abundance of this taxa in the group
            filename = f"{category}_prevalence_plots/{category}_lvl{level}.tsv"

            to_write_df = leveled_df.transpose()
            to_write_df.to_csv(filename, sep="\t")
            
            for column in range(0, col_number):
                # if present in all, change it to 1 (global presence)
                # else change it to 0
                if leveled_df.iloc[row_number-1, column] < 100:                    
                    leveled_df.iloc[row_number-1, column] = 0
                else:
                    leveled_df.iloc[row_number-1, column] = 1
