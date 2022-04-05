from qiime2 import Artifact
import pandas as pd

import sys
import os

def save_long_wide(df, filename):
    """
    Generates tsv for a dataframe and for its transposed
    """
    df.to_csv(f"{filename}_wide.tsv", sep="\t")
    df.transpose().to_csv(f"{filename}_long.tsv", sep="\t")
    return

def relative_abundances(df):
    """
    Obtain the relative abundance of the otus
    """
    df["Total"] = df.sum(axis=1)
    rownum, colnum = df.shape
    for row in range(rownum):
        for col in range(colnum-1):
            df.iloc[row, col] = df.iloc[row, col] * 100 / df.iloc[row, colnum-1]

    df.drop("Total", axis=1, inplace=True)   
    
    return df
    
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

def prevalences(df, metadata):
    """
    Calculate the prevalence for each group
    """
    
    df_prev = pd.concat([df, metadata], axis=1)
    category_dict, category_names_list = create_category_dict(metadata)
    
    for category, values in category_dict.items():
    
        prevalence_per_value = []

        for value in values:

            # Drop metadata columns
            sub_df = df_prev[df_prev[category] == value].drop(category_names_list, axis=1)
            
            # Normalize (0: absence, 1: presence)
            norm_df = normalize_dataframe(sub_df, criteria=1)
            norm_df.loc["Prevalence"] = norm_df.sum(axis=0)

            row_number, col_number = norm_df.shape
        
            for column in range(0, col_number):
                # Get the relative abundance of each taxon on each group
                norm_df.iloc[row_number-1, column] = norm_df.iloc[row_number-1, column]*100/(row_number-1)
                norm_df.rename({"Prevalence":value}, axis=0, inplace=True)

            prevalence_per_value.append(norm_df.loc[value].to_frame().transpose())        
        
        prevalence_df = pd.concat(prevalence_per_value)
        
        save_long_wide(prevalence_df, f"prevalence")  
    

def clean_dataframe(df):
    """
    Remove the columns ending with ;__
    """
    df = df.loc[:,~df.columns.str.endswith(";__")]

    return df


def artifact_from_df(df_in, filename):
    
    clean_qza = Artifact.import_data("FeatureTable[Frequency]", df_in)
    clean_qza.save(f"{filename}.qza")
    
    return

# Input 1: feature table in qza format
# Input 2: metadata in tsv format
# Input 3: name of the output directory
# Input 4: level (for naming purposes)

qza_in = sys.argv[1]
metadata_file = sys.argv[2]
outdir = sys.argv[3]
level = sys.argv[4]

# df from qza
df = Artifact.load(qza_in).view(pd.DataFrame)

# Save the absolute numbers
save_long_wide(df, f"{outdir}/raw/absolute_numbers_lvl_{level}")

# Save the relative numbers
rel_df = relative_abundances(df)
save_long_wide(rel_df, f"{outdir}/raw/relative_numbers_lvl_{level}")

# Read metadata
metadata = pd.read_csv(
    metadata_file,
    sep='\t',
    header=0,
    index_col=0
    )

prevalences(df, metadata)

# Clean metadata
clean_df = clean_dataframe(df)

# save into an artifact
artifact_name = qza_in.replace("raw", "clean")

artifact_from_df(clean_df, artifact_name)
save_long_wide(clean_df, "_clean")

rel_clean_df = relative_abundances(clean_df)
save_long_wide(rel_clean_df, "_clean")

prevalences(clean_df, metadata)

