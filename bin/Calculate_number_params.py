from qiime2 import Artifact
import pandas as pd

import sys
import os
import argparse

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
    df.loc[:,"Total"] = df.sum(axis=1)
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

def prevalences(df, metadata, outdir, mode, level=None):
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

        if level is None:
            save_long_wide(prevalence_df, f"{outdir}/Prevalence_{category}_{mode}")
        else:
            save_long_wide(prevalence_df, f"{outdir}/Prevalence_{category}_lvl_{level}_{mode}")  

def clean_dataframe(df):
    """
    Remove the columns ending with ;__
    """
    df = df.loc[:,~df.columns.str.endswith(";__")]

    return df


def artifact_from_df(df_in, filename):
    
    clean_qza = Artifact.import_data("FeatureTable[Frequency]", df_in)
    clean_qza.save(f"{filename}")
    
    return

# Input 5: mode (collapsed to clean, original not to)

parser = argparse.ArgumentParser(description='Generate the abundances and prevalences from the qza feature table')

parser.add_argument('--qza-in', help='feature table in qza format', dest="qza_in", required=True)
parser.add_argument('--metadata-file', help='metadata in tsv format', dest="metadata", required=True)
parser.add_argument('--outdir', help='name of the output directory', dest="outdir", required=True)
parser.add_argument('--mode', help='mode to choose from', choices= ["collapsed", "full", "filt"], dest="mode", required=True)
parser.add_argument('--level', help='level (for naming purposes)', dest="level", default=None)

args = parser.parse_args()

if args.mode == "full" or args.mode == "filt":
    # df from qza
    df = Artifact.load(args.qza_in).view(pd.DataFrame)
    
    # Save the absolute numbers
    save_long_wide(df, f"absolute_numbers_{args.mode}")

    # Save the relative numbers
    rel_df = relative_abundances(df)
    save_long_wide(rel_df, f"relative_numbers_{args.mode}")

    # Prevalences
    # Read metadata
    metadata = pd.read_csv(
        args.metadata,
        sep='\t',
        header=0,
        index_col=0
        )

    os.mkdir(f"{args.outdir}/Prevalences_{args.mode}")
    prevalences(df=df, metadata=metadata, outdir=f"{args.outdir}/Prevalences_{args.mode}", level=None, mode=args.mode)

    print(f"Stats on {args.mode} un-collapsed table ended successfully!")
    sys.exit(0)

else:

    # df from qza
    df = Artifact.load(qza_in).view(pd.DataFrame)

    # Save the absolute numbers
    save_long_wide(df, f"{outdir}/raw/absolute_numbers_lvl_{level}_raw")

    # Save the relative numbers
    rel_df = relative_abundances(df)
    save_long_wide(rel_df, f"{outdir}/raw/relative_numbers_lvl_{level}_raw")

    # Read metadata
    metadata = pd.read_csv(
        args.metadata,
        sep='\t',
        header=0,
        index_col=0
        )

    # Create directory for prevalences
    os.mkdir(f"{outdir}/raw/Prevalences/")
    prevalences(df=df, metadata=metadata, outdir=f"{outdir}/raw/Prevalences/", level=args.level, mode="raw")
   
    # Clean dataframe
    clean_df = clean_dataframe(df)

    # Save into an artifact
    artifact_name = qza_in.replace("raw", "clean")
    artifact_from_df(clean_df, artifact_name)

    # Save the asbsolute and relative numbers for the clean table
    save_long_wide(clean_df, f"{outdir}/clean/absolute_numbers_lvl_{level}_clean")

    rel_clean_df = relative_abundances(clean_df)
    save_long_wide(rel_clean_df, f"{outdir}/clean/relative_numbers_lvl_{level}_clean")

    # Prevalences for the clean table
    os.mkdir(f"{outdir}/clean/Prevalences")
    prevalences(df=clean_df, metadata=metadata, outdir=f"{outdir}/clean/Prevalences/", level=args.level, mode="clean")

    print(f"Stats on level {args.level} collapsed table ended successfully!")
    sys.exit(0)