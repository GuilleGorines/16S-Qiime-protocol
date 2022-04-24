import argparse
import uuid
import os
import shutil
import pandas as pd
from qiime2 import Visualization
import seaborn as sns; sns.set_theme(color_codes=True)

def save_long_wide(df, filename):
    """
    Generates tsv for a dataframe and for its transposed
    """
    df.to_csv(f"{filename}_wide.tsv", sep="\t")
    df.transpose().to_csv(f"{filename}_long.tsv", sep="\t")
    return

def export_qzv(qzv_in):
    # generate the tmp dir
    # export qzv to a tmp dir
    # get the needed dataframes
    # remove the tmp dir

    tmpdir = uuid.uuid4()

    while tmpdir in os.listdir():
        tmpdir = uuid.uuid4()


    qzv = Visualization.load(qzv_in)
    qzv.export_data(tmpdir)
    
    # ancom table
    df_ancom = pd.read_csv(f"{tmpdir}/ancom.tsv", sep="\t", index_col=0)
    
    # Data table
    df_data = pd.read_csv(f"{tmpdir}/data.tsv", sep="\t", index_col=0)
    
    # Add extra row to avoid NAs
    df_data.loc["Group"] = 2 * ["-"]
    
    # remove the "w"
    df_ancom.drop(["W"], axis=1, inplace=True)

    # Percent abundances
    df_percent_abundances = pd.read_csv(f"{tmpdir}/percent-abundances.tsv", sep = "\t", index_col=0)

    shutil.rmtree(tmpdir)

    return df_ancom, df_data, df_percent_abundances

def get_significative_taxa(df):
    # Get differentially expressed taxa
    # Those with "Reject null hypothesis" set as True
    significative_taxa = df[df["Reject null hypothesis"] == True].index

    if len(significative_taxa) == 0:
        print("No significative data found.")
        return None
    else:
        return list(significative_taxa)

parser = argparse.ArgumentParser(description='Generate the abundances and prevalences from the qza feature table')
parser.add_argument('--qzv-in'         , help='feature table in qza format', dest="qzv_in", required=True)
parser.add_argument('--metadata-file'  , help='metadata in tsv format', dest="metadata", required=True)
parser.add_argument('--metadata-column', help='column to organize the metadata by', dest="metadata_column", required=True)
parser.add_argument('--outdir'         , help='name of the output directory', dest="outdir", required=True)
parser.add_argument('--mode'           , help='mode to choose from (full/filt)', choices= ["full", "filt"], dest="mode", required=True)
parser.add_argument('--state'          , help='state to choose from', choices= ["raw", "clean"], dest="state", required=True)
parser.add_argument('--level'          , help='level (for naming purposes)', dest="level", default=None)
parser.add_argument('--rel-freq-file'  , help='relative frequences file', dest="relfreq_in", required=True)
args = parser.parse_args()

# generate the path where the relative abundances will be
rel_abundances_path = f"../../09-qiime2_collapse_numbers/{args.mode}/{args.state}/lvl_{args.level}/{args.state}/relative_numbers_lvl_{args.level}_{args.state}_long.tsv"

df_ancom, df_data, df_percent_abundances = export_qzv(args.qzv_in)

df_out_1 = pd.concat([df_ancom, df_data, df_percent_abundances], axis=1)

# get the significative data 
significative_taxa = get_significative_taxa(df_out_1)

if significative_taxa is not None:
    rel_abs_df = pd.read_csv(args.relfreq_in, header=0, index_col=0, delimiter="\t")
    sig_tax_abundances = rel_abs_df.loc[:, significative_taxa]

    metadata_df = pd.read_csv(args.metadata, header=0, index_col=0, delimiter="\t").loc[args.metadata_column]

    color_codes = dict(zip(metadata_df.pop(args.metadata_column).unique()), ["green", "red", "blue", "purple", "grey"])
    


