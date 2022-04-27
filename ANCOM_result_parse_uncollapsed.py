'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 2
CREATED: 02-02-2022
REVISED: 27-04-2022
DESCRIPTION: 
    Ad-hoc script to identify the differentially abundant taxa 
INPUT:
    -taxonomy: qiime2 TSV relating taxonomy and hash
    -data: ANCOM data.tsv file
    -ancom: ANCOM ancom.tsv file
    -percent-abundances: ANCOM percent-abundances.tsv file
OUTPUT:
    TSV containing Unified data from ANCOM with the identified hashes
'''
import argparse
import uuid
import os
import sys
import shutil
import pandas as pd
from qiime2 import Visualization
from qiime2 import Artifact
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)

def save_long_wide(df, filename, rowname, colname):
    """
    Generates tsv for a dataframe and for its transposed
    """
    df.to_csv(f"{filename}_row{rowname}_col{colname}.tsv", sep="\t")
    df.transpose().to_csv(f"{filename}_row{colname}_col{rowname}.tsv", sep="\t")
    return

def export_qzv(qzv_in, argument):
    # generate the tmp dir
    # export qzv to a tmp dir
    # get the needed dataframes
    # remove the tmp dir

    tmpdir = str(uuid.uuid4())

    while tmpdir in os.listdir():
        tmpdir = str(uuid.uuid4())

    qzv = Visualization.load(qzv_in)
    qzv.export_data(tmpdir)

    # ANCOM table
    df_ancom = pd.read_csv(f"{tmpdir}/ancom.tsv", sep="\t", index_col=0)
    
    # Data table
    df_data = pd.read_csv(f"{tmpdir}/data.tsv", sep="\t", index_col=0)
    newrow = pd.DataFrame.from_dict({argument : ["-","-"]}, orient="index", columns=["W","clr"])
    
    # Generate extra row to avoid NAs    
    df_data = pd.concat([newrow, df_data], axis=0)
    
    # remove the "w"
    df_ancom.drop(["W"], axis=1, inplace=True)
    
    # Generate extra row to avoid NAs
    newrow = pd.DataFrame.from_dict({argument : ["-"]}, orient="index", columns=["Reject null hypothesis"])
    df_ancom = pd.concat([newrow, df_ancom], axis=0)
    
    # Percent abundances
    df_percent_abundances = pd.read_csv(f"{tmpdir}/percent-abundances.tsv", sep = "\t", index_col=0)
    df_percent_abundances = df_percent_abundances.rename(index={"Group" : argument})

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

def digest_name(string):
    name_as_list = string.replace("d__","domain: " ).replace("p__", "phylum: ").replace("c__", "class: " ).replace("o__", "order: " ).replace("f__","family: " ).replace("g__","genus: "  ).replace("s__","species: ").split(";")
    return name_as_list


parser=argparse.ArgumentParser(description="")
parser.add_argument("--qza-taxonomy"   , help="Taxonomy file in qza format", dest="tax_qza", required=True, )
parser.add_argument("--qzv-ancom"      , help="ANCOM results file in qzv format", dest="ancom_qzv", required=True )
parser.add_argument("--mode"           , help='mode to choose from (full/filt)',  dest="mode", required=True, choices= ["full", "filt"],)
parser.add_argument("--metadata"       , help='metadata file', dest="metadata", required=True)
parser.add_argument("--metdata-column" , help="metadata column for the analysis", dest="metadata_column", required=True)
parser.add_argument("--rel-freq-file"  , help='relative frequences file', dest="relfreq_in", required=True)
args = parser.parse_args()

# Import ANCOM result
df_ancom, df_data, df_percent_abundances = export_qzv(args.ancom_qzv, args.metadata_column)
df_out_1 = pd.concat([df_ancom, df_data, df_percent_abundances], axis=1)

# Import taxonomy
tax_df = Artifact.load(args.tax_qza).view(pd.DataFrame)

# Generate new row so sample-origin goes first
# rowname = data
newrow = pd.DataFrame(["-"]*2).transpose()
newrow.columns = ["Consensus", "Taxon"]
newrow.index=[args.metadata_column]

# Add the new row
tax_df = tax_df.loc[list(df_out_1.index[1:])]
tax_df_meta = pd.concat([newrow, tax_df], axis=0)

# First result: complete ancom results
df_out_1 = pd.concat([tax_df, df_out_1], axis=1)
save_long_wide(df_out_1, f"1-Complete_ancom_result_{args.metadata_column}_uncollapsed_raw_{args.mode}","taxa","ancom-full")

# SECOND RESULT
# Import relative abundances
rel_freq_df = pd.read_csv(args.relfreq_in, header=0, index_col=0, delimiter="\t")

# Import metadata
column_df = pd.read_csv(args.metadata, header=0, index_col=0, delimiter="\t").pop(args.metadata_column)

# Generate second dataframe
# (ANCOM results with relative frequency for each sample)
df_out_2 = pd.concat([tax_df, df_ancom, df_data, pd.concat([column_df, rel_freq_df], axis=0)], axis=1)
new_index = [args.metadata_column] + list(df_out_2["Taxon"])[1:]
ids = ["-"] + list(df_out_2.index)[1:]
df_out_2.index = new_index
df_out_2["Taxon"] = ids
df_out_2 = df_out_2.rename(columns={"Taxon" : "ID"})
save_long_wide(df_out_2, f"2-Ancom_result_w_rel_freq_{args.metadata_column}_uncollapsed_raw_{args.mode}", "tax-id", "tax-ancom-relfreq")

# THIRD RESULT
# (Only significative results)
significative_taxa = get_significative_taxa(df_out_2)

if significative_taxa is not None:
    
    rel_freq_df = pd.concat([tax_df, rel_freq_df], axis=1)
    rel_freq_df.index = list(rel_freq_df["Taxon"])
    rel_freq_df = rel_freq_df.drop(["Taxon","Consensus"], axis=1)
    
    sig_tax_abundances = rel_freq_df.loc[significative_taxa, :]
    
    # change the headers of the table
    # get current names
    rownames = sig_tax_abundances.index
    
    # change current name into new future name
    # Get the values from the wanted columns
    newnames = [digest_name(item) for item in rownames]
    newnames = [f"{item[-1]}; {item[-2]}" if ("uncultured" in item[-1] and len(item) > 2) or ("__" in item[-1] and len(item) > 2) else item[-1] for item in newnames]
    
    namedict = { row : newname for row, newname in zip(rownames, newnames)}
    
    figure_df = sig_tax_abundances.rename(index=namedict)
    
    # associate color code to metadata
    color_codes = dict(zip(column_df.squeeze().unique(), ["#00AA5A", "#C0AB52", "#e16a86",  "#00A6CA",  "#C699E7", "grey"]))
    col_colors = column_df.squeeze().map(color_codes)
    
    figure = sns.clustermap(figure_df,
                  col_colors=col_colors,
                  row_cluster=False,
                  dendrogram_ratio=(0, .15),
                  cbar_pos=(0.9, 0.1, .05, .25),
                  cmap="Greens",
                  figsize=(15,10),
                  )

    handles = [Patch(facecolor=color_codes[name]) for name in color_codes]
    plt.legend(handles, color_codes, title=args.metadata_column,
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure)

    figure.savefig(f"hmap_{args.metadata_column}_unleveled_raw_{args.mode}_xsamples_ytaxa.png")
    
    reverse_figure_df = figure_df.transpose()
    
    figure = sns.clustermap(reverse_figure_df,
                            row_colors=col_colors,
                            col_cluster=False,
                            dendrogram_ratio=(0.15, 0),
                            cbar_pos=(0.9, 0.1, .05, .20),
                            cmap="Greens",
                            figsize=(10,15),
                           )
    plt.legend(handles, color_codes, title=args.metadata_column,
        bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure)
    
    
    
    figure.savefig(f"hmap_{args.metadata_column}_uncollapsed_raw_{args.mode}_xtaxa_ysamples.png")

    # Third file
    # Only significative taxa involved
    relevant_rows = [args.metadata_column] + significative_taxa
    df_out_3 = df_out_2.loc[relevant_rows, :]

    # New df for Mode, level and state columns
    newcols = pd.DataFrame([
        ["-"] + [args.mode] * (df_out_3.shape[0]-1),
        ["-"] + ["uncollapsed"] * (df_out_3.shape[0]-1),
        ["-"] + ["raw"] * (df_out_3.shape[0]-1)
        ]).transpose()
    newcols.columns = ["Mode", "Level", "State"]
    newcols.index = df_out_3.index

    # Add the new columns
    df_out_3 = pd.concat([newcols, df_out_3], axis=1).drop(["Consensus", "ID"], axis=1)
    save_long_wide(df_out_3, f"3-Significative_results_{args.metadata_column}_uncollapsed_raw_{args.mode}", "taxa", "ancom-samples")


