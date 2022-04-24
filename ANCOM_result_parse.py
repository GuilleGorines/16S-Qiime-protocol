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
    df_percent_abundances.rename({"Group" : argument},axis=0)
    
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

df_ancom, df_data, df_percent_abundances = export_qzv(args.qzv_in, args.metadata_column)

# generate first output
df_out_1 = pd.concat([df_ancom, df_data, df_percent_abundances], axis=1)
save_long_wide(df_out_1, "_")

# get the significative data 
significative_taxa = get_significative_taxa(df_out_1)

# import relative abundances
rel_abs_df = pd.read_csv(rel_freq_in, header=0, index_col=0, delimiter="\t")

# Import metadata
column_df = pd.read_csv(args.metadata, header=0, index_col=0, delimiter="\t").pop(args.metadata_column)
df_out_2 = pd.concat([df_ancom, df_data, pd.concat([pd.DataFrame(column_df).transpose(), rel_abs_df], axis=0)], axis=1)

# if there are any significative taxa
# generate the heatmap with dendrogram plot
if significative_taxa is not None:
    sig_tax_abundances = rel_abs_df.loc[significative_taxa, :]
    
    # change the headers of the table
    # get current names
    rownames = sig_tax_abundances.index
    
    # change current name into new future name
    # Get the values from the wanted columns
    newnames = [digest_name(item) for item in rownames]
    newnames = [f"{item[-1]}; {item[-2]}" if "uncultured" in item[-1] and len(item) > 2 else item[-1] for item in newnames]
    
    namedict = { row : newname for row, newname in zip(rownames, newnames)}
    
    figure_df = sig_tax_abundances.rename(index=namedict)
    
    # associate color code to metadata
    color_codes = dict(zip(column_df.unique(), ["green", "red", "blue", "purple", "grey"]))
    col_colors = column_df.map(color_codes)
    
    figure = sns.clustermap(figure_df,
                  col_colors=col_colors,
                  row_cluster=False,
                  dendrogram_ratio=(0, .15),
                  cbar_pos=(0.9, 0.1, .05, .25),
                  cmap="Greens",
                  figsize=(15,10),
                  )

    figure.savefig("hmap_xsamples_ytaxa.png")
    
    reverse_figure_df = figure_df.transpose()
    
    figure = sns.clustermap(reverse_figure_df,
                            row_colors=col_colors,
                            col_cluster=False,
                            dendrogram_ratio=(0.15, 0),
                            cbar_pos=(0.9, 0.1, .05, .20),
                            cmap="Greens",
                            figsize=(10,15),
                           )
    
    figure.savefig("hmap_xtaxa_ysamples.png")
        
    


