import argparse

import pandas as pd
from qiime2 import Visualization

 def export_qzv(qzv_in, tmpdir):
    # export qzv to a tmp dir
    # get the needed dataframes
    # remove the tmp dir
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


parser = argparse.ArgumentParser(description='Generate the abundances and prevalences from the qza feature table')

parser.add_argument('--qzv-in', help='feature table in qza format', dest="qzv_in", required=True)
parser.add_argument('--metadata-file', help='metadata in tsv format', dest="metadata", required=True)
parser.add_argument('--outdir', help='name of the output directory', dest="outdir", required=True)
parser.add_argument('--mode', help='mode to choose from (full/filt)', choices= ["full", "filt"], dest="mode", required=True)
parser.add_argument('--state', help='state to choose from', choices= ["raw", "clean"], dest="state", required=True)
parser.add_argument('--level', help='level (for naming purposes)', dest="level", default=None)
parser.add_argument('--rel-freq-file', help='relative frequences file', dest="relfreq_in", required=True)

args = parser.parse_args()


rel_abundances_path = f"../../09-qiime2_collapse_numbers/{mode_directory}/{args.state}/lvl_{args.level}/{args.state}/relative_numbers_lvl_{args.level}_{args.state}_long.tsv"