import argparse
import pandas as pd
from qiime2 import Artifact

parser = argparse.ArgumentParser(description='Generate the abundances and prevalences from the qza feature table')

parser.add_argument('--qza-in', help='feature table in qza format', dest="qza_in", required=True)
parser.add_argument('--mode', help='mode to choose from', choices= ["full", "filt"], dest="mode", required=True)
args = parser.parse_args()

# Open feature table
feature_table_df = Artifact.load(args.qza_in).view(pd.DataFrame)

# Sum the number of features (this is, the depth)
feature_table_df["Total depth"] = feature_table_df.sum(axis=1).astype(int)
depths_df = pd.DataFrame(feature_table_df.loc[:,"Total depth"])

# Output the depth to 
depths_df.sort_values("Total depth").transpose().to_csv(f"depths_{args.mode}.tsv", sep="\t")