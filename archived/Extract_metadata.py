'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 1
CREATED: 23-03-2022
REVISED: 23-02-2022
DESCRIPTION: 
INPUT:
OUTPUT:
'''
import shutil
import os
import sys

import pandas as pd
from qiime2 import Visualization

infile = sys.argv[1] # "collapsed_table_lvl6.qzv"
outdir = sys.argv[2] # "lvl6"
outfile = sys.argv[3] # 

table_name = f"{outdir}/{outfile}"

# Open visualization
qzv_file = Visualization.load(infile)

# Export the data from the visualization onto a dir
qzv_file.export_data(outdir)

# Open the metadata tsv, remove the 
# type of data (they are all numeric)
# generate the long-table (transposed)
df = pd.read_csv(f"{outdir}/metadata.tsv", sep="\t", header=0, index_col=0)
df = df.drop("#q2:types")
df.to_csv(f"{table_name}_wide.tsv", sep="\t")
df.transpose().to_csv(f"{table_name}_long.tsv", sep="\t")

# Delete unwanted dirs & files
# Hardcoded but its always the same so
dirs_to_del = ["css", "js", "q2templateassets"]

for folder in dirs_to_del:
    shutil.rmtree(f"{outdir}/{folder}")

files_to_del = ["index.html", "metadata.tsv"]
for file in files_to_del:
    os.remove(f"{outdir}/{file}")