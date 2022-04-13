import sys

from qiime2 import Artifact
import pandas as pd

qza_in = sys.argv[1]
out_dir = sys.argv[2]
lvl = sys.argv[3]

# Open and get the metadata
qza_artifact = Artifact.load(qza_in)
df = qza_artifact.view(pd.DataFrame)

# Generate the tsv with the absolute data
# (normal and transposed)
df.to_csv(
    f"{out_dir}/absolute_numbers_lvl_{lvl}_wide",
    sep="\t",
)

df.transpose().to_csv(
    f"{out_dir}/absolute_numbers_lvl_{lvl}_long",
    sep="\t",
)

# add a "Total" line
df["Total"] = df.sum()

# obtain relative abundances (%)
rownum, colnum = df.shape
for row in rownum-1:
    for col in colnum-1:
        df.iloc[row, col] = df.iloc[row, col] * 100 / df.iloc[row, col-1]

# remove the "Total" column
df.drop("Total", axis=1)

# Generate the tsv with the relative abundances
# (normal and transposed)
df.to_csv(
    f"{out_dir}/relative_numbers_lvl_{lvl}_wide",
    sep="\t",
)

df.transpose().to_csv(
    f"{out_dir}/relative_numbers_lvl_{lvl}_long",
    sep="\t",
)