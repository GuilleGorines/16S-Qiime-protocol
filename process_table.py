#!/usr/bin/env python
import sys

taxonomy = sys.argv[1]
biom = sys.argv[2]
outfile_name = sys.argv[3]

# read taxonomy and extract first and second lines, ID and taxonomy, respectively
with open(taxonomy, "r") as infile:
    taxonomy = infile.readlines()
    taxonomy = [line.replace("\n","").split("\t") for line in taxonomy[1:]]
    ids_list = [line[0] for line in taxonomy]

# read the (readable) biom table, extract the header
with open(biom, "r") as infile:
    biom = infile.readlines()
    biom = [line.replace("\n","").split("\t") for line in biom[1:]]

    header = biom[0]

# initialize the changed file
absolute_feature_table = [header]

# change the ids in the biom for the species in the taxonomy
for line in biom[1:]:
    print(line)
    identification = line[0]
    index_of_species = ids_list.index(identification)
    # generate new line with index as name

    taxonomic_name = taxonomy[index_of_species][1]

    new_line = [taxonomic_name] + line[1:]
    print(new_line)
    absolute_feature_table.append(new_line)

# to transpose it, we'll put it in a dict first
absolute_feature_table_dict = {}

for index in range(0,len(header)-1):
    absolute_feature_table_dict[index]=[item[index] for item in absolute_feature_table]

with open(outfile_name,"w") as outfile:
    for item in absolute_feature_table_dict.values():

        line_to_write = "\t".join(item)

        outfile.write(line_to_write)
        outfile.write("\n")