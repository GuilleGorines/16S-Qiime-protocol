#!/usr/bin/env python
import sys

biom = sys.argv[1]
outfile_name = sys.argv[2]

# read the (readable) biom table, extract the header
# organism: row, number of sublists
# samples: columns, number of items in any sublist
with open(biom, "r") as infile:
    biom = infile.readlines()
    biom = [line.replace("\n","").split("\t") for line in biom[1:]]

initial_organisms = len(biom)-1

print("------------- BEGIN -------------")
print("Organism entries: ", initial_organisms)
print("Sample number: ", len(biom[0])-1, "\n")

print("------------- DEDUPLICATION BEGIN -------------")
# avoid identical taxids
# 1-save the header with the sample names
absolute_feature_table = [biom[0]]

absolute_feature_table_dict = {}

for item in biom[1:]:
    tax_name = item[0]
    numbers = [float(number) for number in item[1:]]

    if tax_name in absolute_feature_table_dict.keys():
        print("Found duplicate identification: ", tax_name)
        for index in range(0,len(numbers)-1):
            absolute_feature_table_dict[tax_name][index] += numbers[index]
    else:
        absolute_feature_table_dict[tax_name] = numbers

print("\n")
# organism: number of keys in absolute_feature_table_dict -1, corresponding to the header presentation
# samples: number of entries in any sublist of absolute_feature_table_dict values

values_dict = [ item for item in absolute_feature_table_dict.values()]

print("------------- DEDUPLICATION ENDED -------------")
print("Organism entries: ", len(absolute_feature_table_dict.keys()))
print("Sample number: ", len(values_dict[1]))
duplicate_numbers = initial_organisms - len(absolute_feature_table_dict.keys())
print(duplicate_numbers," duplicated organisms were found.\n")
print("------------- END -------------\n")

for key,values in absolute_feature_table_dict.items():

    file_row = [key]
    file_row.extend(values)
    absolute_feature_table.append(file_row)

# transpose it through dict
absolute_feature_table_dict = {}

for index in range(0,len(absolute_feature_table[0])):
    absolute_feature_table_dict[index]=[item[index] for item in absolute_feature_table]

absolute_feature_table = [item for item in absolute_feature_table_dict.values()]

with open(outfile_name,"w") as outfile:
    for item in absolute_feature_table_dict.values():       

        line_to_write = "\t".join([str(subitem) for subitem in item])

        outfile.write(line_to_write)
        outfile.write("\n")