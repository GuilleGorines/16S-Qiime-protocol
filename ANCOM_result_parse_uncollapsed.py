'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 1
CREATED: 02-02-2022
REVISED: 06-02-2022
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

parser=argparse.ArgumentParser(description="")
parser.add_argument("--taxonomy", required=True, dest="taxfile", help="Taxonomy file")
parser.add_argument("--data", required=True, dest="datafile", help="")
parser.add_argument("--ancom", required = True, dest="ancomfile", help="")
parser.add_argument("--percent-abundances", required = True, dest="abundances", help="")
args = parser.parse_args()

def create_taxonomy_dict(taxfile):
    # Create a dictionary with the hash as key, the identification as value

    taxdict = dict()

    with open(taxfile) as taxfile:
        taxfile = taxfile.readlines()
        taxfile = [item.replace("\n","").split("\t") for item in taxfile]


    for code, identification, confidence in taxfile[1:]:        
        taxdict[code] = identification

    return taxdict

def join_files(abundances_file, datafile, ancomfile):
    # join abundances.tsv, data.tsv & ancom.tsv files into a single list

    # open abundances file
    with open(abundances_file) as abundances_file:
        abundances_file = abundances_file.readlines()
        abundances_file = [line.replace("\n","").split("\t") for line in abundances_file]
    
    # headers first
    joined_list = [abundances_file[0] + ["-","-","-"], abundances_file[1] + ["clr","W","Reject null Hypothesis"]]
    
    # dict with key = identificator; value = whole line
    datadict = dict()

    # avoid first 2 lines (headers)
    for line in abundances_file[2:]:
        datadict[line[0]] = line
                
    # open datafile
    with open(datafile) as datafile:
        datafile = datafile.readlines()
        datafile = [line.replace("\n","").split("\t") for line in datafile]

    # avoid first line (header) and append to the dict the extra lines in the file
        for line in datafile[1:]:
            datadict[line[0]] = datadict[line[0]] + line[1:]
            
    # open ancomfile
    with open(ancomfile) as ancomfile:
        ancomfile = ancomfile.readlines()
        ancomfile = [line.replace("\n","").split("\t") for line in ancomfile]

    # avoid header and append to the dict the extra 1 line in the file
        for line in ancomfile[1:]:        
            datadict[line[0]] = datadict[line[0]] + [line[2]]
    
    # add the dict values (whole line w/ extra info) to the headers
    joined_list = joined_list + list(datadict.values())

    return joined_list

def identify_taxa(in_list, taxdict):
    # changes hashes with their identified taxa
    
    # generate the extra headers ("-" for empty slot to correspond each column with its value)
    changed_list = [ ["-"] + in_list[0], ["code"] + in_list[1]]
    
    # for item in the list, add the identification in the second column
    for position, _ in enumerate(in_list[2:]):
        identified_line = [in_list[position+2][0]] + [taxdict[in_list[position+2][0]]] + in_list[position+2][1:]
        changed_list.append(identified_line)
                        
    return changed_list

def tsv_from_list(inserted_list, filename):
    # Generate a TSV file from a list 
    # Eeach item in the list will be a line

    with open(filename,"w") as outfile:
        for line in inserted_list:
            for column in line:
                outfile.write(column)
                
                if line.index(column) != len(line) -1:
                    outfile.write("\t")
            outfile.write("\n")
    return 


# taxonomy dict from taxonomy file
taxonomy_dict = create_taxonomy_dict(args.taxfile)

# generate filename
joined_list_filename = str(args.datafile).split("/")[0].replace("ancom_","").replace("_dir","") + "_ancom_result_identified.tsv"

# join the files
joined_list = join_files(args.abundances, args.datafile, args.ancomfile)

# identify the codes
join_list_identified = identify_taxa(joined_list, taxonomy_dict)

# save the list to tsv
tsv_from_list(join_list_identified,joined_list_filename)