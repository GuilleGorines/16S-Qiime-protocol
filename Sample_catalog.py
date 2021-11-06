import sys
import os
import re

# Datadir is typically the RAW 
datadir = sys.argv[1]


def find_longest_match(list_string1,list_string2):
    
    punctuation_dict = {}
    
    for string1 in list_string1:
        
        punctuation_dict[string1] = []
        dissected_1 = [letter for letter in string1]
        
        for string2 in list_string2:
            
            dissected_2 = [letter for letter in string2]
            samplename = ""
            
            for position in range(0, min(len(dissected_1),len(dissected_2))):
            
                if dissected_1[position] == dissected_2[position]:
                    samplename += dissected_1[position]
                else:
                    punctuation_dict[string1].append([string2,samplename])
                    break
    
    return punctuation_dict

def find_best_match(punctuation_dict):
    
    final_groups = []
    
    for key,value in punctuation_dict.items():
        
        secondstring = ""
        samplename = ""
        punctuation = 0
        
        for second,shared in value:
            if len(shared) > punctuation:
                secondstring = second
                samplename = shared
                punctuation = len(shared)
                
            else:
                pass
            
        final_groups.append([key,secondstring,samplename])
        
    return final_groups

def create(final_groups, truepath):

    samplesheet = truepath + "/samplesheet.tsv"
    
    with open(samplesheet,"w") as outfile:
            
       outfile.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    
        for file1,file2,samplename in final_groups:
            
            final_samplename = samplename.replace("_L001_R","")
            final_samplename = re.sub(r"_S\d{2,3}$","", final_samplename)

            targetdir = f"{truepath}/{final_samplename}"

            os.mkdir(targetdir)

            file_1_to_be_replaced = truepath + "/" + file1
            file_1_replacement = targetdir + "/" + file1

            file_2_to_be_replaced = truepath + "/" + file2
            file_2_replacement = targetdir + "/" + file2

            os.replace(file_1_to_be_replaced, file_1_replacement)
            os.replace(file_2_to_be_replaced, file_2_replacement)
            
            line = final_samplename + "\t" + file_1_replacement + "\t" + file_2_replacement + "\n"
            outfile.write(line)
       
R1_files = []
R2_files = []

truepath = os.path.realpath(datadir)

for item in os.listdir(datadir):
    if "R1" in item:
        R1_files.append(item)
    elif "R2" in item:
        R2_files.append(item)

punctuation_dict = find_longest_match(R1_files,R2_files) 
final_groups = find_best_match(punctuation_dict)
create(final_groups,truepath)
