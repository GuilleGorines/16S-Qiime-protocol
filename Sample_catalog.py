import sys
import os
import re

# Datadir is typically the RAW 
datadir = sys.argv[1]

# destiny dir is typically ANALYSIS/00-reads
destiny_dir = sys.argv[2]

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

def create(final_groups, destiny_dir):          
    
    for file1,file2,samplename in final_groups:
        
        final_samplename = samplename.replace("_L001_R","")
        final_samplename = re.sub(r"_S\d{2,3}$","", final_samplename)

        targetdir = os.path.realpath(destiny_dir) + "/" + final_samplename

        extension_list = [".fastq",".fastq.gz",".fq", ".fq.gz"]

        for extension in extension_list:
            if file1.endswith(extension) and file2.endswith(extension):
                file_extension = extension
                break

        file1_nicename = final_samplename + "_R1" + file_extension 
        file2_nicename = final_samplename + "_R2" + file_extension 

        os.mkdir(targetdir)

        file1_truepath = os.path.realpath(datadir + "/" + file1)
        file2_truepath = os.path.realpath(datadir + "/" + file2)

        file1_destinypath = targetdir + "/" + file1_nicename
        file2_destinypath = targetdir + "/" + file2_nicename

        os.symlink(file1_truepath, file1_destinypath)
        os.symlink(file2_truepath, file2_destinypath)

    
R1_files = []
R2_files = []

truepath = os.path.realpath(destiny_dir)


for item in os.listdir(datadir):
    if "R1" in item:
        R1_files.append(item)
    elif "R2" in item:
        R2_files.append(item)

punctuation_dict = find_longest_match(R1_files,R2_files) 
final_groups = find_best_match(punctuation_dict)
create(final_groups,truepath)
