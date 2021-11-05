import sys
import os

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
    
     with open(f"{truepath}/samplesheet.tsv","w") as outfile:
            
        outfile.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    
        for file1,file2,samplename in final_groups:
            
            final_samplename = samplename.replace("_R","")

            targetdir = f"{truepath}/{final_samplename}"

            os.mkdir(targetdir)

            os.replace(f"{truepath}/{file1}",f"{targetdir}/{file1}")
            os.replace(f"{truepath}/{file2}",f"{targetdir}/{file2}")
            
            outfile.write(f"{final_samplename}\t{targetdir}/{file1}\t{targetdir}/{file2}\n")
       
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
