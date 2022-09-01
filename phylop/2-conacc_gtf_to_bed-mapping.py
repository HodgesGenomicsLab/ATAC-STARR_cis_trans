import argparse
import glob
import os
import pandas as pd
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")

import split_filename

# %% ARGPARSE arguments.

arg_parser = argparse.ArgumentParser(description="turn gtf output to .bed file format")

arg_parser.add_argument("-m", "--mapping", help='mapping file w/ branches, model and input file')

# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

MAPPING_FILE = args.mapping  # the mapping file to read from. 

def drop_gtf_cols(F):
    
    path, file, sample_id = split_filename.split_filename(F)
    
    os.chdir(path)  # go to the path
    
    newF = os.path.join(path, f"{sample_id}-clean.bed") # path to write to 
   
    cmd = f"cut -f 1,4,5,6 {F} > {newF}" # cut the files w/ info
   

    # check that you haven't dropped these columns yet
    if os.path.exists(newF) is False:
        subprocess.call(cmd, shell = True)  # run cut command in shell. 
    
    else:
        print("already dropped cols", F)
  
    return newF

    
def format_cols(F):
    
    f = pd.read_csv(F, sep = '\t', header = None)  # open file
    
    ncols = len(list(f))


    if ncols == 4:
        
        
        cols = ["chr", "start", "end", "phylop"]
        
        f.columns = cols
        
        f["region_id"] = f["chr"] +":"+ f["start"].map(str) +"-"+ f["end"].map(str)
        
        f.to_csv(F, sep = '\t', index = False)  # save results 

    
def main(argv):    
    
    """
    format cleaned phyloP results into .bed file
    """
    
    with open(MAPPING_FILE, "r") as mapping_file:
        
        first_line = mapping_file.readline()  
        
        """
        example - first_line
        0       hg38    full    /data/hodges_lab/ATAC-STARR_B-cells/results/results_human-evolution/regions/chr-all_uniq_diffAct_regions2/chr8/chr8-ad     30way   /data/hodges_lab/ATAC-STARR_B-cells/results/results_human-evolution/regions/chr-all_uniq_diffAct_regions2/chr8/multiz30way_br-hg38_mod-full/chr8_ad_conacc.bed
        """
        
        path = "/".join((first_line.split("\t")[3]).split("/")[:-2])
        msaway = first_line.split("\t")[4]
        
        FS = glob.glob(os.path.join(path, "chr*", f"multiz{msaway}_**/chr*_conacc.bed"))

        for F in FS:
                       
            newF = drop_gtf_cols(F)

            format_cols(newF)

if __name__ == "__main__":
    main(sys.argv[1:])