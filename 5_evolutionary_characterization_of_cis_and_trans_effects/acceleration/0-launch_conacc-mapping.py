#!/usr/bin/env python
# coding: utf-8



import argparse
import glob
import itertools 
import numpy as np
import os, sys
import subprocess
import time


sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import count_lines as cl

sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")
import config_readwrite as crw
import chr_functions
import split_filename


###
# ARGS
###
arg_parser = argparse.ArgumentParser(description=" describe argparse")


arg_parser.add_argument(
                        "-b", "--bedfile", 
                        help='bed file to run phyloP on.\
                        Can also be a redo mapping file.\
                        e.g. redo.txt'
                        )

arg_parser.add_argument(
                        "-br", "--branches", 
                        choices = ['hg38', 'rheMac8', 'hg38-rheMac8'], 
                        help='str, target branch to measure conservation/acceleration on',
                        default = "hg38",
                        type=str,
                        )

arg_parser.add_argument(
                        "-msa", "--multiz", 
                        choices =["20way", "30way", "100way"],
                        help='str, 20way, 30way, 100way multiz alignments in hg38',
                        default="30way",
                        type=str
                        )

arg_parser.add_argument(
                        "-mod", "--model", 
                        choices=['full', 'rheMac8_noOWM', 'hg38_noAPES'],
                        help='str, neutral model used to estimate conservation/acceleration', 
                        default="full",
                        type=str
                        )

arg_parser.add_argument(
                        "-o", "--outdir", 
                        help='str, path to output directory',
                        type=str
                        )

arg_parser.add_argument("-l", "--lines", 
                        help='int, n lines to split files into',
                        default=100, 
                        type=int
                       )

###
# PARSE THE ARGUMENTS
###

args = arg_parser.parse_args()

BEDF = args.bedfile  # the regions to test for acceleration 
BRANCHES = args.branches  # the branches to test.
MSAWAY = args.multiz  # multiple sequence alignment, "30way"
MODELS = args.model
DATA_PATH = args.outdir
LINES = args.lines  # nlines to split mapping file into. 
PATH, FILENAME, SAMPLE_ID = split_filename.split_filename(BEDF)

BRANCHES = list(BRANCHES) #['hg38', 'rheMac8', 'hg38-rheMac8'] # code is written to accommodate lists, but argparse handles single str, so make a list with one element

MODELS =list(MODELS)#['full', 'rheMac8_noOWM', 'hg38_noAPES']  # code is written to accommodate lists, but argparse handles single str, so make a list with one element


###
# FUNCTIONS
###

def split_by_chr(f):
    
    """
    return path containing bedfile into chromosomes .bed files
    
    1. get file name
    2. make chr_path name
    3. make chr_path dir
    4. split the bed file by chromosome, 
    5. go to the directory
    6. return chr_path, list of split files
    
    """
    
    #1
    path, file_name, sample = split_filename.split_filename(f)
    
    #2
    chr_path = os.path.join(path, f"chr-{sample}") # make a path
    
    #3
    if os.path.exists(chr_path) is False:
        os.mkdir(chr_path)  # make the chr dir
    
    # 4
    chr_functions.split_into_chr_bed(f, chr_path)  # split by chromosome
    
    #5
    os.chdir(chr_path)
    return chr_path, glob.glob("chr*.bed"), sample



# in case you need to split file on size before getting started
def split_by_line(f, data_path, chr_num, nlines):
    
    """
    return list of files split up into n equal lines
    
    1. make a chr-specific directory
    2. check that you haven't already split
    3. if not, split the files by n lines
    4. move files to chr-specific directory
    5. return a list of the splits. 
    
    """
    
    #1
    chr_path = os.path.join(data_path, chr_num) # make dir for chromosome splits
    
    try:
        os.mkdir(chr_path)
    except FileExistsError:
        pass
    
    # change dir to the output chr path (not the original CHR_PATH, 
    # where file is split on CHR, but not line number)
    
    #2
    split_fs = glob.glob(os.path.join(chr_path, f"{chr_num}-*"))

    #3 split the file in command line into sizes of nlines
    os.chdir(data_path)
    
    cmd = f"split -l {nlines} {f} {chr_num}-"
    
    if len(split_fs) ==0:
        print("splitting")
        subprocess.call(cmd, shell = True)
        
        #4
        mv_cmd = f"mv {chr_num}-* {chr_path}"
        subprocess.call(mv_cmd, shell = True)
        
        split_fs = glob.glob(f"{chr_path}/{chr_num}-*")

    else:
        print("already split")
    
    #5
    return split_fs



def make_run_list(branches, models, chrs):
    
    """
    return list of branch, model, and chr pairs to write in mapping file, run
    
    exclude specific branch, model tuples that don't make sense to test. 
    
    inputs 
        branches (list)
        models (list)
        chr_list (list)
        
    output 
        runs (list of lists)
    """
    runs = []
    no_runs = [('hg38', 'rheMac8_noOWM') , 
               ('hg38-rheMac8', 'rheMac8_noOWM'), 
               ('hg38-rheMac8', 'hg38_noAPES'),
              ('rheMac8', 'hg38_noAPES') 
              ] # don't run these tuples. Not interested yet in these results
    for b in branches:
        for m in models:
            for c in chrs:           
                combo = [b, m, c]
                if combo not in runs and (b,m) not in no_runs:
                    runs.append(combo)
                    
    return runs



def write_mapping_file(mapping_file, run_list, n, msaway):
    
    """
    write to mapping file for slurm array. Line includes index, branch, model, and file to run
    
    0. make the outfile that the results will be written to
    1. make the row as a tab-delimited string
    2. open the mapping file 
    3. append the row to the file
    4. increase the n counter by the number of rows,
    5. return n to count for next set of rows. 
    
    inputs: mapping file (.txt), run_list (list), n (int)
    outputs: n (int)
    """
    
    for r in run_list:
        
        #0 make the outfile that results will be written to
        branch, model, chr_file = r  # make variables from list
        chr_, chr_n = (chr_file.split("/")[-1]).split("-")  # chr3-aj
        
        phylop_path = "/".join(chr_file.split("/")[:-1])  # path to chr_file
        
        outpath = os.path.join(
        phylop_path, f"multiz{msaway}_br-{branch}_mod-{model}") # path to write to
        
        outf = os.path.join(outpath, f"{chr_}_{chr_n}_conacc.bed")  # file to write during phylop run 


        #1 make the row as tabbed-string (index, branch, model, and file to run)
        row =  str(n) + "\t"+ "\t".join(r) +"\t"+ msaway + "\t"+ outf + "\n"

        #2 open the mapping file
        with open(mapping_file, "a") as writer:

            #3 write the row 
            writer.write(row)
            writer.close()
        #4
        n+=1 # for counting the next row. 
        
    #5
    return n



def run_conacc_slurm_array(mapping_file, num_files):
    
    """
    return slurm script command with array argument
    
    inputs 
        mapping file - branch, model, run file, and msaway info in tab-seperated format (txt)
        num_files - length of mapping file/ number of files to run (int)

    hardcoded
        slurm script path
        memory for array run
        limit to n files run in array (%10)
        
    output
        slurm command (str) with mapping file, array format, and memory specified
    """

    script = os.path.join("/data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/phylop", "conacc_array-mapping.slurm")


    mem = "--mem=36GB"
    array = f"--array [0-{num_files}]%30"  # only run ten files at a time

    # make the command

    cmd = f"sbatch {array} {mem} {script} {mapping_file}"  #{chrnum} {branches} {msaway} {mod}"
    
    print(cmd)

    # return it
    return cmd



def main(argv):
    """
    launch a slurm array with .bed file split on chromosome number and equal lines
    
    inputs - 
        bedfile (.bed) 
        data_path to save the split up files (dir)
        n lines - number of lines to write per file
        BRANCHES - list of branches to test
        MODELS - list of neutral models to test acceleration against. 
        MSAWAY - multiz alignment to use. 
        
    outputs - 
        mapping_file with all parameters for running phyloP (branch, model, multiz)
        launched slurm command
    """
    
    if ".txt" not in BEDF:

        #1 split bedfile by chr
        chr_path, chrs, sample_id = split_by_chr(BEDF)

        #2 go to the chr dir
        os.chdir(chr_path)

        #3 exclude these chromosomes by taking set difference.
        excl_chr = set(['chrX.bed', 'chrY.bed', 'chrM.bed', 
                        'chr14_KI270726v1_random.bed', 'chr16_KI270728v1_random.bed'
                        'chr14_KI270722v1_random.bed', 'chr16_KI270728v1_random.bed', 
                        'chr14_KI270722v1_random.bed'
                       ])

        chrs_ = list(set(chrs).difference(excl_chr))


        mapping_file = os.path.join(chr_path, f"mapping_{sample_id}.txt")


        # iterate through chromosomes
        for chr_ in chrs_: 

            chr_num = chr_.split(".bed")[0]

            # split chr_files into 1000 line files.
            split_fs = split_by_line(chr_, chr_path, chr_num, LINES)

            run_list = make_run_list(BRANCHES, MODELS, split_fs)
            #print(run_list)

            # write to mapping file. 
            n = write_mapping_file(mapping_file, run_list, n, MSAWAY)

            print(len(run_list),"phylop chr x branch x model comparisons to run")

            # launch the slurm job
            cmd = run_conacc_slurm_array(mapping_file, n)


    else:
        mapping_file = BEDF
        n = cl.count_lines(mapping_file)
        cmd = run_conacc_slurm_array(mapping_file, n)
    
    
    subprocess.call(cmd, shell = True)

if __name__ == "__main__":
    main(sys.argv[1:])
