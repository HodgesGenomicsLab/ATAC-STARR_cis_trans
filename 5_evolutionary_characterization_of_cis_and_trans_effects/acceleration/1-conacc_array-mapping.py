#!/usr/bin/env python
# coding: utf-8


import argparse
import glob
import os
import sys
import subprocess


sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")

import chr_functions
import split_filename

    

# %% ARGPARSE arguments.

arg_parser = argparse.ArgumentParser(description="run phylop estimates of human sequence conservation/ acceleration")

arg_parser.add_argument("-i", "--index_number", help='array number', type=int)
arg_parser.add_argument("-m", "--mapping", help='mapping file w/ branches, model and input file')

# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

IDX = args.index_number  # the index
MAPPING_FILE = args.mapping  # the mapping file to read from. 


# FUNCTIONS


def read_mappingfile_index(i, mapping_file):
    
     with open(mapping_file, "r") as reader:
        lines = reader.readlines()
        
        for n, line in enumerate(lines):

            if n == i:
                index, branch, model, file, msaway, outfile = line.split("\t")
                outfile = outfile.strip("\n")
                chrnum, n_ = (file.split("/")[-1]).split("-")

                phylop_path = "/".join(file.split("/")[:-1])
                outpath = "/".join(outfile.split("/")[:-1])


                return index, branch, model, chrnum,n_, phylop_path, msaway, file, outfile, outpath
            
            
def get_model(model):
    
    model_dict = {
        "full":"/dors/capra_lab/data/ucsc/hg38/multiz30way/hg38.phastCons30way.mod",
        "hg38-rhemac8":"/dors/capra_lab/data/ucsc/hg38/multiz30way/hg38.phastCons30way_hg38-rheMac8.mod",
        "rhemac8_noowm":"/dors/capra_lab/data/ucsc/hg38/multiz30way/hg38.phastCons30way_rheMac8_noOWM.mod",
        "hg38_noapes":"/dors/capra_lab/data/ucsc/hg38/multiz30way/hg38.phastCons30way_hg38_noAPES.mod",
        }
    
    return model_dict[model]


def get_maf(msaway, chrnum):
    
    maf_zipped = os.path.join("/dors/capra_lab/data/ucsc/hg38", f"multiz{msaway}", "maf", f"{chrnum}.maf.gz")
    maf_unzipped = maf_zipped.split(".gz")[0]
    
    return maf_zipped, maf_unzipped


def get_phylop_bin():
    
    return "/dors/capra_lab/bin/./phyloP"



def run_phylop(msaway, ocr, n, chrnum, path, random_seed, branch, model):

    print(ocr, chrnum)

    #msaway = str(msa) + "way"

    # neutral tree
    mod = get_model(model)  # get the dictionary of the models

    # multiple sequence alignment file
    maf_zipped, maf_unzipped = get_maf(msaway, chrnum)

    # maf needs to be unzipped?

    if os.path.exists(maf_unzipped) is False:
        cmd = f"gunzip {maf_zipped}"
        subprocess.call(cmd, shell=True)

    # make outpath
    outpath = os.path.join(
        path, f"multiz{msaway}_br-{branch}_mod-{model}")

    try:
        os.mkdir(outpath)
    except FileExistsError:
        pass

    # make outfile
    outf = os.path.join(outpath, f"{chrnum}_{n}_conacc.bed")

    # Already done phyloP analysis on this file?
    if os.path.exists(outf) is False or os.path.getsize(outf) == 0:
        phylop = get_phylop_bin()
        
        # run phyloP!
        cmd = f"{phylop} --features {ocr} --msa-format MAF --method LRT --branch {branch} --mode CONACC -d {random_seed}         -g {mod} {maf_unzipped}> {outf}"

        # write run to log
        runlog_f = os.path.join(outpath, "runlog.txt")
        with open(runlog_f, "a") as runlog:
            runlog.write(cmd + "\n\n")

        # print(cmd)
        subprocess.call(cmd, shell=True)

        # check results
        if os.path.getsize(outf) > 0:

            # delete temp
            temp = os.path.join(path, f"temp_{chrnum}.bed")
            if os.path.exists(temp) is True:
                os.remove(temp)
                print("removed", temp)

        else:
            print("this didn't run", ocr)
    else:
        print("already processed", outf)

        # delete temp
        temp = os.path.join(path, f"temp_{chrnum}.bed")
        if os.path.exists(temp) is True:
            os.remove(temp)
            print("removed", temp)
        
    return outf


# # MAIN

def main(argv):
    index, branch, model, chrnum, n, phylop_path, msaway, file, outf, outpath = read_mappingfile_index(IDX, MAPPING_FILE)

    os.chdir(phylop_path)  # change directory
    
    concat = os.path.join(outpath, f"{chrnum}_conacc.bed")
        
    print(IDX, MAPPING_FILE, "\n\n running", index, "\n", branch, model, chrnum, msaway,  "\n", file, n, "\n\n")
    
    run_phylop(msaway, file, n, chrnum, phylop_path, RANDOM_SEED, branch, model)   
    

if __name__ == "__main__":
    main(sys.argv[1:])






