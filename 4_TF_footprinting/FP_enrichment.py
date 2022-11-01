
from functools import partial
from itertools import product
import glob
import numpy as np
import os
import pandas as pd
import statsmodels
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")

import config_readwrite as crw
import fet
import chr_functions
import split_filename

COL = sys.argv[1]
CL = sys.argv[2]


name = "/data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/config"
config, configfile_name = crw.read_config(name)

REGIONS = config["CIS_TRANS"]["regions"]
REGION_ANNOT = config["CIS_TRANS"]["regions_annotations"]


ID_TAG = config["TF_FOOTPRINTING_JASPAR"]["ID_TAG"]
FP_RE = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["FP"]
SHUF = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["FP_SHUF"]


FP_SHUF_OR_RE = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["FP_SHUF_OR"]  # write

RE = config["TF_FOOTPRINTING"]["results"]


###
# FUNCTIONS
###

def get_fp_counts(test, df):
    
    """
    Return sets of region_ids in activity category (e.g. cis-only region_ids)
    """
    
    test_region_ids = set(df.loc[df[test]==1, "region_id"])
    nottest_region_ids = set(df.loc[df[test]!=1, "region_id"])
    
    return test_region_ids, nottest_region_ids


def make_shuf_vector(shuf_file, tf, ytest):
    """
    Return sum of FP overlaps, non-overlap regions for a single TF across the shuffle files.
    
    1. per shuf_file, read only region_id column and tf column
    2. include regions that match the original ATAC-STARR regions in size. 
    3. count the number of TF FP overlaps
    4. sum across shuffled files
    5. return sum of FP overlaps, non-overlaps
    
    Note 
    - limit the overlap to the matching regions that shuffles were generated from. 
    
    """
    #1
    cols = ["region_id", tf]
    n_overlaps, n_nooverlaps = 0,0

    shuf = pd.read_csv(shuf_file, sep = '\t', low_memory = False, usecols = cols)
    
    #2
    test_shuf = shuf.loc[shuf['region_id'].isin(ytest), shuf.columns[1]].value_counts().reset_index()
    
    #3
    if len(test_shuf)>1:  # if some shuffles overlap the footprint
        no, yes = test_shuf.iloc[:,1]

    elif False in test_shuf.iloc[:,0]:  # no shuffles overlap the footprint
        no,yes = len(test_shuf), 0

    else:  # all shuffles overlap the footprint
        no,yes = 0, len(test_shuf)
    #4
    n_overlaps +=yes
    n_nooverlaps +=no
    
    #5
    return n_overlaps, n_nooverlaps


def make_2x2_shufs(test_tf, tf, ytest, shuf_file):
    
    """
    return counts of FP that overlaps an activity category (e.g. cis-only) 
    and counts of FP that overlaps the matched shuffled regions to that activity category (e.g. cis-only)
    """
    
    test = len(set(test_tf[tf]))
    if len(set(test_tf[tf]))>1:  # check that TF footprints at all. 

        # differential FP in yes test set?
        first_set = set(test_tf.loc[test_tf.region_id.isin(ytest), tf])
        
        if len(first_set)>1:

            b,a = test_tf.loc[test_tf.region_id.isin(ytest)].groupby(tf).count().reset_index().iloc[:, 1]
            
        elif False in first_set:
            b,a = len(ytest), 0  # all values are False in y set
        
        else:
            b,a = 0, len(ytest)
        
        c,d = make_shuf_vector(shuf_file, tf, ytest)
            
    else:
        a,b,c,d = 0,0,0,0
    
    return a,b,c,d


def run_OR(cl, col, tfs, fp, regions, re, shuf_file):
    
    out = os.path.join(re, f"{cl}-{col}_x_shuf_OR.tsv")

    ytest, ntest = get_fp_counts(col, regions)
    
    if os.path.exists(out) is True:
        
        comp = pd.read_csv(out, sep = '\t', usecols = ["comparison"])
        comparisons = list(set(comp["comparison"]))

    else:
        comparisons = []

                 
    for i, tf in enumerate(tfs):
        
        comparison = f"{col}_x_{tf}-shuf"
        
        if comparison not in comparisons: # check that you haven't run this alreayd

            test_tf = fp[["region_id", tf]].drop_duplicates()

            a,b,c,d = make_2x2_shufs(test_tf, tf, ytest, shuf_file)

            check_fp = a+b+c+d  # did footprint overlap any element?
            check_fp_n = a+c

            if check_fp>0 and check_fp_n>50:

                results = fet.get_2x2(a,b,c,d, comparison)
                results["arch_id"] = tf
                results["tested"] = col

                # write the results to a file. 
                if i==0:
                    results.to_csv(out, sep = '\t', mode="a", index = False)
                else:
                    results.to_csv(out, sep = '\t', mode="a", index = False, header = False)
        else:
            print('already compared', comparison)

    return out

def fdr_results(out):
    
    df = pd.read_csv(out, sep = '\t', nrows = 2)

    if "FDR_P" not in list(df):
    
        df = pd.read_csv(out, sep = '\t') # load full dataframe
        # get vector of p values
        df = df.loc[df["P"] !="P"]
        pvals = df["P"].astype(float)
        print(pvals)
    
        # run FDR correction at alpha = 0.05
        df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

        # add an asterisks column
        df["asterisks"] = None
        df.loc[df["reject_null"]== True, "asterisks"] = "*"
    
        # add negative log ten p
        df['-log10p'] = np.log10(df["FDR_P"])*-1
        

        
        # re-write file
        df.to_csv(out, sep = '\t', index = False)
        
    else:
        df = pd.read_csv(out, sep = '\t') # load full dataframe
        df = df.loc[df["P"] !="P"]
        df = df.drop_duplicates()
        df.to_csv(out, sep = '\t', index = False)
        print("already FDR corrected")
    
                
###    
#MAIN    
###


def main(argv):


    # get sample names
    path, region_file, region = split_filename.split_filename(REGIONS)

    # open FP and region data
    fp = pd.read_csv(FP_RE, sep = '\t')

    regions = pd.read_csv(REGION_ANNOT, sep ='\t',low_memory=False)


    tfs = list(fp)[4:]  #list of TF names to test

    print("\n\nN TF FPs to test", len(tfs), "\n\n")
    
    out = run_OR(CL, COL, tfs, fp, regions, RE, SHUF) # run enrichment

    fdr_results(out)

if __name__ == "__main__":
    main(sys.argv[1:])