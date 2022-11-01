#!/usr/bin/env python
# coding: utf-8

# In[6]:


import glob
import numpy as np
import os
import pandas as pd
import pybedtools as pbt
import subprocess
import sys

sys.path.append(os.getcwd())
import UKBB_traits

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import config_readwrite as crw

import split_filename

name = "/data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/config"
config, configfile_name = crw.read_config(name)
PATH =config["UKBB"]["path"] 

LD_PY = config["GWAS_CATALOG"]['bin_ld_expand']
ANC = config["GWAS_CATALOG"]["ANC"]
NOMINAL_P = float(config["GWAS_CATALOG"]["NOMINAL_P"])


# # functions 

# In[9]:


def get_meta(f):
    
    """
    return meta variable (bool) reflecting whether the GWAS has cross-population meta-analysis data
    
    input
        f (str) - GWAS file name
        
    method
        read columns as pandas dataframe
        assess whether pval_meta field is in df columns
        assign meta variable as True or False
        
    """
    df_ = pd.read_csv(f, sep = '\t', nrows=2)
    if "pval_meta" in list(df_):
        META = True
    else:
        META=False

    return META


# In[10]:


def keep_rearrange_cols(meta):
    """
    return lists of columns to keep and how to order for .bed file
    
    META or EUR? if there is a meta analysis across populations, prioritize meta p-value. 
        - some trait only have EUR p-values, so I have to use those in some instances. 
    """
    if meta is True:
        POP = "meta"
    else:
        POP = "EUR"
    
    print("pvalue to use", POP)
    
    keep = [
            'chr',
            'pos',
            'ref',
            'alt',

            f'beta_{POP}',
            f'se_{POP}',
            f'pval_{POP}',
            
            ]
    
    rearrange =[
                '#CHR',
                'pos',
                "POS+1",
                'ref',
                'alt',
                f'beta_{POP}',
                f'se_{POP}',
                f'pval_{POP}',
            
                ]
    PVAL_COL =  f'pval_{POP}'
    
    return keep, rearrange, PVAL_COL


# In[11]:


def clean_pval(gwas, pval_col, nominal_p, outfile):
    """
    return rows w genome-wide significant p-values
    
    input 
        gwas: pandas dataframe, gwas catalog bed file
        nominal_p: (float) nominal genome-wide significance p-value
        
    output
        cleanp: pandas dataframe filtered on genome-wide significant p-value threshold (smaller than nominal p)
    """
    print('cleaning GWAS. Keep var w/ pval<', pval_col, nominal_p)

    cleanp = gwas.loc[gwas[pval_col].astype(float)<=np.log10(nominal_p)]
    
    cleanp.drop_duplicates().to_csv(outfile, sep='\t', index=False) # write file
    
    print(gwas.shape, cleanp.shape, "n insig", gwas.shape[0]-cleanp.shape[0])
    
    return cleanp                  


# In[12]:


def ld_expand(clean_gwas, anc, ld_py, outfile):
    """
    write LD expand using 1Kg ancestry references
    
    inputs
        clean_gwas (str)- path to clean gwas file
        anc (str) - ancestry to LD expand on in 1KG
        ld_py (str)- python LD expand script
        outfile (str) - file to write LD expanded results to
        
    """
    print('ld expanding', clean_gwas)
    
    cmd = f"python {ld_py} {anc} {clean_gwas} > {outfile}"
    subprocess.call(cmd, shell = True)


# In[13]:


def drop_duplicates(out_ld, out_uniq):
    
    """
    drop duplicates from LD-expanded GWAS file
    
    input
        out_ld (str) - path to original gwas file, before cleaning, before LD expanding. 
        out_unq (str) - file to write w/ uniq entries. 
        
    method
        load file as panda dataframe, drop duplicates. Probably could do this in bash. 
    """
    print("dropping duplicates from", out_ld)
    
    # load the LD-expanded dataframe
    df= pd.read_csv(out_ld, sep ='\t', header =None, low_memory=False)
    uniq = df.drop_duplicates() # drop duplicates
    
    print(df.shape[0], uniq.shape[0])
    
    # save the dataframe
    uniq.to_csv(out_uniq, sep = '\t', header=False, index=False)

    
def format_df(f, keep, rearrange):
    df = pd.read_csv(f, sep = '\t', low_memory=False)


    path, filename, sid = split_filename.split_filename(f)

    # pick columns to keep
    df = df[keep]

    # remove chr formatting issues
    df = df.loc[~df["chr"].isna()]
    df = df.loc[~df["chr"].map(str).str.contains(";")]
    df = df.loc[~df["chr"].map(str).str.contains("x")]

    # make into a bed file
    df["POS+1"] = df["pos"].map(int) +1
    df["#CHR"] = "chr"+df["chr"].map(str)


    # rearrange columns
    df = df[rearrange].drop_duplicates()
    
    return df

def main(argv):

    dl_dict = UKBB_traits.dl2_dict() # second round of traits


    for f in dl_dict.values():

        f=(f.strip(".bgz")).split("/")[-1]
        path, filename, sid = split_filename.split_filename(f)
        path = os.getcwd()
        f = os.path.join(PATH, f)


        # outfiles to write
        out_clean = os.path.join(path, (sid + "clean.bed"))
        out_ld = os.path.join(path, (sid + "clean_LD_exp.bed"))
        out_uniq = os.path.join(path, (sid + "clean_LD_exp_uniq.bed")) # write this

        if os.path.exists(out_uniq) is False:

            if os.path.exists(out_ld) is False:

                if os.path.exists(out_clean) is False:

                    META = get_meta(f)  # is there a meta analysis for this GWAS?

                    keep, rearrange, PVAL_COL = keep_rearrange_cols(META)
                    
                    #format the dataframe
                    df=format_df(f, keep, rearrange)

                    # clean p-value
                    clean_pval(df, PVAL_COL, NOMINAL_P, out_clean)
                    
                    # LD expand
                    ld_expand(out_clean, ANC, LD_PY, out_ld)
                    
                    # drop duplicates
                    drop_duplicates(out_ld, out_uniq)
                    
                    os.remove(out_ld)
                    os.remove(out_clean)
                    os.remove(f)
                   
                else:
                    print("cleaned, need to LD expand")

                    # LD expand
                    ld_expand(out_clean, ANC, LD_PY, out_ld)
                    
                    drop_duplicates(out_ld, out_uniq)
                    
                    os.remove(out_clean)
                    os.remove(out_ld)
            else:
                print("LD expanded, but need to remove duplicates")
                
                # remove duplicates
                drop_duplicates(out_ld, out_uniq)
                
                os.remove(out_ld)
                
        else:
            print("clean, LD expanded, no duplicates, let's go!")

if __name__ == "__main__":
    main(sys.argv[1:])




