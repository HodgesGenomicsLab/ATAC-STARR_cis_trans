#!/usr/bin/env python
# coding: utf-8

# In[1]:


from functools import partial
import glob
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import pybedtools as pbt
import re
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")
import config_readwrite as crw
import chr_functions
import split_filename


# In[2]:


name = "/data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/config"
config, configfile_name = crw.read_config(name)

SHUF = config["SHUFFLES"]["shuf-all"]
REGIONS = config["CIS_TRANS"]["regions"]

ID_TAG = config["TF_FOOTPRINTING_JASPAR"]["ID_TAG"]

PATH = config["TF_FOOTPRINTING"]["PATH"]
RE = config["TF_FOOTPRINTING"]["results"]

path, region_file, region = split_filename.split_filename(REGIONS)


# 0. string split to get TF name
# 1. intersect regions w footprint files
# 2. create matrix w/ region_id (should be 17605 regions x 693 archetypes) 

# # functions

# ## get the TF id

# In[5]:


def get_arch_id(filename):
    """
    return archetype ID w/ chr/int split on first two values
    
    do string formatting
    """
    
    match = re.match(r"([a-z]+)([0-9]+)", filename.split("_")[0], re.I)  
    
    if match:
        items = match.groups()
    arch_id = ("".join(items))
    
    return arch_id


# ## get TF FP file to intersect 

# In[6]:


def get_tf_bound_file(f, filename, cell_line):

    tf_bound_file = os.path.join(f, "beds", f"{filename}_{cell_line}_1.8-filter_bound.bed")
    
    return tf_bound_file


# ## build a reference dict w/ all TF FP files for that species

# In[7]:


def get_tf_paths(path, cell_line, id_tag):
    """
    return dictionary w/ key = archetype id and value = tuple of TF name and full file path.  
    
    1. get all the file names
    2. get arch_id name - string split function
    3. get TF name - str split
    4. get TF FP bound file
    5. make dictionary
    6. return dictionary
    """
    tf_dict ={}
    
    #1
    fs = glob.glob(os.path.join(path, f"*{id_tag}*"))
    print(len(fs))
    for f in fs:
        #2
        filename = os.path.split(f)[1]  # get the file name
        #3
        if "MA" in id_tag:
            arch_id = filename.upper()  
            tf_name = filename.split("_")[0]
        else:
            arch_id=get_arch_id(filename) #FOR ARCHTYPES ONLY
            
            #4
            tf_name = "_".join((filename.split(arch_id)[1]).split("_")[:-1])  # get TF name by splitting string. 

        
        #5
        tf_bound_file = get_tf_bound_file(f, filename, cell_line)  # get path to TF bound.bed file
        #6
        tf_dict[arch_id] = (tf_name, tf_bound_file)
    
    return tf_dict


# ## Intersect regions w FP

# In[8]:


def intersect_w_fp(regions, fp_name, fp_bed, outdir, cell_line, shuf):
    
    if shuf is True:
        cell_line = f"shuf-{cell_line}"
    
    
    if os.path.exists(os.path.join(outdir, cell_line)) is False:  # make a directory
        os.mkdir(os.path.join(outdir, cell_line))
        
    out = os.path.join(outdir, cell_line, f"{cell_line}-{fp_name}.bed")
    
    if os.path.exists(out) is False:  # if intersection hasn't already been done. 
        
        bedregions = pbt.BedTool(regions)
    
        bedregions.intersect(pbt.BedTool(fp_bed), wo = True).saveas(out)
        print(out)
    
    return out
    


# ## combine all intersections 

# In[9]:


def concat_intersections(path, cell_line, out, id_tag, shuf):

    if shuf is True:
        cell_line = f"shuf-{cell_line}"

    # make query to concatenate all individual TF FP intersections
    cat_query = os.path.join(path, cell_line, f"{cell_line}*{id_tag}*.bed" )
    
    if os.path.exists(out) is False:
        cmd = f"cat {cat_query} > {out}"
        subprocess.call(cmd, shell=True)

    return out


# ## make matrix from combined intersections 

# In[10]:


"""
load region x TF FP intersection
pivot table to make matrix w/ counts of FP overlap
"""

def make_matrix(concat, matrix_outfile):
    
    if os.path.exists(matrix_outfile) is False:
    
        cols= [
            "#chr",
            "start", "end", "region_id",
            "#chr_tf",
            "start_tf",
            "end_tf",
            "tfid",
            "score", 
            "strand",
            "overlap"
        ]
        
        df = pd.read_csv(concat, sep='\t', header =None, 
                         usecols=[0,1,2,3,4,5,6,7,8,9,21],
                         names = cols
                        )
        
        df=df.drop_duplicates().reset_index()
        
        # calculate the length of the TF footprint
        df["len"] = df["end_tf"] - df['start_tf']
        
        # CLEAN region must overlap at least 1/2 of the TF footprint motif len
        df = df.loc[df["overlap"]>=(df["len"]/2)]
        
        # write cleaned up concat
        concat_out = os.path.splitext(concat)[0] + "_clean.bed"
        
        #write file to liftover
        liftover_out = os.path.splitext(concat)[0] + "_for_liftover.bed"
    
        # add an id for the liftover
        df["tfid2"] = df["#chr_tf"] + ":" + df["start_tf"].map(str) + ":" + df["end_tf"].map(str) + "_" + df["tfid"]

        # pick cols for liftover bed
        liftover_cols= ["#chr_tf","start_tf","end_tf", "tfid2", "region_id"]
            
        df.to_csv(concat_out, sep='\t', index=False)
        df[liftover_cols].to_csv(liftover_out, sep='\t', index=False)
        
        # groupby for matrix

        grouped = df.groupby(["region_id", "tfid"])['overlap'].count().reset_index()

        matrix = grouped.pivot(index="region_id", columns = "tfid", values = "overlap")
        matrix.to_csv(matrix_outfile, sep='\t')
    
    else:
        matrix = pd.read_csv(matrix_outfile, sep='\t')
    
    return matrix


# # pipeline function 

# In[11]:


def pipeline(outfile, cl, id_tag, regions, re, path, fp_path, shuf, matrix_outfile):
    
#if os.path.exists(outfile) is False:
    tfs = get_tf_paths(fp_path, cl, id_tag) ## get all the tf.bed paths

    for k, vs in tfs.items():  # per TF footprint file
        tfname, FP_bed = vs[0], vs[1]

        out = intersect_w_fp(regions, k, FP_bed, re, cl, shuf)  # intersect regions w fp

    concat = concat_intersections(path, cl, outfile, id_tag, shuf)
    
    matrix = make_matrix(concat, matrix_outfile)
    
    return matrix


# # main

# In[ ]:


results = {}
CLS =["GM12878", "LCL8664"]
for CL in CLS:
   
    # get config for that CL
    FP_PATH = config["TF_FOOTPRINTING_JASPAR"][f"{CL}_bindetect"]

    FP_RE = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["FP"] # write
    MATRIX = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["matrix"] # write

    FP_RE_SHUF = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["FP"] # write
    SHUF_MATRIX = config[f"TF_FOOTPRINTING_JASPAR_{CL}"]["matrix_shuf"] # write
    
    SHUFFLED = False
    cl_matrix = pipeline(FP_RE, CL, ID_TAG, REGIONS, RE, PATH, FP_PATH, SHUFFLED, MATRIX)
    results[CL] = cl_matrix # add to results dict
    
    SHUFFLED = True
    cl_matrix_shuf = pipeline(FP_RE_SHUF, CL, ID_TAG, SHUF, RE,PATH, FP_PATH, SHUFFLED, SHUF_MATRIX)
    results[f"shuf-{CL}"] = cl_matrix_shuf # add to results dict


# reference this:
# https://stackoverflow.com/questions/44578571/intersect-two-boolean-arrays-for-true

# # subtract hu frm rh matrix to get differential FP matrix

# In[ ]:


hu = results["GM12878"]
rh = results["LCL8664"]


# In[ ]:




