"""
20220711
SarahFong

Returns
    (1) .tsv file w fold-change between observed counts 
        and median of expected counts from i shuffled iterations. 
        Also calculates fold-change distribution w/ bootstrap. 

    (2) .csv w expected and observed counts

Used for GWAS and eQTL enrichments.

Inputs
    Bedfile (str) - file w/ full path, quantify overlap w/ variant file
    Variant file (str) - file w full path, target for measuring overlap in .bed file
    Background file (str) - file w/ full path of regions to include in expected background. 
    Iters (int) - number of iterations to shuffle background file. 
    Number of threads (int) - n threads to run shuffles
    Outfile (str) - file w/ full path to write results to, write expected counts to 
    
"""

import os
import sys, traceback
import argparse
import csv
import datetime
import numpy as np
from functools import partial
from multiprocessing import Pool
import pandas as pd
sys.path.append(os.getcwd())
sys.path.append("/home/fongsl/.conda/envs/sfenv") # add environment w/ pybedtools
from pybedtools import BedTool


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enrichment of bed file w/ target file overlap.")

arg_parser.add_argument("bedfile", help='bed file to calculate enrichment in')

arg_parser.add_argument("varfile", help='variant file')

arg_parser.add_argument("bkgd", help='bkgd file')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of shuffled iterations; default=100')

arg_parser.add_argument("-n", "--num_threads", type=int,
                        help='number of threads; default=SLURM_CPUS_PER_TASK or 1')

arg_parser.add_argument("-o", "--outfile", type=str, default=None,
                        help="print expected counts to file")

args = arg_parser.parse_args()

# save parameters
BED_FILENAME = args.bedfile
VAR_FILENAME = args.varfile
BKGD_FILENAME = args.bkgd
OUT_FILENAME = args.outfile
ITERATIONS = args.iters


# calculate the number of threads
if args.num_threads:
    num_threads = args.num_threads
else:
    num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))

    
###
# functions to intersect bed files 
###


def calculateObserved(bed, gwas):
    
    """
    return count of overlaps between bed and target file. 
    
    input
        .bed file (pybedtools object) 
        target file (pybedtools object)
    
    
    """
    obs_sum = 0

    obs_intersect = bed.intersect(gwas, wo=True) # intersect obs regions w/ test file (e.g. enh x gwas snps)

    for line in obs_intersect:
        obs_sum += int(line[-1]) # sum snps in obs regions

    return obs_sum


def calculateExpected(bed, gwas, incl_file, i):
    """
    return expected count of overlaps of shuffled bed file (from a defined background) and target file

        - shuffled regions match on chromosome, non-overlapping

    input
        .bed file (pybedtools object) 
        target file (pybedtools object)
        incl_file (str?) - background file to include in shuffle
        
    method
        1. shuffle bed file in included background file
        2. intersect shuffle w/ target file
        3. count number of overlaps in shuffle. 
        
    output expected counts
    """
    
    exp_sum = 0

    #1
    rand_file = bed.shuffle(genome='hg38', incl=incl_file, chrom=True, noOverlapping=True) # shuffle obs regions
    
    #2
    exp_intersect = rand_file.intersect(gwas, wo=True) # intersect shuffled obs regions w/ test file
    
    #3
    for line in exp_intersect:
        exp_sum += int(line[-1]) # sum snps in shuffled exp regions

    print("exp-", exp_sum)
    
    return exp_sum


def calculateEmpiricalP(obs, exp_sum_list):
    """
    return two lists
        (1) info - vector w/  
                n_obs, 
                median_exp, 
                std, 
                fold-change  # calculated from the median of expected shuffle 
                p_val
                
        (2) fold_changes- vector expected fold changes (to calculate confidence interval)
        
    input
        observed overlap count (int)
        list of expected overlap counts (list of ints)
    
    method
        1. get median of expected overlap counts
        2. get standard deviation of expected overlap counts
        3. center expected overlap counts at median
        4. Sum the number of centered expected counts greater than observed centered count
            This is two tailed because it evaluates both sides of the distribution (w/ abs value). 
        5. calculate fold change as observed/ median expected w/ pseudo count
        6. calculate fold change of each "obs"/ expected w/ pseudo count
        7. calculate the p-value as count of equal or more extreme values than observed value
        8. return list of empirical info + fold changes
        
        
    
    """
    #1
    mu = np.median(exp_sum_list)  # median of exp.dist
    
    #2
    sigma = np.std(exp_sum_list)  # std
    
    #3
    dist_from_mu = [exp - mu for exp in exp_sum_list] # center the distribution 
    
    #4
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu)) # count values >= centered obs

    #5
    fold_change = (obs + 1.0) / (mu + 1.0) # fold change obs from median expected w pseudo count
    
    #6
    fold_changes = list((obs + 1.0) / (m + 1.0) for m in exp_sum_list) # fold change obs from /each exp w pseudo count
    
    #7
    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0)  # probability of observing obs-like value equal or more extreme in expected distribution
    
    #8
    info = [
            obs, 
            mu, 
            sigma, 
            fold_change, 
            p_val, 
            str(datetime.datetime.now())
            ]
    
    return info, fold_changes


def bootstrapCI(fold_changes_list): # bootstrap C.I.s
    
    """
    return the confidence intervals for list of obs/expected fold changes. 
    
    input
        list of floats - obs/exp fold changes, length N. Each index is the obs/each expected overlap) 
        
    method 
        1. get the length of the list
        2. get the median of the list
        3. set bootstrap parameters
        4. per iteration, randomly choose elements from the fold changes list w replacement
        5. append the median to the list of bootstrapped_medians
        6. turn medians into a dataframe
        7. calculate the delta distances from the population median. 
        8. sort from largest to smallest difference
        9. get 0.025 adn 0.975 quantile values of the delta median distribution. 
        10. calculate confidence intervals and the population median - quantile values

    """
    
    #1
    n = len(fold_changes_list) # size of distribution to bootstrap
    #2
    xbar = np.median(fold_changes_list) # get median fold-change from n shuffles (if using mean for empirical obs, use mean here)
    #3
    nboot = 10000 # resample 10000 times
    val = 0
    bs_medians = []
    #4
    while val < nboot:

        bs_dist = np.random.choice(fold_changes_list, replace = True, size = n)
        
        #5
        bsmedian = np.median(bs_dist)
        bs_medians.append(bsmedian)
        val +=1
    #6
    bs = pd.DataFrame(data = bs_medians, index = np.arange(nboot), columns = ["bs_medians"]) # make dataframe of bootstraps

    #7
    bs["deltas"] = bs.bs_medians - xbar

    #8
    bs = bs.sort_values(by = "deltas", ascending= False)
    
    #9 # get 95th CI
    low = bs.deltas.quantile(0.025) 
    high = bs.deltas.quantile(0.975)

    #10
    ci = xbar - [high, low]
    print("CI", low, high, ci, xbar)
    return ci


def write_csv(file, info_list):
    """
    write list to csv file
    
    input 
        csv filename (str)
        info_list (list)

    method
        1. open file
        2. call csv_writer
        3. write row
        4. close csv file
    """
    #1
    csvfile = open(file, "a", newline='')  # open the .csv in append mode
    #2
    csv_writer = csv.writer(csvfile, delimiter = '\t') # call the writer function
    #3
    csv_writer.writerow(info_list) # write the info list to the csv 
    #4
    csvfile.close() # close the csv

    
    
def main(argv):
    
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())))

    sid = (BED_FILENAME.split("/")[-1]).split(".bed")[0] # name
    target = (VAR_FILENAME.split("/")[-1]).split(".bed")[0]
    ntotal_target = sum(1 for line in open(VAR_FILENAME))
    print("\n\noutfile:", OUT_FILENAME)

    # observed overlap counts
    obs_sum = calculateObserved(BedTool(BED_FILENAME), BedTool(VAR_FILENAME))
    print("\nobs_sum", obs_sum)

    # create pool and run expected shuffle counts in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, BedTool(BED_FILENAME), BedTool(VAR_FILENAME), BKGD_FILENAME)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])

    # wait for results to finish before calculating p-value
    pool.close()
    pool.join()

    # calculate empirical p value
    obs_p, fold_changes= calculateEmpiricalP(obs_sum, exp_sum_list)

    # bootstrap CI
    ci = bootstrapCI(fold_changes)

    # gather, write data
    obs_p.extend(ci)
    info =[ITERATIONS, sid, target, ntotal_target]
    obs_p.extend(info)

    header = ['Observed', 'Expected', 'StdDev',
              'FoldChange', 'p-value','date_time',
              "ci_975", "ci_025", "iters", "sid", "target", "ntotalvars_in_target"
             ]
    print(header)
    print(obs_p)

    
    # if zero lines, write header
    if os.path.exists(OUT_FILENAME) is False:
        write_csv(OUT_FILENAME, header)

    write_csv(OUT_FILENAME, obs_p)
    
    
    # make an expected file
    exp_file = os.path.join(os.path.splitext(OUT_FILENAME)[0] + "-exp.csv")

    exp_header = ["iters", "sid", "target", "ntotalvars_in_target", "obs", "expn"]
    
    # copy the important information about the run
    exp_info = info.copy()
    exp_info.extend([obs_sum])
    
    # extend w/ fold changes
    exp_info.extend(exp_sum_list)
    
    # write the expected fold_changes 
    if os.path.exists(exp_file) is False:

        write_csv(exp_file, exp_header)
    write_csv(exp_file, exp_info)

    
if __name__ == "__main__":
    main(sys.argv[1:])