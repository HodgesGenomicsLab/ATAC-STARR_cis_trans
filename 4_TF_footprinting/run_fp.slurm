#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=05:30:00
#SBATCH --mem=16GB

#SBATCH --output=/data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/SF_hu-rhe/activity_bin/

#SBATCH --job-name=CON_ACC
#SBATCH --account=capra_lab
#SBATCH --output=%J.stderr


ml load Anaconda3
source activate /home/fongsl/.conda/envs/sfenv/

python /data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/TF_FP/FP_enrichment.py ${1} ${2}