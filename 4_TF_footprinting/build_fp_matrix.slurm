#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=64GB

#SBATCH --output=/data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/SF_hu-rhe/activity_bin/

#SBATCH --job-name=CON_ACC
#SBATCH --account=capra_lab
#SBATCH --output=%J.stderr


ml load Anaconda3

source /home/${USER}/.bashrc
source activate /dors/capra_lab/users/fongsl/sfenv

python /data/hodges_lab/ATAC-STARR_B-cells/bin_human-evolution/TF_FP/build_fp_matrix.py 