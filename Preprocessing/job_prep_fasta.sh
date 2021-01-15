#!/usr/bin/env bash
#SBATCH --partition=single
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=60gb
#SBATCH --job-name=vcfpreprocessing
#SBATCH --mail-type=ALL

work_dir=/pfs/work7/workspace/scratch

source ${work_dir}/db0052-conda-0/conda/bin/activate ma_deepsv

python ${work_dir}/db0052-ma_deepsv-0/data/preprocess_fasta.py
