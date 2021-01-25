#!/usr/bin/env bash
#SBATCH --partition=gpu_8
#SBATCH --ntasks=1
#SBATCH --time=30:00:00
#SBATCH --mem=10gb
#SBATCH --gres=gpu:1
#SBATCH --job-name=modeltraining
#SBATCH --mail-type=ALL

work_dir=/pfs/work7/workspace/scratch

source ${work_dir}/db0052-conda-0/conda/bin/activate ma_deepsv

python ${work_dir}/db0052-ma_deepsv-0/Model/train.py --val_epoch $1 --lr $2 --model $3

