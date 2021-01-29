#!/usr/bin/env bash
#SBATCH --partition=gpu_8
#SBATCH --ntasks=1
#SBATCH --time=30:00:00
#SBATCH --mem=40gb
#SBATCH --gres=gpu:2
#SBATCH --job-name=modeltraining
#SBATCH --mail-type=ALL

work_dir=/pfs/work7/workspace/scratch

source ${work_dir}/db0052-conda-0/conda/bin/activate ma_deepsv

python ${work_dir}/db0052-ma_deepsv-0/Model/train_ray_tune.py --model $1

