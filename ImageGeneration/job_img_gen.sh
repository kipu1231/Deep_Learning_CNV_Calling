#!/usr/bin/env bash
#SBATCH --partition=single
#SBATCH --ntasks=1
#SBATCH --time=15:00:00
#SBATCH --mem=60gb
#SBATCH --job-name=preprocessing
#SBATCH --mail-type=ALL

work_dir=/pfs/work7/workspace/scratch

source ${work_dir}/db0052-conda-0/conda/bin/activate ma_deepsv

python ${work_dir}/db0052-ma_deepsv-0/data/Image_Generation.py --sample_id $1

#change for different preprocessing steps
#python ${work_dir}//db0052-ma_deepsv-0/data/preprocess_bam.py --sample_id $1 --chromosome $2
#python ${work_dir}//db0052-ma_deepsv-0/data/preprocess_vcf.py --prepro_type $1
#python ${work_dir}//db0052-ma_deepsv-0/data/preprocess_fasta.py