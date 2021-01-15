import os
import numpy as np
import subprocess
import shlex
import argparse


def get_bam_name(sample_id):
    path_bam = None

    info_samples = np.array([['NA19625', 'NA19625.wgs.ILLUMINA.bwa.ASW.high_cov_pcr_free.20140203.bam'],
                            ['NA19648',  'NA19648.wgs.ILLUMINA.bwa.MXL.high_cov_pcr_free.20140203.bam'],
                            ['NA18525', 'NA18525.wgs.ILLUMINA.bwa.CHB.high_cov_pcr_free.20140203.bam'],
                            ['NA19017', 'NA19017.wgs.ILLUMINA.bwa.LWK.high_cov_pcr_free.20140203.bam'],
                            ['NA18511', 'NA18511_chr_1_sorted.bam'],
                            ['NA20502',  'NA20502.wgs.ILLUMINA.bwa.TSI.high_cov_pcr_free.20140203.bam']])

    for data in info_samples:
        if sample_id == data[0]:
            path_bam = data[1]

    return path_bam


def preprocess_bam_per_chromosome(args, chr_id):
    'create directory to save bam file'
    save_dir_bam = os.path.join(args.sample_id, 'chr_' + chr_id)

    if not os.path.exists(save_dir_bam):
        os.makedirs(save_dir_bam)

    name_bam_out = os.path.join(save_dir_bam, args.sample_id + "_chr_" + chr_id + ".bam")

    'extract bam file for given chromosome and index it'
    p1 = subprocess.call(shlex.split("%s view -b %s %s > %s" % (args.samtools, args.bam_fn, chr_id, name_bam_out)),
                         stdout=open(name_bam_out, 'w'))
    p2 = subprocess.call(shlex.split("%s index %s" % (args.samtools, name_bam_out)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Preprocess Bam Files")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools'")

    parser.add_argument('--work_dir_bwcluster', type=str, default="/pfs/work7/workspace/scratch",
                        help="Path of work space directory on bwcluster")
    parser.add_argument('--ws_dir', type=str, default="db0052-ma_deepsv-0/data",
                        help="Path of workspace directory")
    parser.add_argument('--work_dir', type=str, default="",
                        help="Path of working directory")

    parser.add_argument('--sample_id', type=str, default="NA18525",
                        help="ID of sample to be processed")
    parser.add_argument('--bam_fn', type=str, default='default.bam',
                        help="Sorted bam file input")
    parser.add_argument('--bam_bai_fn', type=str, default='default.bam.bai',
                        help="Index for sorted bam file input")

    parser.add_argument('--chromosome', type=str, default="all",
                        help="Specify if bam files should be generated for all or a specific chromosome")

    args = parser.parse_args()

    'create directory to save data ./db0052-ma_deepsv-0/data/NA20502'
    args.work_dir = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)

    if not os.path.exists(args.work_dir):
        os.makedirs(args.work_dir)

    'get name of bam and bai file for sample id'
    name_bam = get_bam_name(args.sample_id)

    if name_bam is not None:
        args.bam_fn = os.path.join(args.sample_id, name_bam)
        args.bam_bai_fn = os.path.join(args.sample_id, name_bam + ".bai")

    #print(args.bam_bai_fn)
    #print(args.bam_fn)

    'Start Preprocessing'
    l1 = np.arange(1, 23, 1)
    l1 = l1.astype(str)

    if args.chromosome == 'all':
        for i in l1:
            preprocess_bam_per_chromosome(args, i)

    elif args.chromosome in l1:
        preprocess_bam_per_chromosome(args, args.chromosome)

    else:
        print('chromosome does not exist')
