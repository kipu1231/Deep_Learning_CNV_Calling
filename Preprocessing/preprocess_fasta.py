import argparse
import subprocess
import shlex
import os

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

    args = parser.parse_args()

    args.work_dir = os.path.join(args.work_dir_bwcluster, args.ws_dir)
    fasta_fn = os.path.join(args.work_dir, "Homo_sapiens.GRCh37.dna.primary_assembly.fa")

    'Index Fasta'
    p2 = subprocess.call(shlex.split("%s faidx %s" % (args.samtools, fasta_fn)))
