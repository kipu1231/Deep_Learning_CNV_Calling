import os
import numpy as np
import subprocess
import shlex
import argparse
from scipy.stats import poisson, nbinom
import matplotlib.pyplot as plt
import pysam
import math


def get_bam_seq_info(bam_path):
    #p2 = subprocess.call(shlex.split("samtools view -c -f 1 %s" % (bam_path)))
    p2 = subprocess.call(shlex.split("samtools idxstats %s" % (bam_path)))

    # popen_bam = os.popen('samtools idxstats ' + str(bam_path))
    # sys_line = popen_bam.readlines()
    # popen_bam.close()
    #
    # for dep in sys_line:
    #     print(dep)
    #     chr = dep.strip('\n').split('\t')[0]
    #     chr_len = dep.strip('\n').split('\t')[1]
    #     mapped_reads = dep.strip('\n').split('\t')[2]
    #     print(chr)
    #     print(chr_len)
    #     print(mapped_reads)



def get_mapped_reads(bam_path):
    popen_bam = os.popen('samtools depth -a ' + str(bam_path))
    sys_line = popen_bam.readlines()
    popen_bam.close()

    dep_count = 0
    gen_len = 0
    h = 0
    for dep in sys_line:
        dep_count = dep_count + int(dep.strip('\n').split('\t')[2])
        h = int(dep.strip('\n').split('\t')[1])
        if h > gen_len:
            gen_len = h
    print(dep_count)
    print(gen_len)
    return dep_count, gen_len


def get_chr_stats(dep_count, gen_len, bin_size):
    num_rd_per_base = dep_count
    bin_size = bin_size
    len_genome = gen_len
    d = 3

    mu = (num_rd_per_base / len_genome) * bin_size
    var = mu / (d - 1)

    med = poisson.median(mu) / bin_size
    l_p, u_p = poisson.interval(0.9, mu)
    lower_percentile = l_p/bin_size
    upper_percentile = u_p/bin_size

    return med, lower_percentile, upper_percentile


if __name__ == '__main__':
    bam_path = '/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/Preprocessing/NA18511/chr_1/NA18511_chr_1.bam'
    #get_mapped_reads(bam_path)
    get_bam_seq_info(bam_path)


    bamfile = pysam.AlignmentFile(bam_path, 'rb')

    # depth = bamfile.count_coverage('1',60006,60100)
    # np_depth = np.array(list(depth))
    # sum_depth = np_depth.sum(axis=0)
    # print(depth)
    # ave_depth = sum_depth.sum() / (60100 - 60006 + 1)
    # print(ave_depth)
    # if ave_depth > 0: ave_depth = math.log2(ave_depth / float(cova))

    # num_reads = 11758289
    # num_rd_per_base = 1128005051
    #
    # bin_size = 50
    # len_genome = 249250621
    # d = 3
    #
    # mu = (num_rd_per_base/len_genome)*bin_size
    # var = mu/(3-1)
    # med = poisson.median(mu)/50
    # u_p, o_p = poisson.interval(0.9, mu)
    #
    # print(u_p / 50)
    # print(med)
    # print(o_p/50)


