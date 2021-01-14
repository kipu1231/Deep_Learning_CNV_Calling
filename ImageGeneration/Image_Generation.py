import argparse
import subprocess
import shlex
import os
import pysam
#from ImageGeneration.Variant_Filtering import *
#from ImageGeneration.Input_Image import *
#from ImageGeneration.Input_Image_cnnCNV import *
from Image_Visualization import *
import random


def prepare_variants(args, vcf_path, chr_id):
    'open csv and retrieve necessary information'
    cnv_file = open(vcf_path, 'r') #csv with preprocessed CNVs
    vcf_detail = []
    for temp in cnv_file:
        tmp = temp.split('\t')
        chr_id = tmp[0]
        if chr_id == chr_id:
            scan_l_pos = int(tmp[1])
            scan_r_pos = int(tmp[2])
            sv_type = tmp[3]
            vcf_tuples = (chr_id, scan_l_pos, scan_r_pos, sv_type)
            vcf_detail.append(vcf_tuples)
    return vcf_detail


def prepare_non_variants(args, fasta_path, chr_id):
    'add non CNV sites to vcf'
    vcf_detail_non_cnv = []

    'get required information for random sampling'
    #avg_chr_len = get_avg_len_of_chr(chr)
    img_save_dir = get_img_save_dir(args, chr_id)

    'calculate number of images to be generated to balance training set'
    img_to_be_gen = len(os.listdir(img_save_dir))
    print("Goal image number: " + str(img_to_be_gen))
    img_to_be_gen = int((img_to_be_gen/10)*1.25)

    sample_l_pos = random.sample(range(0, int(args.chr_len)), img_to_be_gen)
    #sample_l_pos = random.sample(range(0, int(args.chr_len)), 100)

    for pos in sample_l_pos:
        scan_l_pos = pos
        scan_r_pos = pos + 500
        if check_n(fasta_path, chr_id, scan_l_pos, scan_r_pos):
            if check_vcf(args, chr_id, scan_l_pos, scan_r_pos):
                sv_type = "NOCNV"
                vcf_tuples = (chr_id, scan_l_pos, scan_r_pos, sv_type)
                vcf_detail_non_cnv.append(vcf_tuples)

    print(len(vcf_detail_non_cnv))

    return vcf_detail_non_cnv


def generate_training_img(vcf_tuples, bam_path):
    'filter variant using kmeans --> method needs to be reworked'
    # filter_cnv(vcf_detail, sam_file)

    # if args.mode == 'deepsv':
    #     sam_file = pysam.AlignmentFile(bam_path, "rb")
    #     num_cnvs = len(vcf_tuples)
    #     cnv_np_ar = np.array(vcf_tuples)
    #
    #     for i in range(num_cnvs):
    #         svtype = cnv_np_ar[i, 3]
    #         chr = cnv_np_ar[i, 0]
    #         clip_record = get_clip_num(sam_file, cnv_np_ar[i, 0], int(cnv_np_ar[i, 1]), int(cnv_np_ar[i, 2]), cnv_np_ar[i, 3])
    #         clip_dict_record = dict(clip_record)
    #
    #         gene_pic = gene_point_pic(cnv_np_ar[i, 0], int(cnv_np_ar[i, 1]), int(cnv_np_ar[i, 2]))
    #
    #         for each_pic in gene_pic:
    #             pile_record = pipeup_column(sam_file, each_pic[0], each_pic[1], each_pic[2])
    #             deletion_length = each_pic[2] - each_pic[1]
    #             draw_pic(clip_dict_record, pile_record, each_pic[1], deletion_length, args, chr, svtype)
    #
    # elif args.mode == 'cnncnv':
    #     generate_cnn_cnv_images(args, vcf_tuples, bam_path)

    if args.mode == 'newvis':
        visual_training_img(args, vcf_tuples, bam_path)


def check_vcf(args, chr_id, l_pos, r_pos):
    vcf_fn = get_vcf_path(args)
    vcf = pysam.VariantFile(vcf_fn, 'r')

    f = vcf.fetch(str(chr_id), int(l_pos), int(r_pos))
    for read in f:
        return False
    return True


def get_vcf_csv_path(args):
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        vcf_fn = os.path.join(bw_path, str(args.sample_id) + "_cnv_list.csv")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        vcf_fn = os.path.join(local_path, str(args.sample_id) + "_cnv_list.csv")

    return vcf_fn


def get_vcf_path(args):
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        vcf_fn = os.path.join(bw_path, str(args.sample_id) + "_all_sv.vcf.gz")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        vcf_fn = os.path.join(local_path, str(args.sample_id) + "_all_sv.vcf.gz")

    return vcf_fn


def get_bam_path(args, chr_id):
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        bam_fn_loc = os.path.join(bw_path, "chr_" + str(chr_id))
        bam_fn = os.path.join(bam_fn_loc, str(args.sample_id) + "_chr_" + str(chr_id) + ".bam")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        bam_fn_loc = os.path.join(local_path, "chr_" + str(chr_id))
        bam_fn = os.path.join(bam_fn_loc, str(args.sample_id) + "_chr_" + str(chr_id) + ".bam")

    return bam_fn


def get_fasta_path(args):
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir)
        fasta_fn = os.path.join(bw_path, "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
    else:
        fasta_fn = os.path.join(args.data_dir, "Homo_sapiens.GRCh37.dna.primary_assembly.fa")

    return fasta_fn


def get_img_save_dir(args, chr_id):
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        bam_fn_loc = os.path.join(bw_path, "chr_" + str(chr_id))
        img_path = os.path.join(bam_fn_loc, "img_vis_50")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        img_path = os.path.join(local_path,"chr_" + str(chr_id), "img_deepsv")

    return img_path


def get_avg_len_of_chr(chr_id):
    info_avg_chr_len = np.array([['1', 248956422], ['2', 242193529], ['3', 198295559],['4', 190214555], ['5', 181538259],
                                ['6', 170805979], ['7', 159345973], ['8', 145138636], ['9', 138394717], ['10', 133797422],
                                 ['11', 135086622], ['12', 133275309], ['13', 114364328], ['14', 107043718],
                                 ['15', 101991189], ['15', 90338345], ['17', 81195210], ['18', 78077248], ['19', 58617616],
                                 ['20', 63025520], ['21', 46709983], ['22', 50818468]])

    for data in info_avg_chr_len:
        if chr_id == data[0]:
            avg_len = data[1]

    return avg_len


def get_mapped_reads(bam_path):
    open_bam = os.popen('samtools depth -a ' + str(bam_path))
    sys_line = open_bam.readlines()
    open_bam.close()

    dep_count = 0
    gen_len = 0
    h = 0
    for dep in sys_line:
        dep_count = dep_count + int(dep.strip('\n').split('\t')[2])
        h = int(dep.strip('\n').split('\t')[1])
        if h > gen_len:
            gen_len = h

    return dep_count, gen_len


def get_cov_info(bam_path, chr_id):
    popen_bam = os.popen('samtools idxstats ' + str(bam_path))
    sys_line = popen_bam.readlines()
    popen_bam.close()

    mapped_reads = 0
    chr_len = 0

    for dep in sys_line:
        if str(dep.strip('\n').split('\t')[0]) == str(chr_id):
            chr_len = int(dep.strip('\n').split('\t')[1])
            mapped_reads = int(dep.strip('\n').split('\t')[2]) * 100

    return mapped_reads, chr_len


def check_n(fasta_path, chr_id, l_pos, r_pos):
    fastafile = pysam.FastaFile(fasta_path)
    seq = fastafile.fetch(str(chr_id), l_pos, r_pos)

    if 'N' in seq:
        return False

    else:
        return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate images")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools'")

    parser.add_argument('--bwcluster', type=bool, default=True,
                        help="Define if script is run on bwcluster or not")
    parser.add_argument('--work_dir_bwcluster', type=str, default="/pfs/work7/workspace/scratch",
                        help="Path of work space directory on bwcluster")
    parser.add_argument('--ws_dir', type=str, default="db0052-ma_deepsv-0/data",
                        help="Path of workspace directory")
    parser.add_argument('--work_dir', type=str, default="/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/ImageGeneration",
                        help="Path of working directory")
    parser.add_argument('--data_dir', type=str, default="/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/Preprocessing",
                        help="Path of data directory")

    parser.add_argument('--sample_id', type=str, default="NA18511",
                        help="ID of sample to be processed")
    parser.add_argument('--chr_id', type=str, default="all",
                        help="Specify chr ID or generate for all")
    parser.add_argument('--chr_len', type=int, default=1,
                        help="Length of chromosome")
    parser.add_argument('--rd_count', type=int, default=1,
                        help="Read depth count per base of chromosome")

    parser.add_argument('--image_gen_type', type=str, default='del',
                        help="Generate images only for: del, dub or non_cnv")
    parser.add_argument('--mode', type=str, default='newvis',
                        help="Generate training images according to deepsv or cnncv")
    parser.add_argument('--img_size', type=int, default=50,
                        help="Determine image size")

    args = parser.parse_args()

    'get Path to CNV file'
    vcf_fn = get_vcf_csv_path(args)

    'get Path to fasta file'
    fast_fn = get_fasta_path(args)

    'get Path to BAM file and generate images'
    if args.chr_id == 'all':
        for chr_i in range(1, 23):
            bam_fn = get_bam_path(args, chr_i)
            'get chromosome information'
            args.rd_count, args.chr_len = get_cov_info(bam_fn, chr_i)
            #args.rd_count, args.chr_len = 1128005051, 249250621
            'generate variant training images'
            vcf_var = prepare_variants(args, vcf_fn, chr_i)
            generate_training_img(vcf_var, bam_fn)

            'generate non variant training images'
            vcf_non_var = prepare_non_variants(args, fast_fn, chr_i)
            generate_training_img(vcf_non_var, bam_fn)
    else:
        bam_fn = get_bam_path(args, args.chr_id)
        'get chromosome information'
        args.rd_count, args.chr_len = get_cov_info(bam_fn, args.chr_id)
        #args.rd_count, args.chr_len = 1128005051, 249250621
        'generate variant training images'
        vcf_var = prepare_variants(args, vcf_fn, args.chr_id)
        generate_training_img(vcf_var, bam_fn)

        'generate non variant training images'
        vcf_non_var = prepare_non_variants(args, fast_fn, args.chr_id)
        #print("vcf generation done")
        generate_training_img(vcf_non_var, bam_fn)
