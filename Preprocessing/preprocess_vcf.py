import argparse
import subprocess
import shlex
import os
import vcf
import pandas as pd


def get_path(args):
    if args.bwcluster:
        path = os.path.join(args.work_dir_bwcluster, args.ws_dir)
    else:
        path = args.data_dir

    return path


def preprocess_vcf(args, sample_id, stat):
    if args.prepro_type == 'split':
        work_dir = get_path(args)
        vcf_all_samples_path = os.path.join(work_dir, args.vcf_fn)

        save_dir_vcf = os.path.join(work_dir, sample_id)
        if not os.path.exists(save_dir_vcf):
            os.makedirs(save_dir_vcf)
        name_vcf_out = os.path.join(save_dir_vcf, sample_id + "_all_sv.vcf.gz")

        p3 = subprocess.call(shlex.split("bcftools view -Oz -s %s -o %s %s" % (sample_id, name_vcf_out,  vcf_all_samples_path )))

    elif args.prepro_type == 'tabix':
        #vcf_one_sample_path = os.path.join(save_dir_vcf, sample_id + "_all_sv.vcf.gz")
        p2 = subprocess.call(shlex.split("tabix -p vcf %s " % (args.vcf_fn)))

    elif args.prepro_type == 'format':
        # p4 = subprocess.call(shlex.split("bcftools stats -s %s %s" % (sample_id, args.vcf_fn)))
        prepare_vcf(args, sample_id, stat)



def prepare_vcf(args, sample_id, stat):
    vcf_reader = vcf.Reader(open(args.vcf_fn))

    sample = vcf_reader.samples[0] if vcf_reader.samples else None

    l2 = list(map(str, range(1, 23, 1)))
    cnv_list = []

    for vcf_record in vcf_reader:
        # Ignore SV in chromosome X
        if vcf_record.CHROM not in l2: continue

        #check quality of variant
        if vcf_record.FILTER and "PASS" not in vcf_record.FILTER:continue

        # No change (0)
        if vcf_record.is_indel: continue

        # Ignore MNPs (0)
        if len(vcf_record.REF) != 1 and len(vcf_record.ALT[0]) != 1: continue

        # Ignore imprecise SVs (around 365)
        if "IMPRECISE" in vcf_record.INFO: continue

        #Ignoore triallelic
        #if len(vcf_record.ALT) > 1: continue

        # Filter for CNVs and get
        sv_type = vcf_record.INFO["SVTYPE"]

        'handle variants with CNV tag'
        if sv_type == "CNV":
            # Bi-allelic dublications --> include in Set
            if vcf_record.INFO['CS'] == 'DUP_gs':
                sv_type = "DUP"
            # deletions, dublications and other --> exclude (only 5-8% of all CNVs)
            elif vcf_record.INFO['CS'] == 'DUP_uwash':
                continue

        #only focus on CNVs
        if sv_type == "DEL" or sv_type =="DUP" or sv_type =="CNV":
            #check if variant is present in sample
            if vcf_record.genotype(args.sample_id).is_variant:

                'assign correct end position for dublication'
                if sv_type == "DUP":
                    cnv_r_pos = vcf_record.INFO["END"]
                    cnv_r_pos_sup = vcf_record.INFO["END"]
                    cnv_len = cnv_r_pos - vcf_record.POS

                elif sv_type == "DEL":
                    cnv_r_pos = vcf_record.INFO["END"]
                    cnv_r_pos_sup = vcf_record.POS + len(vcf_record.REF) - 1
                    cnv_len = cnv_r_pos - vcf_record.POS

                    #if cnv_len >= 10000:
                    #   continue

            else:continue
        else:continue
        cnv_list.append([vcf_record.CHROM, vcf_record.POS, cnv_r_pos, sv_type, cnv_len, sample])

    if stat:
        print("Statistic of :" + str(sample_id))
        get_cnv_statistics(cnv_list)

    'store the cnv_list in respective folder'
    cnv_df = pd.DataFrame(cnv_list)

    cnv_list_file = os.path.join(sample_id, sample_id + "_cnv_list.csv")
    cnv_df.to_csv(cnv_list_file, index=False, sep='\t', header=False)


def get_cnv_statistics(cnv_list):
    l3 = list(map(str, range(1, 23, 1)))

    'Distribution of CNV length'
    del_50, del_50_200, del_200_500, del_500_1000, del_1000_5000, del_5000_10000, del_10000 = 0, 0, 0, 0, 0, 0, 0
    dub_50, dub_50_200, dub_200_500, dub_500_1000, dub_1000_5000, dub_5000_10000, dub_10000 = 0, 0, 0, 0, 0, 0, 0

    for row in cnv_list:
        if row[0] in l3:
            l = int(row[2]) - int(row[1])
            if l <50:
                if row[3] == 'DEL': del_50 =+ 1
                elif row[3] == 'DUP': dub_50 += 1
            elif 50 >= l < 200:
                if row[3] == 'DEL': del_50_200 += 1
                elif row[3] == 'DUP': dub_50_200 += 1
            elif 200 >= l < 500:
                if row[3]== 'DEL': del_200_500 += 1
                elif row[3] == 'DUP': dub_200_500 += 1
            elif 500 >= l < 1000:
                if row[3] == 'DEL': del_500_1000 += 1
                elif row[3] == 'DUP': dub_500_1000 += 1
            elif 1000 >= l < 5000:
                if row[3] == 'DEL': del_1000_5000 += 1
                elif row[3] == 'DUP': dub_1000_5000 += 1
            elif 5000 >= l < 10000:
                if row[3] == 'DEL': del_5000_10000 += 1
                elif row[3] == 'DUP': dub_5000_10000 += 1
            else:
                if row[3] == 'DEL': del_10000 += 1
                elif row[3] == 'DUP': dub_10000 += 1

    del_total = del_50 + del_50_200 + del_200_500 + del_500_1000 + del_1000_5000 + del_5000_10000
    dub_total = dub_50 + dub_50_200 + dub_200_500 + dub_500_1000 + dub_1000_5000 + dub_5000_10000 + dub_10000

    print("Statistics of Deletions")
    print("DEL smaller than 50: " + str(del_50) + " / in percent: " + str(del_50/del_total if del_total != 0 else 0))
    print("DEL between 50 and 200: " + str(del_50_200) + " / in percent: " + str((del_50_200/del_total)*100 if del_total != 0 else 0)+"%")
    print("DEL between 200 and 500: " + str(del_200_500) + " / in percent: " + str((del_200_500/del_total)*100 if del_total != 0 else 0)+"%")
    print("DEL between 500 and 1000: " + str(del_500_1000) + " / in percent: " + str((del_500_1000/del_total)*100 if del_total != 0 else 0)+"%")
    print("DEL between 1000 and 5000: " + str(del_1000_5000) + " / in percent: " + str((del_1000_5000/del_total)*100 if del_total != 0 else 0)+"%")
    print("DEL between 5000 and 10000: " + str(del_5000_10000) + " / in percent: " + str((del_5000_10000/del_total)*100 if del_total != 0 else 0)+"%")
    print("DEL larger than 10000: " + str(del_10000))
    print("DEL in total: " + str(del_total))

    print("Statistics of Dublications")
    print("DUB smaller than 50: " + str(dub_50) + " / in percent: " + str(dub_50 / dub_total*100 if dub_total != 0 else 0)+ "%")
    print("DUB between 50 and 200: " + str(dub_50_200) + " / in percent: " + str((dub_50_200 / dub_total) * 100 if dub_total != 0 else 0) + "%")
    print("DUB between 200 and 500: " + str(dub_200_500) + " / in percent: " + str((dub_200_500 / dub_total) * 100 if dub_total != 0 else 0) + "%")
    print("DUB between 500 and 1000: " + str(dub_500_1000) + " / in percent: " + str((dub_500_1000 / dub_total) * 100 if dub_total != 0 else 0) + "%")
    print("DUB between 1000 and 5000: " + str(dub_1000_5000) + " / in percent: " + str((dub_1000_5000 / dub_total) * 100 if dub_total != 0 else 0) + "%")
    print("DUB between 5000 and 10000: " + str(dub_5000_10000) + " / in percent: " + str((dub_5000_10000 / dub_total) * 100 if dub_total != 0 else 0) + "%")
    print("DUB larger than 10000: " + str(dub_10000))
    print("DUB in total: " + str(dub_total))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Preprocess VCF Files")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools'")

    parser.add_argument('--bwcluster', type=bool, default=True,
                        help="Define if script is run on bwcluster or not")
    parser.add_argument('--work_dir_bwcluster', type=str, default="/pfs/work7/workspace/scratch",
                        help="Path of work space directory on bwcluster")
    parser.add_argument('--ws_dir', type=str, default="db0052-ma_deepsv-0/data",
                        help="Path of workspace directory")
    parser.add_argument('--work_dir', type=str, default="/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/Preprocessing",
                        help="Path of working directory")
    parser.add_argument('--data_dir', type=str, default="/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/Preprocessing",
                        help="Path of data directory")

    parser.add_argument('--sample_id', type=str, default="NA20502",
                        help="ID of sample to be processed")
    parser.add_argument('--sample_all', type=bool, default=True,
                        help="Preprocess vcf for all samples")
    parser.add_argument('--with_stat', type=bool, default=True,
                        help="Show statistics of CNVs")
    parser.add_argument('--prepro_type', type=str, default='format',
                        help="Preprocess types: split, tabix, format")
    parser.add_argument('--vcf_fn', type=str, default='ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz',
                        help="VCF file to be processed")

    args = parser.parse_args()

    'Prepare path to vcf file for subset preprocessing'
    #if args.prepro_type == 'split' and args.vcf_fn is not None:
    #    args.vcf_fn = os.path.join(args.work_dir, "VCF", args.vcf_fn)

    'Get sample id and split vcf by sample'
    l1 = ['NA18525', 'NA19625', 'NA19648', 'NA20502']
    cnv_list = []

    if args.sample_all:
        for sample in l1:
            args.work_dir = get_path(args)
            if args.prepro_type == 'tabix':
                args.vcf_fn = os.path.join(args.work_dir, sample, sample + "_all_sv.vcf")
            elif args.prepro_type == 'format':
                args.sample_id = sample
                args.vcf_fn = os.path.join(args.work_dir, args.sample_id, args.sample_id + "_all_sv.vcf")
            preprocess_vcf(args, sample, args.with_stat)
    elif args.sample_all == False and args.sample_id is not None:
        args.work_dir = get_path(args)
        if args.prepro_type == 'tabix':
            args.vcf_fn = os.path.join(args.work_dir, args.sample_id, args.sample_id + "_all_sv.vcf.gz")
        elif args.prepro_type == 'format':
            args.vcf_fn = os.path.join(args.work_dir, args.sample_id, args.sample_id + "_all_sv.vcf")
        preprocess_vcf(args, args.sample_id, args.with_stat)



