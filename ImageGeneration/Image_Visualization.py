import numpy as np
import pysam
import subprocess
import shlex
import os
from PIL import Image, ImageDraw
from scipy.stats import poisson


def visual_training_img(args, vcf_tuples, bam_path):
    """ loop through tuples of CNVs/NoCNVs (chr_id, l_pos, r_pos, sv_type) """
    num_cnvs = len(vcf_tuples)
    cnv_np_ar = np.array(vcf_tuples)
    add_on_side = args.img_size * 0.6

    for i in range(num_cnvs):
        scan_l_pos = int(cnv_np_ar[i, 1]) - add_on_side
        scan_r_pos = int(cnv_np_ar[i, 2]) + add_on_side
        svLen = scan_r_pos - scan_l_pos
        svtype = cnv_np_ar[i, 3]
        chrom = cnv_np_ar[i, 0]

        #print("------")
        #print("New Variant: " + str(svtype) + ": " + str(scan_l_pos) + " - " + str(scan_r_pos))

        """ divide region into 50 bp long bins """
        gene_pic = gene_point_pic(chrom, scan_l_pos, scan_r_pos, args.img_size)

        """ get list of discordant read pairs """
        disc_list = get_disc_read(bam_path, chrom, scan_l_pos,scan_r_pos)

        #pile_record = pileup_column(bam_path, chrom, scan_l_pos,scan_r_pos, disc_list)
        #deletion_length = scan_r_pos - scan_l_pos
        #draw_pic(pile_record, scan_l_pos, deletion_length, args, chr, svtype)

        for each_pic in gene_pic:
           pile_record = pileup_column(bam_path, each_pic[0], each_pic[1], each_pic[2], disc_list)

           if svtype == "NOCNV":
                if len(pile_record) < 1:
                    print(each_pic[1])
                    continue
           deletion_length = each_pic[2] - each_pic[1]
           draw_pic(pile_record, each_pic[1], deletion_length, args, chrom, svtype, args.img_size)

    return "hello"


# returns start and end gene locus for each training image
def gene_point_pic(chr_id, pos_l, pos_r, img_size):
    gene_pic = []
    num_img = (pos_r - pos_l) // int(img_size)
    every_len = pos_l

    for i in range(int(num_img) + 1):
        if i < num_img:
            gene_tuple = (chr_id, every_len, every_len + img_size)
            gene_pic.append(gene_tuple)
            every_len = every_len + img_size

        elif i == num_img:
            gene_tuple = (chr_id, pos_r - img_size, pos_r)
            gene_pic.append(gene_tuple)

    return gene_pic


# returns list with discordant read pairs
def get_disc_read(bam_path, chr_id,  scan_l_pos, scan_r_pos):
    p2 = subprocess.Popen(shlex.split("samtools view -F 1294 %s %s:%d-%d" % (bam_path, chr_id, scan_l_pos, scan_r_pos)),
                          stdout=subprocess.PIPE)
    disc_list = []
    for l in p2.stdout:
        l = l.split()
        if l[0][0] == "@":
            continue
        disc_list.append(l[0].decode('utf-8'))

    return disc_list


# generates pileup sequence
def pileup_column(bam_path, chr_id, pos_l, pos_r, discordant_reads):
    sam_file = pysam.AlignmentFile(bam_path, "rb")
    pile_record = []
    standard = 0
    discordant = 0
    split = 0

    for pileupcolumn in sam_file.pileup(chr_id, pos_l, pos_r, stepper="all", flag_filter=8):
        if pos_l <= pileupcolumn.pos <= pos_r:
            for pileupread in pileupcolumn.pileups:
                is_clipped = None
                hard_clipped = None
                if None != pileupread.alignment.cigarstring and None != pileupread.query_position:
                    pileupread_re = rearrange_string(pileupread.alignment)

                    if pileupread.alignment.is_paired:  # read not paired，do not consider this kind of read
                        if 'S' in pileupread.alignment.cigarstring:
                            is_clipped = True
                        else:
                            is_clipped = False
                        if 'H' in pileupread.alignment.cigarstring:
                            hard_clipped = True
                        else:
                            hard_clipped = False
                    else:
                        continue

                    pile_result = (pileupcolumn.pos, pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position],
                                   pileupread_re[pileupread.query_position], pileupread.alignment.is_paired, pileupread.alignment.is_proper_pair,
                                   pileupread.alignment.is_read1, pileupread.alignment.mapping_quality, is_clipped, hard_clipped)

                    #count reads
                    if pile_result[5] is False:
                        discordant = discordant + 1
                    if pile_result[8] is True:
                        split = split + 1
                    if pile_result[5] is True and pile_result[8] is False:
                        standard = standard + 1

                    pile_record.append(pile_result)

    # print(standard)
    # print(discordant)
    # print(split)

    return pile_record


# adds cigarstring information to base sequence of read
def rearrange_string(read):
    bases = read.query_sequence
    bases_checked = 0
    read_len_count = 0
    new_bases = ''

    for cigar_portion in read.cigartuples:
        # match
        if cigar_portion[0] == 0:
            cigar_base = bases[bases_checked:(cigar_portion[1] + bases_checked)]
            new_bases = new_bases + cigar_base
            for M_num in range(cigar_portion[1]):
                read_len_count = read_len_count + 1
            bases_checked = bases_checked + cigar_portion[1]

        # insertion
        elif cigar_portion[0] == 1:
            cigar_base = ''
            bases_checked = bases_checked + cigar_portion[1]
            for I_num in range(cigar_portion[1]):
                read_len_count = read_len_count + 1
                cigar_base = cigar_base + 'i'
            new_bases = new_bases + cigar_base

        # deletion
        elif cigar_portion[0] == 2:
            cigar_base = ''
            for i in range(cigar_portion[1]):
                cigar_base = cigar_base + 'd'
            new_bases = new_bases + cigar_base

        # soft-clipped
        elif cigar_portion[0] == 4:
            cigar_base = ''
            for i in range(cigar_portion[1]):
                cigar_base = cigar_base + 's'
                read_len_count = read_len_count + 1
            new_bases = new_bases + cigar_base
            bases_checked = bases_checked + cigar_portion[1]

        # hard-clipped
        elif cigar_portion[0] == 5:
            cigar_base = ''
            for i in range(cigar_portion[1]):
                cigar_base = cigar_base + 'h'
            new_bases = new_bases + cigar_base

    return new_bases


# function to draw training image
def draw_pic(pile_record, del_pos_np_start, deletion_length, args, chr_id, svtype, img_size):
    # draw empty image
    size = img_size * 5
    size = size + 15
    blank = Image.new("RGB", [size, size], "white")
    drawObject = ImageDraw.Draw(blank)

    # get values
    pile_record_len = len(pile_record) #number of reads within 50bp regions
    y_start_index = 0
    old_x_start = 5

    # loop over every entry in pile_record, get rectangle pos and rgb value
    for j in range(pile_record_len):
        x_start = (pile_record[j][0] - del_pos_np_start) * 5 + 5
        if old_x_start == x_start:
            old_x_start = x_start
            y_start = 5 + y_start_index * 5
            y_start_index += 1
            x_end = x_start + 5
            y_end = y_start + 5
            #print("rectangle pos: " + str(x_start) + " " + str(x_end) + " / " + str(y_start) + " " + str(y_end))
            base_rgb = get_rgb_value(pile_record[j])
            if base_rgb == None:
                print(pile_record[j])
            drawObject.rectangle((x_start, y_start, x_end, y_end), fill=tuple(base_rgb))

        elif old_x_start != x_start:
            old_x_start = x_start
            y_start_index = 0
            y_start = 5 + y_start_index * 5
            y_start_index += 1
            x_end = x_start + 5
            y_end = y_start + 5

            #print("rectangle pos: " + str(x_start) + " " + str(x_end) + " / " + str(y_start) + " " + str(y_end))
            base_rgb = get_rgb_value(pile_record[j])
            if base_rgb == None:
                print(pile_record[j])
            drawObject.rectangle((x_start, y_start, x_end, y_end), fill=tuple(base_rgb))

    med, lower_percentile, upper_percentile = get_chr_stats(args.rd_count, args.chr_len, 50)
    drawObject.line((0, med*5+5, size, med*5+5), fill=(0,0,0))

    'Save images with label'
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        bam_fn_loc = os.path.join(bw_path, "chr_" + str(chr_id))
        img_path = os.path.join(bam_fn_loc, "img_vis_50")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        img_path = os.path.join(local_path,"chr_" + str(chr_id), "img_test3")

    if svtype == "DEL":
        img_name = 'del_1_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if svtype == "DUP":
        img_name = 'dup_2_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if svtype == "NOCNV":
        img_name = 'no_0_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if not os.path.exists(img_path):
        os.makedirs(img_path)
    final_img_path = os.path.join(img_path, img_name)
    blank.save(final_img_path, "PNG")


#define coler for pile_record
def get_rgb_value(pile_record):
    base = pile_record[2]
    base_re = pile_record[3]
    is_paired = pile_record[4]
    is_proper_pair = pile_record[5]
    is_read1 = pile_record[6]
    map_quality = pile_record[7]
    is_clipped = pile_record[8]
    #hard_clipped = pile_record[9] No need, in bwa hard clipping is used for supplementary reads

    # handle discordant reads
    if not is_proper_pair:
        if is_read1:
            if base_re == 's' or base_re == 'd':
                return [255, 255, 255]
            else:
                # lightblue corrected with map quality
                add_qual = get_map_qual_color(map_quality)
                red = 0 + add_qual
                green = 255
                blue = 255
                return [red, green, blue]

        else:
            if base_re == 's' or base_re == 'd':
                return [255, 255, 255]
            else:
                # pink corrected with map quality
                add_qual = get_map_qual_color(map_quality)
                red = 255
                green = 0 + add_qual
                blue = 255
                return [red, green, blue]

    # handle split reads
    # True if cigarstring contains S
    elif is_clipped:
        if base_re == "s":
            # yellow corrected with map quality
            add_qual = get_map_qual_color(map_quality)
            red = 255
            green = 255
            blue = 0 + add_qual
            return [red, green, blue]
        else:
            bas_col = get_base_color(base, map_quality)
            return bas_col

    # handle all other reads
    elif is_paired and is_proper_pair:
        if map_quality >= 1:
            bas_col = get_base_color(base, map_quality)
            return bas_col
        elif map_quality == 0 or map_quality == -1:
            bas_col = get_base_color(base, map_quality)
            return bas_col

    else:
        return [255, 255, 255]


def get_base_color(pile_record_base, map_qual):
    if pile_record_base == "A" or pile_record_base == "a":
        # red corrected with map quality
        add_qual = get_map_qual_color(map_qual)
        red = 255
        green = 0 + add_qual
        blue = 0 + add_qual
        return [red, green, blue]

    elif pile_record_base == "C" or pile_record_base == "c":
        # blue corrected with map quality
        add_qual = get_map_qual_color(map_qual)
        red = 0 + add_qual
        green = 0 + add_qual
        blue = 255
        return [red, green, blue]

    elif pile_record_base == "T" or pile_record_base == "t":
        # green corrected with map quality
        add_qual = get_map_qual_color(map_qual)
        red = 0 + add_qual
        green = 255
        blue = 0 + add_qual
        return [red, green, blue]

    elif pile_record_base == "G" or pile_record_base == "g":
        # black corrected with map quality
        add_qual = get_map_qual_color(map_qual)
        red = 0 + add_qual
        green = 0 + add_qual
        blue = 0 + add_qual
        return [red, green, blue]

    else:
        print(pile_record_base)
        # white
        return [255, 255, 255]


def get_map_qual_color(map_qual):
    if map_qual >= 30:
        return 0
    elif 20 <= map_qual <= 29:
        return 60
    elif 10 <= map_qual <= 19:
        return 120
    elif 1 <= map_qual <= 9:
        return 180
    elif map_qual <= 0:
        return 240


def get_chr_stats(dep_count, gen_len, bin_size):
    num_rd_per_base = dep_count
    bin_size = bin_size
    len_genome = gen_len
    d = 3

    mu = (num_rd_per_base / len_genome) * bin_size
    var = mu / (d - 1)

    med = poisson.median(mu) / bin_size
    l_p, u_p = poisson.interval(0.9, mu)
    lower_percentile = l_p / bin_size
    upper_percentile = u_p / bin_size

    return med, lower_percentile, upper_percentile







