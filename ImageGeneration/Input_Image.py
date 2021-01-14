from PIL import Image, ImageDraw
import os
import numpy as np
from scipy.stats import poisson


# returns start and end gene locus for each training image
def gene_point_pic(chr_id, pos_l, pos_r):
    gene_pic = []
    every_len = pos_l
    while every_len < pos_r:
        gene_tuple = (chr_id, every_len, every_len + 50)
        gene_pic.append(gene_tuple)
        every_len = every_len + 50
    if every_len >= pos_r:
        gene_tuple = (chr_id, every_len - 50, pos_r)
        gene_pic.append(gene_tuple)
    return gene_pic


# generates pileup sequence
def pipeup_column(sam_file, chr_id, pos_l, pos_r):
    pile_record = []
    for pileupcolumn in sam_file.pileup(chr_id, pos_l, pos_r):
        if pos_l <= pileupcolumn.pos <= pos_r:
            for pileupread in pileupcolumn.pileups:
                if None != pileupread.alignment.cigarstring and None != pileupread.query_position:
                    index = 0
                    map_type = -1
                    for cigar in pileupread.alignment.cigartuples:
                        if pileupread.query_position > index:
                            index = cigar[1] + index
                            map_type = cigar[0]

                    pile_result = (
                        pileupcolumn.pos, pileupread.alignment.is_paired, pileupread.alignment.is_proper_pair,
                        pileupread.alignment.mapping_quality, map_type,
                        pileupread.alignment.query_sequence[pileupread.query_position])
                    pile_record.append(pile_result)
    return pile_record


# divide image
'''TODO: Add line with information about stats'''
#https://stackoverflow.com/questions/13053443/drawing-a-line-on-an-image-with-pil
def draw_pic(clip_dict_record, pile_record, del_pos_np_start, deletion_length, args, chr, svtype):
    blank = Image.new("RGB", [256, 256], "white")
    pile_record_len = len(pile_record)
    drawObject = ImageDraw.Draw(blank)
    y_start_index = 0
    old_x_start = 5
    for j in range(pile_record_len):
        #print("-- %d" % pile_record[j][0])
        #print("--- %d" % del_pos_np_start)
        x_start = (pile_record[j][0] - del_pos_np_start) * 5 + 5
        #print("x_start %d " % x_start)
        if old_x_start == x_start:
            old_x_start = x_start
            #print("x_pic_start %d" % x_start)
            y_start = 5 + y_start_index * 5
            #print("y_pic_start %d" % y_start)
            #print("y_index %d" % y_start_index)
            y_start_index += 1
            x_end = x_start + 5
            y_end = y_start + 5
            if pile_record[j][0] in clip_dict_record:
                base_rgb = get_rgb(-clip_dict_record[pile_record[j][0]], pile_record[j])
            else:
                base_rgb = get_rgb(0, pile_record[j])
            #print("rgb")
            #print(base_rgb)
            drawObject.rectangle((x_start, y_start, x_end, y_end), fill=base_rgb)
        elif old_x_start != x_start:
            old_x_start = x_start
            #print("x_pic_start %d" % x_start)
            y_start_index = 0
            y_start = 5 + y_start_index * 5
            #print("y_pic_start %d" % y_start)
            y_start_index += 1
            x_end = x_start + 5
            y_end = y_start + 5

            if pile_record[j][0] in clip_dict_record:
                base_rgb = get_rgb(-clip_dict_record[pile_record[j][0]], pile_record[j])
            else:
                base_rgb = get_rgb(0, pile_record[j])
            #print("rgb")
            #print(base_rgb)
            drawObject.rectangle((x_start, y_start, x_end, y_end), fill=base_rgb)

    med, lower_percentile, upper_percentile = get_chr_stats(args.rd_count, args.chr_len, 50)
    drawObject.line((0, med, drawObject.size[0], med), fill=(0,0,0))


    'Save images with label'
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        bam_fn_loc = os.path.join(bw_path, "chr_" + str(chr))
        img_path = os.path.join(bam_fn_loc, "img_deepsv")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        img_path = os.path.join(local_path,"chr_" + str(chr), "img_deepsv")

    if svtype == "DEL":
        img_name = 'del_1_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if svtype == "DUP":
        img_name = 'dup_2_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if svtype == "NOCNV":
        img_name = 'no_0_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if not os.path.exists(img_path):
        os.makedirs(img_path)
    blank.save(str(img_path) + "/" + str(img_name), "PNG")


#get correct color for image
def get_rgb(clip_value, every_pile_record):
    if every_pile_record[5] == 'A':

        # Red
        base_A = [255, 0, 0]
        if every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 0, 0]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 60, 60]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 70, 70]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 80, 80]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 90, 90]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 100, 100]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 110, 110]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 120, 120]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 130, 130]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 140, 140]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 150, 150]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 160, 160]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 170, 170]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 180, 180]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_A = [255, 190, 190]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_A = [255, 200, 200]

        base_A[1] = base_A[1] + clip_value
        base_A[2] = base_A[2] + clip_value
        base_A = tuple(base_A)
        return base_A

    elif every_pile_record[5] == 'T':
        # green
        base_T = [0, 255, 0]
        if every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_T = [0, 255, 0]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_T = [60, 255, 60]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_T = [70, 255, 70]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_T = [80, 255, 80]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_T = [90, 255, 90]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_T = [100, 255, 100]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_T = [110, 255, 110]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_T = [120, 255, 120]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_T = [130, 255, 130]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_T = [140, 255, 140]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_T = [150, 255, 150]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_T = [160, 255, 160]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_T = [170, 255, 170]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_T = [180, 255, 180]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_T = [190, 255, 190]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_T = [200, 255, 200]

        base_T[0] = base_T[0] + clip_value
        base_T[2] = base_T[2] + clip_value
        base_T = tuple(base_T)
        return base_T

    elif every_pile_record[5] == 'C':
        # blue
        base_C = [0, 0, 255]
        if every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_C = [0, 0, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_C = [60, 60, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_C = [70, 70, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_C = [80, 80, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_C = [90, 90, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_C = [100, 100, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_C = [110, 110, 255]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_C = [120, 120, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_C = [130, 130, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_C = [140, 140, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_C = [150, 150, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_C = [160, 160, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_C = [170, 170, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_C = [180, 180, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_C = [190, 190, 255]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_C = [200, 200, 255]

        base_C[0] = base_C[0] + clip_value
        base_C[1] = base_C[1] + clip_value
        base_C = tuple(base_C)
        return base_C

    elif every_pile_record[5] == 'G':
        # black
        base_G = [0, 0, 0]
        if every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_G = [0, 0, 0]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_G = [60, 60, 60]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_G = [70, 70, 70]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_G = [80, 80, 80]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_G = [90, 90, 90]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_G = [100, 100, 100]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_G = [110, 110, 110]
        elif every_pile_record[1] == 'True' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_G = [120, 120, 120]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_G = [130, 130, 130]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_G = [140, 140, 140]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_G = [150, 150, 150]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'True' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_G = [160, 160, 160]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] != 4:
            base_G = [170, 170, 170]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] >= 20 and \
                every_pile_record[4] == 4:
            base_G = [180, 180, 180]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] != 4:
            base_G = [190, 190, 190]
        elif every_pile_record[1] == 'False' and every_pile_record[2] == 'False' and every_pile_record[3] < 20 and \
                every_pile_record[4] == 4:
            base_G = [200, 200, 200]

        base_G = np.array(base_G)
        base_G = base_G + clip_value
        base_G = tuple(base_G)
        return base_G


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
