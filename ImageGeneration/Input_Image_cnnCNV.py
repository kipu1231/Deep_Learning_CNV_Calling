# coding:utf-8
import pysam
import os
from PIL import Image
import numpy as np
from PIL import Image, ImageDraw
from scipy.stats import poisson

'''TODO: Add line with information about stats'''
#https://stackoverflow.com/questions/13053443/drawing-a-line-on-an-image-with-pil
def generate_cnn_cnv_images(args, vcf_tuples, bam_path):
    num_cnvs = len(vcf_tuples)
    cnv_np_ar = np.array(vcf_tuples)

    for i in range(num_cnvs):
        scan_l_pos = int(cnv_np_ar[i, 1]) - 200
        scan_r_pos = int(cnv_np_ar[i, 2]) + 200
        svLen = scan_r_pos - scan_l_pos
        svtype = cnv_np_ar[i, 3]
        chrom = cnv_np_ar[i, 0]

        'divide region into 50 bp long bins'
        gene_pic = gene_point_pic_cnncnv(chrom, scan_l_pos, scan_r_pos)
        for each_pic in gene_pic:
            pile_record = pipeup_column_cnncnv(bam_path, each_pic[0], each_pic[1], each_pic[2])
            del_len = each_pic[2] - each_pic [1]
            draw_pic(pile_record, each_pic[1], del_len, args, chr, svtype)

        # popen_bam = os.popen(
        #     'samtools depth -r ' + chrom + ':' + str(scan_l_pos) + '-' + str(scan_r_pos) + ' ' + bam_path)
        # sys_line = popen_bam.readlines()
        # popen_bam.close()
        # depth_lines = []
        # for dep in sys_line:
        #     depth_lines.append(int(dep.strip('\n').split('\t')[2]))
        #
        # if len(depth_lines) == 0:
        #     height = 1
        # else:
        #     height = int(max(depth_lines)) * 2
        #
        # width = scan_r_pos - scan_l_pos
        #
        # 'save image in correct folder'
        # if args.bwcluster:
        #     bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        #     bam_fn_loc = os.path.join(bw_path, "chr_" + str(chrom))
        #     img_path = os.path.join(bam_fn_loc, "img_cnncnv")
        # else:
        #     local_path = os.path.join(args.data_dir, args.sample_id)
        #     img_path = os.path.join(local_path, "chr_" + str(chrom), "img_cnncnv")
        #
        # if not os.path.exists(img_path):
        #     os.makedirs(img_path)
        #
        # if svtype == "DEL":
        #     img_name = 'del_1_' + str(scan_l_pos) + '_' + str(scan_r_pos) + ".png"
        #
        # if svtype == "DUP":
        #     img_name = 'dup_2_' + str(scan_l_pos) + '_' + str(scan_r_pos) + ".png"
        #
        # if svtype == "NOCNV":
        #     img_name = 'no_0_' + str(scan_l_pos) + '_' + str(scan_r_pos) + ".png"
        #
        # save_img_name = os.path.join(img_path, img_name)
        #
        # read_package = []
        # c = 0
        # bam_file = pysam.AlignmentFile(bam_path, "rb")
        # if scan_l_pos < 0:
        #     scan_l_pos = 1
        #     scan_r_pos = scan_r_pos - scan_l_pos + 1
        # for read in bam_file.fetch(chrom, scan_l_pos, scan_r_pos):
        #     if (read.cigarstring is not None) and read_can_shown(read, scan_l_pos,
        #                                                          scan_r_pos) and read.mapping_quality > 10:
        #         new_base_string = rearrange_string(read)
        #         if read_corner_shown(read, scan_l_pos, scan_r_pos, new_base_string):
        #             read_package.append((new_base_string, read))
        #             c = c + 1
        #
        # dic, height_new = read_to_dictionary(read_package, scan_r_pos, height)
        # new_dic = rearrange_read_dictionary(dic)
        # print('ok')
        # if not c == 0:
        #     print(c)
        #     draw_pgn('left', dic, width, height_new, scan_l_pos, scan_r_pos, save_img_name)
        # else:
        #     print("empty")
        #     draw_empty('left', dic, 100, 100, scan_l_pos, scan_r_pos, save_img_name)


# returns start and end gene locus for each training image
def gene_point_pic_cnncnv(chr_id, pos_l, pos_r):
    gene_pic = []
    every_len = pos_l
    while every_len < pos_r:
        gene_tuple = (chr_id, every_len, every_len + 100)
        gene_pic.append(gene_tuple)
        every_len = every_len + 100
    if every_len >= pos_r:
        gene_tuple = (chr_id, every_len - 100, pos_r)
        gene_pic.append(gene_tuple)
    return gene_pic


# generates pileup sequence
def pipeup_column_cnncnv(bam_path, chr_id, pos_l, pos_r):
    sam_file = pysam.AlignmentFile(bam_path, "rb")
    pile_record = []

    for pileupcolumn in sam_file.pileup(chr_id, pos_l, pos_r):
        if pos_l <= pileupcolumn.pos <= pos_r:
            for pileupread in pileupcolumn.pileups:
                is_clipped = None
                hard_clipped = None
                if None != pileupread.alignment.cigarstring and None != pileupread.query_position:
                    index = 0
                    map_type = -1
                    pileupread_re = rearrange_string(pileupread.alignment)

                    for cigar in pileupread.alignment.cigartuples:
                        if pileupread.query_position > index:
                            index = cigar[1] + index
                            map_type = cigar[0]

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

                    # print('----')
                    # print(pileupread.alignment.query_length)
                    # print(len(pileupread_re))
                    # print(pileupread.alignment.cigarstring)
                    # print(pileupread.alignment.query_sequence)
                    # print(pileupread_re)


                    pile_result = (
                        pileupcolumn.pos, pileupread.alignment.is_paired, pileupread.alignment.is_proper_pair,
                        pileupread.alignment.is_read1, is_clipped, hard_clipped, pileupread.alignment.mapping_quality,
                        map_type,
                        pileupread.alignment.query_sequence[pileupread.query_position], pileupread_re[pileupread.query_position])
                    pile_record.append(pile_result)
    return pile_record


#adds cigarstring information to base sequence of read
def rearrange_string(read):
    bases = read.query_sequence
    bases_checked = 0
    read_len_count = 0
    new_bases = ''

    for cigar_portion in read.cigartuples:
        if cigar_portion[0] == 0:
            cigar_base = bases[bases_checked:(cigar_portion[1] + bases_checked)]
            new_bases = new_bases + cigar_base
            for M_num in range(cigar_portion[1]):
                read_len_count = read_len_count + 1
            bases_checked = bases_checked + cigar_portion[1]

        elif cigar_portion[0] == 1:
            cigar_base = ''
            bases_checked = bases_checked + cigar_portion[1]
            for I_num in range(cigar_portion[1]):
                read_len_count = read_len_count + 1
                cigar_base = cigar_base + 'i'
            new_bases = new_bases + cigar_base

        elif cigar_portion[0] == 2:
            cigar_base = ''
            for i in range(cigar_portion[1]):
                cigar_base = cigar_base + 'd'
            new_bases = new_bases + cigar_base

        elif cigar_portion[0] == 4:
            cigar_base = ''
            for i in range(cigar_portion[1]):
                cigar_base = cigar_base + 's'
                read_len_count = read_len_count + 1
            new_bases = new_bases + cigar_base
            bases_checked = bases_checked + cigar_portion[1]

        elif cigar_portion[0] == 5:
            cigar_base = ''
            for i in range(cigar_portion[1]):
                cigar_base = cigar_base + 'h'
            new_bases = new_bases + cigar_base

    return new_bases


#draw image
def draw_pic(pile_record, del_pos_np_start, deletion_length, args, chr, svtype):
    blank = Image.new("RGB", [505, 505], "white")
    draw_object = ImageDraw.Draw(blank)
    y_start_index = 0
    old_x_start = 5
    h = 0
    l = 0

    for j in range(len(pile_record)):
        x_start = (pile_record[j][0] - del_pos_np_start) * 5 + 5

        'check if entry if at beginning'
        if old_x_start == x_start:
            old_x_start = x_start
            y_start = 5 + y_start_index * 5
            y_start_index += 1
            x_end = x_start + 5
            y_end = y_start + 5

            'get rgb value for pile record and draw image'
            red, green, blue = get_RGB(pile_record[j])
            if green == 255:
                h = h +1
            else:
                l = l +1

            draw_object.rectangle((x_start, y_start, x_end, y_end), fill=(red, green, blue))

        elif old_x_start != x_start:
            old_x_start = x_start
            y_start_index = 0
            y_start = 5 + y_start_index * 5
            y_start_index += 1
            x_end = x_start + 5
            y_end = y_start + 5

            'get rgb value for pile record and draw image'
            red, green, blue = get_RGB(pile_record[j])
            if green == 255:
                h = h +1
            else:
                l = l +1
            draw_object.rectangle((x_start, y_start, x_end, y_end), fill=(red, green, blue))

    'Save images with label'
    if args.bwcluster:
        bw_path = os.path.join(args.work_dir_bwcluster, args.ws_dir, args.sample_id)
        bam_fn_loc = os.path.join(bw_path, "chr_" + str(chr))
        img_path = os.path.join(bam_fn_loc, "img_cnncnv")
    else:
        local_path = os.path.join(args.data_dir, args.sample_id)
        img_path = os.path.join(local_path,"chr_" + str(chr), "img_cnncnv")

    if svtype == "DEL":
        img_name = 'del_1_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if svtype == "DUP":
        img_name = 'dup_2_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if svtype == "NOCNV":
        img_name = 'no_0_' + str(del_pos_np_start) + '_' + str(del_pos_np_start + deletion_length) + ".png"

    if not os.path.exists(img_path):
        os.makedirs(img_path)
    blank.save(str(img_path) + "/" + str(img_name), "PNG")

    print(h)
    print(l)


#change color according to mapping quality
def quality_to_vRGB(quality):
    if quality >= 30:
        return 0
    elif 20 <= quality <= 29:
        return 50
    elif 10 <= quality <= 19:
        return 100
    elif 0 <= quality <= 10:
        return 180


#define coler for pile_record
def get_RGB(pile_record):
    base = pile_record[9]
    quality = pile_record[6]
    is_clipped = pile_record[4]
    hard_clipped = pile_record[5]
    is_read1 = pile_record[3]
    is_paired = pile_record[1]
    is_proper_pair = pile_record[2]

    if is_paired and is_proper_pair:
        is_concordant = True
    else:
        is_concordant = False

    if is_clipped:
        if base in "atcgATCG":
            red = 0
            green = 255
            blue = 0
            red = red + quality_to_vRGB(quality)
            blue = blue + quality_to_vRGB(quality)
            return red, green, blue
        elif base == "s":
            print('soft clipped')
            red = 255
            green = 0
            blue = 0
            green = green + quality_to_vRGB(quality)
            blue = blue + quality_to_vRGB(quality)
            return red, green, blue

    elif hard_clipped:
        print('hard clipped')
        red = 255
        green = 0
        blue = 255
        if base == "h":
            return 255, 255, 255
        return red, green, blue

    elif is_concordant:
        red = 0
        green = 255
        blue = 0
        if quality > 1:
            red = red + quality_to_vRGB(quality)
            blue = blue + quality_to_vRGB(quality)
            return red, green, blue
        elif quality == 0 or quality == -1:
            return 255, 255, 255
    elif not is_concordant:
        red = 0
        green = 255
        blue = 0
        if is_read1:
            if quality > 1:
                print('discordant read')
                red = 255
                blue = blue + quality_to_vRGB(quality)
                return red, green, blue
            elif base == 's' or base == 'd':
                return 255, 255, 255
        else:
            if quality > 1:
                print('discordant read')
                blue = 255
                red = red + quality_to_vRGB(quality)
                return red, green, blue
            elif base == 's' or base == 'd':
                return 255, 255, 255
    return 255, 255, 255


def get_chr_stats_cnncnv(dep_count, gen_len, bin_size):
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


def read_can_shown(read, scan_l_pos, scan_r_pos):
    read_pos1 = read.reference_start
    read_pos2 = read.reference_start + read_inferred_len(read)
    if (read_pos2 > scan_l_pos) and (read_pos1 < scan_r_pos):
        res = True
        for cigar_portion in read.cigartuples:
            if not ((cigar_portion[0] == 0) or (cigar_portion[0] == 1) or (cigar_portion[0] == 2) or (
                    cigar_portion[0] == 4) or (cigar_portion[0] == 5)):
                res = False
        return res
    else:
        return False


def read_corner_shown(read, scan_l_pos, scan_r_pos, new_bases):
    read_pos1 = read.reference_start
    read_pos2 = read.reference_start + read_inferred_len(read)
    if (read_pos1 < scan_l_pos) and (read_pos2 > scan_l_pos):
        if ('A' or 'G' or 'C' or 'T' or 'a' or 'g' or 'c' or 't') in new_bases[(scan_l_pos - read_pos1):len(new_bases)]:
            return True
        else:
            return False

    if (read_pos1 < scan_r_pos) and (read_pos2 > scan_r_pos):
        if ('A' or 'G' or 'C' or 'T' or 'a' or 'g' or 'c' or 't') in new_bases[0:(scan_r_pos - read_pos1)]:
            return True
        else:
            return False

    if (read_pos1 >= scan_l_pos) and (read_pos2 <= scan_r_pos):
        return True


def read_inferred_len(read):
    infer_len = 0
    for cigar_portion in read.cigartuples:
        if (cigar_portion[0] == 0) or (cigar_portion[0] == 2) or (cigar_portion[0] == 4):
            infer_len = infer_len + cigar_portion[1]
    return infer_len


def is_empty(read_list):
    tag = True
    for li in read_list:
        if li:
            tag = False
            break
    return tag


def get_shortest_tail_row(read_list, scan_r_pos):
    if is_empty(read_list):
        return 0
    else:
        tail = scan_r_pos
        short_row = 0
        for i in range(len(read_list)):
            if read_list[i]:
                if tail > read_list[i][-1][1]:
                    tail = read_list[i][-1][1]
                    short_row = i
        return short_row


def get_nearest_tail_row(read_list, scan_r_pos, read_pos1, read_pos2):
    near_row = 0
    insert_id = 0
    out = 0
    for i in range(len(read_list)):
        if read_list[i] != [] and out == 0:
            for j in range(len(read_list[i])):

                if j > 0 and (read_pos1 > read_list[i][j - 1][1]) and (read_pos2 < read_list[i][j][0]):
                    near_row = i
                    insert_id = j
                    out = 1
                    break
                elif j == len(read_list[i]) - 1 and read_pos1 > read_list[i][-1][1]:
                    near_row = i
                    insert_id = len(read_list[i])

    if near_row == 0:
        near_row = find_next_empty_row(read_list)
        insert_id = 0
    return near_row, insert_id


def find_next_empty_row(read_list):
    row = 0
    for i in range(len(read_list)):
        if not read_list[i]:
            row = i
            break
    return row


def rearrange_read_dictionary(dictionary):
    soft_clip = []
    hard_clip = []
    read1 = []
    read2 = []
    normal = []

    new_dic = {}

    for line in dictionary.values():
        tag = 0
        for read in line:
            if read[4]:
                tag = 4
                break
            elif read[4] == False and read[6] == True:
                tag = 3
            elif read[4] == False and read[6] == False and read[5] == False and read[7] == True:
                tag = 2
            elif read[4] == False and read[6] == False and read[5] == False and read[7] == False:
                tag = 1
        if tag == 4:
            soft_clip.append(line)
        elif tag == 3:
            hard_clip.append(line)
        elif tag == 2:
            read1.append(line)
        elif tag == 1:
            read2.append(line)
        else:
            normal.append(line)
    i = 0
    for j in range(len(normal)):
        new_dic[i] = normal[j]
        i += 1
    for j in range(len(read1)):
        new_dic[i] = read1[j]
        i += 1
    for j in range(len(read2)):
        new_dic[i] = read2[j]
        i += 1
    for j in range(len(soft_clip)):
        new_dic[i] = soft_clip[j]
        i += 1
    for j in range(len(hard_clip)):
        new_dic[i] = hard_clip[j]
        i += 1
    return new_dic


def read_to_dictionary(read_package, scan_r_pos, height):
    dictionary = {}
    read_list = [0 for i in range(height)]
    for i in range(height):
        read_list[i] = []

    base = []
    read = []
    read_pos1 = []
    read_pos2 = []
    quality = []
    is_concordant = []
    is_clipped = []
    hard_clipped = []
    is_read1 = []

    'check if any read is paired'
    read_package_new = []
    for i in range(len(read_package)):
        if read_package[i][1].is_paired:
            read_package_new.append(read_package[i])

    read_package = read_package_new

    for i in range(len(read_package)):
        base.append(read_package[i][0])
        quality.append(read_package[i][1].mapping_quality)
        read.append(read_package[i][1])

        if read[i].cigartuples[0][0] == 4 or read[i].cigartuples[0][0] == 5:
            read_pos1.append(read[i].reference_start - read[i].cigartuples[0][1])
            read_pos2.append(read_pos1[-1] + read_inferred_len(read[i]))
        else:
            read_pos1.append(read[i].reference_start)
            read_pos2.append(read_pos1[-1] + read_inferred_len(read[i]))
        if read[i].is_paired:  # read not paired，do not consider this kind of read
            is_concordant.append(read[i].is_proper_pair)
            is_read1.append(read[i].is_read1)
            if 'S' in read[i].cigarstring:
                is_clipped.append(True)
            else:
                is_clipped.append(False)

            if 'H' in read[i].cigarstring:
                hard_clipped.append(True)
            else:
                hard_clipped.append(False)
        else:
            continue

    for i in range(height):
        if not read_package:
            break

        for j in range(len(read_package)):
            if read_package == [] or j >= len(read_package):
                break
            Len = len(read_package)
            if not read_list[i]:
                read_list[i].append((read_pos1.pop(0), read_pos2.pop(0),
                                     base.pop(0), quality.pop(0), is_clipped.pop(0),
                                     is_concordant.pop(0), hard_clipped.pop(0), is_read1.pop(0)))
                read.pop(0)
                read_package.pop(0)
                if j >= len(read_package):
                    break
                else:
                    j = j - 1

            elif read_list[i][-1][1] < read_pos1[j] < scan_r_pos:
                read_list[i].append((read_pos1.pop(j), read_pos2.pop(j),
                                     base.pop(j), quality.pop(j), is_clipped.pop(j),
                                     is_concordant.pop(j), hard_clipped.pop(j), is_read1.pop(j)))
                read.pop(j)
                read_package.pop(j)
                if j >= len(read_package):
                    break
                else:
                    j = j - 1

    for i in range(height):
        if read_list[i]:
            dictionary[i] = read_list[i]

    height_new = int(len(dictionary))
    return dictionary, height_new


def draw_empty(which_bp, dic, width, height, scan_l_pos, scan_r_pos, img_name):
    new_im = Image.new("RGB", (width, height), (255, 255, 255))
    new_im.save(img_name, "PNG")
    res_im = new_im.resize((100, 100))
    res_im = new_im.resize((100, 100), Image.ANTIALIAS)
    res_im.save(img_name + '.resize.png', "PNG")


def draw_pgn(which_bp, dic, width, height, scan_l_pos, scan_r_pos, img_name):
    new_im = Image.new("RGB", (width, height), (255, 255, 255))
    for key in range(height):
        for read_tuple in dic[key]:
            read_pos1 = read_tuple[0]
            read_pos2 = read_tuple[1]
            base = read_tuple[2]
            quality = read_tuple[3]
            is_clipped = read_tuple[4]
            is_concordant = read_tuple[5]
            hard_clipped = read_tuple[6]
            is_read1 = read_tuple[7]
            col = read_pos1 - scan_l_pos
            index_in_read = 0
            for i in range(len(base)):
                if 0 <= col < width:
                    row = key
                    red, green, blue = get_RGB(which_bp, base[index_in_read], quality, is_clipped, is_concordant,
                                               hard_clipped, is_read1)
                    new_im.putpixel((col, row), (red, green, blue))
                    index_in_read = index_in_read + 1
                    col = col + 1
                elif col < 0:
                    index_in_read = index_in_read + 1
                    col = col + 1

    # new_im.save(img_name, "PNG")
    res_im = new_im.resize((100, 100))

    if width > 2060:
        new_left = new_im.crop((0, 0, 1000, height)).resize((20, 100), Image.ANTIALIAS)
        new_right = new_im.crop((width - 1000, 0, width, height)).resize((20, 100), Image.ANTIALIAS)
        new_middle = new_im.crop((1000, 0, width - 1000, height)).resize((60, 100), Image.ANTIALIAS)
        res_im.paste(new_left, (0, 0, 20, 100))
        res_im.paste(new_middle, (20, 0, 80, 100))
        res_im.paste(new_right, (80, 0, 100, 100))
        res_im.save(img_name, "PNG")
    else:
        res_im = new_im.resize((100, 100), Image.ANTIALIAS)
        res_im.save(img_name, "PNG")


def get_range(b1, b2):
    tmp = str(abs(b2 - b1))
    high_bit = int(tmp[0]) + 1
    new_tmp = str(high_bit)
    for i in range(1, len(tmp)):
        new_tmp = new_tmp + '0'
    new_tmp = int(new_tmp)
    width = new_tmp + 400
    height = 150
    scan_l_pos = int(b1 - (new_tmp - (b2 - b1)) / 2 - 200)
    scan_r_pos = scan_l_pos + width

    return width, height, scan_l_pos, scan_r_pos
