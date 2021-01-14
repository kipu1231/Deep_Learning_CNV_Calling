from PIL import Image, ImageDraw
import numpy as np
import pandas as pd
from KMeansSource.k_means import *
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from sklearn.cluster import KMeans


# generate candidate deletion: get single, split and discordant reads around del, kmeans and sorting false positives
# + breakpoint estimation
# ehemals call_dell
def filter_cnv(vcf_del, sam_file):
    vcf_len = len(vcf_del)
    del_pos = []
    del_left_pos = 0
    del_right_pos = 0
    dup_left_pos = 0
    dup_right_pos = 0

    for i in range(vcf_len):
        if vcf_del[i][3] == "DUP": continue
        print('------------')
        print(str(vcf_del[i][3]) + " Nummer : " + str(i))

        read_depth = get_depth(sam_file, vcf_del[i][0], int(vcf_del[i][1] - 200), int(vcf_del[i][2] + 200))
        seq_depth = []

        len_deletion = len(read_depth)
        for j in range(len_deletion):
            seq_depth.append((int(vcf_del[i][1] - 200 + j), int(read_depth[j])))

        clip_pic = get_clip_num(sam_file, vcf_del[i][0], vcf_del[i][1] - 200, vcf_del[i][2] + 200, vcf_del[i][3])

        seq_depth_array = np.array(seq_depth)
        seq_clip_array = np.array(clip_pic)

        'Plot for visualizing read depth and discordant read pairs'
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
        fig = plt.figure()

        # ax = init_pic(6, 1, 1, fig, '2d')
        # ax.plot(seq_depth_array[:, 0], seq_depth_array[:, 1], color='r')
        # ax.plot(seq_clip_array[:, 0], seq_clip_array[:, 1] * 2, color='g')
        # plt.show()

        'Preprocess and eliminate noise before clustering'
        'Generates clustering input: (position, read depth, discordant read count and split read count)'
        df_depth = pd.DataFrame(seq_depth_array)
        df_clip = pd.DataFrame(seq_clip_array)

        df_merge = pd.merge(df_depth, df_clip, on=[0], how='left')
        df_merge_fill = df_merge.fillna(0)
        df_merge_fill['mean'] = df_merge_fill['1_x'].rolling(window=61, center=True).mean()
        df_merge_fill['std'] = df_merge_fill['1_x'].rolling(window=61, center=True).std()

        df_merge_fill['result'] = (df_merge_fill['mean'] / df_merge_fill['std']) * 5
        # df_merge_fill = clean_dataset(df_merge_fill)

        if i == 5:
            print(df_merge_fill.iloc[[2027]])
            print(type(df_merge_fill['result'][2027]))

        smooth_kmeans = np.c_[df_merge_fill[30:-30][0], df_merge_fill[30:-30]['result']]
        smooth_kmeans = np.c_[smooth_kmeans, df_merge_fill[30:-30]['1_y']]
        smooth_kmeans = pd.DataFrame(smooth_kmeans)
        ind = smooth_kmeans.index[np.isinf(smooth_kmeans).any(1)]

        if not ind.empty:
            for i in ind:
                print(smooth_kmeans.iloc[[i]])

        smooth_kmeans = clean_dataset(pd.DataFrame(smooth_kmeans))
        # print(smooth_kmeans)
        # print(np.all(np.isfinite(smooth_kmeans)))

        'cluster into three groups: before variant, variant, after variant'
        # result_3D = k_means(smooth_kmeans, 3, 100000)
        kmeans = KMeans(n_clusters=3, random_state=0)
        cluster_pred = kmeans.fit_predict(smooth_kmeans)
        result_3D = np.c_[smooth_kmeans, cluster_pred]
        result_3D = np.nan_to_num(result_3D)
        result_3D = pd.DataFrame(result_3D)

        'Dataframe per cluster'
        cluster_one = result_3D.loc[result_3D[3] == 0.0]
        cluster_two = result_3D.loc[result_3D[3] == 1.0]
        cluster_three = result_3D.loc[result_3D[3] == 2.0]

        x = 5

        'Calculate metric'
        c1_mean = cluster_one[1].mean() // x
        c1_min = cluster_one[1].min() // x
        c1_max = cluster_one[1].max() // x

        c2_mean = cluster_two[1].mean() // x
        c2_min = cluster_two[1].min() // x
        c2_max = cluster_two[1].max() // x

        c3_mean = cluster_three[1].mean() // x
        c3_min = cluster_three[1].min() // x
        c3_max = cluster_three[1].max() // x

        'Input'
        class_array_min = [c1_min, c2_min, c3_min]
        class_array_max = [c1_max, c2_max, c3_max]
        class_array = [c1_mean, c2_mean, c3_mean]
        class_sort = np.sort(class_array)

        # del_left_pos = 0
        # del_right_pos = 0
        # dup_left_pos = 0
        # dup_right_pos = 0

        'filter false positives for deletions'
        'check if read depth mean difference between variant and none variant has a sufficient size'
        if vcf_del[i][3] == "DUP":
            dup_left_pos = dup_left_pos + 1
            if class_sort[2] > 4 * class_sort[1] // 5 and class_sort[2] > 4 * class_sort[0] // 5:
                sufficient_difference = True
                dup_right_pos = dup_right_pos + 1
            else:
                sufficient_difference = False
        elif vcf_del[i][3] == "DEL":
            del_left_pos = del_left_pos + 1
            if class_sort[0] < 4 * class_sort[1] // 5 and class_sort[0] < 4 * class_sort[2] // 5:
                sufficient_difference = True
                del_right_pos = del_right_pos + 1
            else:
                sufficient_difference = False

    #     if sufficient_difference:
    #         class_min_index = np.argmin(class_array)
    #         class_max_index = np.argmax(class_array)
    #         class_mid_index = 3 - class_max_index - class_min_index
    #         seq_depth_dict = dict(seq_depth)
    #
    #         if vcf_del[i][3] == "DUP":
    #             cnv_type = 'DUP'
    #             class_index = class_max_index
    #             class_array_cnv = class_array_max
    #         elif vcf_del[i][3] == "DEL":
    #             cnv_type = 'DEL'
    #             class_index = class_min_index
    #             class_array_cnv = class_array_min
    #
    #         #if cluster with lowest mean RD does not belong to before variant class or after variant class
    #         if class_index == 0 and int(class_one_3D[0, 0]) != smooth_kmeans[0, 0] and int(class_one_3D[-1, 0]) != smooth_kmeans[-1, 0]:
    #             del_left_pos = int(class_one_3D[0, 0])
    #             del_right_pos = int(class_one_3D[-1, 0])
    #
    #             while del_left_pos:
    #                 if seq_depth_dict[del_left_pos] > class_array_cnv[class_mid_index] or seq_depth[0][0] == del_left_pos:
    #                     break
    #                 del_left_pos -= 1
    #             while del_right_pos:
    #                 if seq_depth_dict[del_right_pos] > class_array_cnv[class_mid_index] or del_right_pos == \
    #                         seq_depth[-1][0]:
    #                     break
    #                 del_right_pos += 1
    #
    #             del_left_diff = int(class_one_3D[0, 0] - del_left_pos)
    #             del_right_diff = int(del_right_pos - class_one_3D[-1, 0])
    #             del_pos.append((del_left_pos, del_right_pos, int(class_one_3D[0, 0]), int(class_one_3D[-1, 0]), cnv_type))
    #
    #         #if cluster with lowest mean RD does not belong to before variant class or after variant class
    #         elif class_index == 1 and int(class_two_3D[0, 0]) != smooth_kmeans[0, 0] and smooth_kmeans[-1, 0] != int(class_two_3D[-1, 0]):
    #             del_left_pos = int(class_two_3D[0, 0])
    #             del_right_pos = int(class_two_3D[-1, 0])
    #
    #             while del_left_pos:
    #                 if seq_depth_dict[del_left_pos] > class_array_cnv[class_mid_index] or del_left_pos == seq_depth[0][0]:
    #                     break
    #                 del_left_pos -= 1
    #             while del_right_pos:
    #                 if seq_depth_dict[del_right_pos] > class_array_cnv[class_mid_index] or del_right_pos == \
    #                         seq_depth[-1][0]:
    #                     break
    #                 del_right_pos += 1
    #
    #             del_left_diff = int(class_two_3D[0, 0] - del_left_pos)
    #             del_right_diff = int(del_right_pos - class_two_3D[-1, 0])
    #             del_pos.append((del_left_pos, del_right_pos, int(class_two_3D[0, 0]), int(class_two_3D[-1, 0]), cnv_type))
    #
    #         #if cluster with lowest mean RD does not belong to before variant class or after variant class
    #         elif class_index == 2 and int(class_three_3D[0, 0]) != smooth_kmeans[0, 0] and int(class_three_3D[-1, 0]) != smooth_kmeans[-1, 0]:
    #             del_left_pos = int(class_three_3D[0, 0])
    #             del_right_pos = int(class_three_3D[-1, 0])
    #
    #             while del_left_pos:
    #                 if seq_depth_dict[del_left_pos] > class_array_cnv[class_mid_index] or del_left_pos == seq_depth[0][0]:
    #                     break
    #                 del_left_pos -= 1
    #             while del_right_pos:
    #                 if seq_depth_dict[del_right_pos] > class_array_cnv[class_mid_index] or del_right_pos == \
    #                         seq_depth[-1][0]:
    #                     break
    #                 del_right_pos += 1
    #
    #             del_left_diff = int(class_three_3D[0, 0] - del_left_pos)
    #             del_right_diff = int(del_right_pos - class_three_3D[-1, 0])
    #             del_pos.append((del_left_pos, del_right_pos, int(class_three_3D[0, 0]), int(class_three_3D[-1, 0]), cnv_type))
    #
    #         # ax3 = init_pic(6, 1, 3, fig, '3d')
    #         # ax3.scatter(class_one_3D[:, 0], class_one_3D[:, 1], class_one_3D[:, 2], color='r')
    #         # ax3.scatter(class_two_3D[:, 0], class_two_3D[:, 1], class_two_3D[:, 2], color='g')
    #         # ax3.scatter(class_three_3D[:, 0], class_three_3D[:, 1], class_three_3D[:, 2], color='b')
    #         # ax3.set_xlabel('X')
    #         # ax3.set_zlabel('Z')
    #         # ax3.set_ylabel('Y')
    #
    #     #     ax4 = init_pic(6, 1, 4, fig, '2d')
    #     #     ax4.plot(class_one_3D[:, 0], class_one_3D[:, 1], color='r')
    #     #     ax4.plot(class_two_3D[:, 0], class_two_3D[:, 1], color='g')
    #     #     ax4.plot(class_three_3D[:, 0], class_three_3D[:, 1], color='b')
    #     #     ax4.set_xlabel('gene locus')
    #     #     ax4.set_ylabel('read depth')
    #     #
    #     #     print("left_pos %d" % del_left_pos)
    #     #     print("right_pos %d" % del_right_pos)
    #     #
    #     #     if del_right_pos - del_left_pos != 0 and del_right_pos + del_left_pos != 0:
    #     #         del_len = del_right_pos - del_left_pos + 1
    #     #         del_pic = []
    #     #         for del_i in range(del_len):
    #     #             del_pic.append((del_left_pos + del_i, seq_depth_dict[del_left_pos + del_i]))
    #     #
    #     #         del_pic = np.array(del_pic)
    #     #         ax5 = init_pic(6, 1, 5, fig, '2d')
    #     #         ax5.plot(del_pic[:, 0], del_pic[:, 1], color='r')
    #     #         ax5.set_xlabel('gene locus of variant')
    #     #         ax5.set_ylabel('read depth')
    #     #
    #     # fig.set_size_inches(18.5, 20.5)
    #     # fig.suptitle('Read depth around variation')
    #     # plt.show()
    #     # print("last_i= %d" % i)
    #     # if os.path.exists('your file name'):
    #     #     fig.savefig("your file name" + del_name + "_" + str(i) + '_' + str(del_left_pos) + '_' + str(
    #     #         del_right_pos) + ".png")
    #     # else:
    #     #     os.mkdir('your file name')
    #     #     fig.savefig("your file name" + del_name + "_" + str(i) + '_' + str(del_left_pos) + '_' + str(
    #     #         del_right_pos) + ".png")
    #
    #     #plt.close('all')
    #
    print("Summe DUB: " + str(dup_left_pos))
    print("Davon Diff enough: " + str(dup_right_pos))
    print("Summe DEL: " + str(del_left_pos))
    print("Davon Diff enough: " + str(del_right_pos))
    #
    #
    # return del_pos


# Returns depth at variant position
def get_depth(sam_file, chr_id, pos_l, pos_r):
    read_depth = sam_file.count_coverage(chr_id, pos_l, pos_r)
    depth = np.array(list(read_depth)).sum(axis=0)
    depth = list(depth)
    return depth


def get_clip_num(sam_file, chr_id, pos_l, pos_r, cnv_type):
    clip_temp = []
    for read in sam_file.fetch(chr_id, pos_l, pos_r):
        if read.cigarstring is not None:
            base_pos = read.get_reference_positions(True)
            read_len = len(read.get_reference_positions(True))
            index = 0
            for read_map_pos in range(read_len):
                if base_pos[read_map_pos] is not None:
                    break
                else:
                    index += 1
            read_start = read.reference_start - index
            read_end = read_start + read_len - 1

            for i in range(pos_r - pos_l + 1):
                if read_end >= pos_l + i >= read_start:
                    index_ptr = 0
                    map_type = -1
                    for cigar in read.cigartuples:
                        if pos_l + i > index_ptr:
                            index_ptr = cigar[1] + index_ptr
                            map_type = cigar[0]
                    if cnv_type == "DEL":
                        clip_temp.append((pos_l + i, -map_type))
                    elif cnv_type == "DUP":
                        clip_temp.append((pos_l + i, map_type))

    clip_record_np = np.array(clip_temp)
    df = pd.DataFrame(clip_record_np)
    clip_record_df = df.groupby(0).sum()
    clip_record_df = clip_record_df // 4

    temp = clip_record_df.reset_index()
    clip_record = np.array(temp).tolist()
    return clip_record


# functions needed for plotting
def init_pic(row, col, th, fig, flag):
    if flag == '2d':
        ax = fig.add_subplot(row, col, th)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        return ax
    elif flag == '3d':
        ax = fig.add_subplot(row, col, th, projection='3d')
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        return ax


# clean dataframe before inputing in kmeans
def clean_dataset(df):
    print(len(df))
    assert isinstance(df, pd.DataFrame), "df needs to be a pd.DataFrame"
    df.dropna(inplace=True)
    indices_to_keep = ~df.isin([np.nan, np.inf, -np.inf]).any(1)
    print(len(df[indices_to_keep].astype(np.float64)))
    return df[indices_to_keep].astype(np.float64)


