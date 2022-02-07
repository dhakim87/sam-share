import numpy as np
import matplotlib
from numpy.linalg import LinAlgError

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import sqlite3
import glob
from collections import defaultdict
from os.path import basename
import pandas as pd
from pandas.plotting import parallel_coordinates
import math
from itertools import chain, combinations


# https://stackoverflow.com/questions/1482308/how-to-get-all-subsets-of-a-set-powerset
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def build_priv_mat(gotus, counts_map):
    if counts_map is None:
        return None

    mat = []
    for gotu in gotus:
        row = []
        if gotu in counts_map:
            row.append(counts_map[gotu])
        else:
            row.append(0)

        mat.append(row)
    mat = np.array(mat)
    return mat


def build_shared_mat(gotus, counts_map, gotus_x=None, gotus_y=None):
    mat_all = []
    if gotus_x is None:
        gotus_x = gotus
    if gotus_y is None:
        gotus_y = gotus

    for gotu1 in gotus_y:
        row = []
        for gotu2 in gotus_x:
            count = 0
            if (gotu1, gotu2) in counts_map:
                count = counts_map[(gotu1, gotu2)]
            row.append(count)
        mat_all.append(row)
    mat = np.array(mat_all)

    return mat


def build_overlap_mat_2D(gotus, overlap_map, length_map, gotus_x=None, gotus_y=None):
    mat_all = []
    if gotus_x is None:
        gotus_x = gotus
    if gotus_y is None:
        gotus_y = gotus

    for gotu1 in gotus_y:
        row = []
        for gotu2 in gotus_x:
            if (gotu1, gotu2) in overlap_map:
                gotu1_ranges = overlap_map[(gotu1, gotu2)]
                gotu1_length = length_map[gotu1]
                total_len = 0
                for r in gotu1_ranges:
                    total_len += r[1] - r[0] + 1
                frac_len = total_len / gotu1_length
                row.append(frac_len)
            else:
                row.append(0)
        mat_all.append(row)
    mat = np.array(mat_all)
    return mat

def build_overlap_mat_powerset(gotus, overlap_map, length_map, gotus_x=None, gotus_y=None):
    mat_all = []
    if gotus_x is None:
        gotus_x = gotus
    if gotus_y is None:
        gotus_y = gotus

    for gotu1 in gotus_y:
        row = []
        for gotu2_set in powerset(gotus_x):
            gotu2 = ";".join(sorted(list(gotu2_set)))
            if (gotu1, gotu2) in overlap_map:
                gotu1_ranges = overlap_map[(gotu1, gotu2)]
                gotu1_length = length_map[gotu1]
                total_len = 0
                for r in gotu1_ranges:
                    total_len += r[1] - r[0] + 1
                frac_len = total_len / gotu1_length
                row.append(frac_len)
            else:
                row.append(0)
        mat_all.append(row)
    mat = np.array(mat_all)
    return mat


# Pairwise overlap, makes lots of assumptions about independence.
def calc_distribution_2D(source_index, overlap):
    # Make lots of false assumptions about how overlap is distributed, whee
    sums = [0] * overlap.shape[0]
    for s in powerset(range(overlap.shape[0])):
        if len(s) == 0:
            continue
        p = 1
        for i in range(overlap.shape[0]):
            if i in s:
                p *= overlap[source_index, i] / overlap[source_index, source_index]
            else:
                p *= (1 - overlap[source_index, i] / overlap[source_index, source_index])
        # for ind in s:
        #     p *= overlap[source_index, ind]
        for ind in s:
            sums[ind] += p * 1 / len(s)
    return sums


def calc_distribution_powerset(source_index, overlap):
    # rows should sum to 1.
    overlap[source_index,:] /= sum(overlap[source_index,:])
    # print("BUT WHY!?")

    sums = [0] * overlap.shape[0]
    all_sets = list(powerset(range(overlap.shape[0])))
    # print(all_sets)
    # print(overlap[source_index,:])
    # print(sum(overlap[source_index,:]))
    source_set_index = all_sets.index((source_index,))

    for j in range(overlap.shape[1]):
        s = all_sets[j]
        if len(s) == 0:
            continue
        p = overlap[source_index, j]  #/ overlap[source_index, source_set_index]
        for ind in s:
            sums[ind] += p * 1 / len(s)
    return sums


class GOTUData:
    def __init__(self):
        self.all_gotus_set = set()
        self.gotu_to_priv = {}
        self.gotu_to_split = {}
        self.gotu_to_shared = {}
        self.gotu_to_overlap = defaultdict(list)
        self.gotu_to_coverage = None
        self.gotu_to_length = None

    def set_coverage(self, gotu_to_coverage):
        self.gotu_to_coverage = gotu_to_coverage

    def set_length(self, gotu_to_length):
        self.gotu_to_length = gotu_to_length

    # https://github.com/ucsd-cmi/zebra_filter/blob/master/cover.py
    def compress_ranges(self):
        for gotu1, gotu2 in self.gotu_to_overlap:
            gotu1_ranges = self.gotu_to_overlap[(gotu1,gotu2)]
            sorted_ranges = sorted(gotu1_ranges, key=lambda x:x[0])

            new_ranges = []
            start_val = -1
            end_val = -1

            for range in sorted_ranges:
                if end_val == -1:
                    # case 1: no active range, start active range.
                    start_val = range[0]
                    end_val = range[1]
                elif end_val >= range[0] - 1:
                    # case 2: active range continues through this range
                    # extend active range
                    if end_val < range[1]:
                        end_val = range[1]
                else:
                    # if end_val < r.start - 1:
                    # case 3: active range ends before this range begins
                    # write new range out, then start new active range
                    new_ranges.append((start_val,end_val))
                    start_val = range[0]
                    end_val = range[1]
            if end_val != -1:
                new_ranges.append((start_val, end_val))
            self.gotu_to_overlap[gotu1, gotu2] = new_ranges


def plot(gotu_data, title, gotus, names=None, **kwargs):
    mat_akk_priv = build_priv_mat(gotus, gotu_data.gotu_to_priv)
    mat_akk_split = build_priv_mat(gotus, gotu_data.gotu_to_split)
    mat_akk_cover = build_priv_mat(gotus, gotu_data.gotu_to_coverage)
    mat_akk = build_shared_mat(gotus, gotu_data.gotu_to_shared)

    if np.max(mat_akk, initial=0) == 0:
        return

    numplots = 3
    if mat_akk_cover is not None:
        numplots = 4

    plt.subplots(1, numplots, gridspec_kw={'width_ratios': [8] + [1] * (numplots-1)})
    plt.subplot(1, numplots, 2)
    plt.xticks([])
    plt.imshow(mat_akk_priv, cmap='hot', interpolation='nearest', **kwargs)
    plt.title("Private Reads")

    plt.colorbar()
    plt.subplot(1, numplots, 3)
    plt.xticks([])
    plt.imshow(mat_akk_split, cmap='hot', interpolation='nearest', **kwargs)
    plt.title("Split Reads")

    if mat_akk_cover is not None:
        plt.colorbar()
        plt.subplot(1, numplots, 4)
        plt.xticks([])
        plt.imshow(mat_akk_cover, cmap='hot', interpolation='nearest', vmin=0, vmax=1)
        plt.title("Coverage")

    if names is not None:
        for i in range(len(names)):
            priv = mat_akk_priv[i][0]
            split = mat_akk_split[i][0]
            shared_in_table = sum(mat_akk[i])
            split_ratio = split / (priv + split)
            in_table = shared_in_table / split
            print(i, names[i], mat_akk_priv[i], mat_akk_split[i], sum(mat_akk[i]), "SPLIT:", split_ratio, "IN_TABLE:", in_table)

    plt.colorbar()
    plt.subplot(1, numplots, 1)
    plt.imshow(mat_akk, cmap='hot', interpolation='nearest', **kwargs)
    # plt.spy(mat_akk, precision=100, markersize=4)
    if names is not None and len(names) <= 50:
        plt.yticks(range(len(names)), names)
    plt.title("Share Matrix")
    plt.colorbar()
    plt.suptitle(title)
    plt.show()


def plot_MOM(gotu_data, title, gotus, names=None):
    mat_akk_overlap = build_overlap_mat_powerset(gotus, gotu_data.gotu_to_overlap, gotu_data.gotu_to_length)
    print(mat_akk_overlap)
    plt.imshow(mat_akk_overlap, cmap='hot', interpolation='nearest')
    if names is not None and len(names) <= 50:
        plt.yticks(range(len(names)), names)
    plt.title(title)
    plt.colorbar()
    plt.show()

def akkermansia_parallel_plots(gotu_data_by_file, merged_gotu_data):
    mpl_colors = ['r','g','b','c','m','y','k', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive']
    akk_priv = defaultdict(list)
    akk_split = defaultdict(list)
    akk_total = defaultdict(list)
    c_disease_status = []
    c_locale = []
    c_fit_line = []
    akk_ids = [
        # "G001683795",
        # "G900097105",
        "G000020225",
        "G000723745",
        "G000436395",
        "G001917295",
        "G000437075",
        "G001647615",
        "G001578645",
        "G001580195",
        "G001940945",
        "G000980515"
    ]

    for fname in gotu_data_by_file:
        first_part = basename(fname)[:7]
        if first_part.startswith("Q.") or first_part.startswith("S."):
            disease_status = 'g' if first_part[6] == '1' else 'r'
            locale = first_part[2:5]
            locale_lookup = {"714":'r', "715":'g', "716":'b', "717":'c', "718":'m'}
            c_disease_status.append(disease_status)
            c_locale.append(locale_lookup.get(locale, "k"))
        else:
            continue

        privs = gotu_data_by_file[fname].gotu_to_priv
        splits = gotu_data_by_file[fname].gotu_to_split
        shareds = gotu_data_by_file[fname].gotu_to_shared
        for gid in akk_ids:
            if gid in privs:
                priv = privs[gid]
            else:
                priv = 0
            if gid in splits:
                split = splits[gid]
            else:
                split = 0
            total_share = 0
            for (g1, g2) in shareds:
                if g1 == gid:
                    total_share += shareds[(g1,g2)]

            total = priv + total_share
            akk_priv[gid].append(priv)
            akk_split[gid].append(split)
            akk_total[gid].append(total)

    akk_priv_angles = {}
    for akk_id in akk_ids[1:]:
        x = akk_priv[akk_ids[0]]
        y = akk_priv[akk_id]
        angles = []
        for i in range(len(x)):
            angle = math.degrees(math.atan2(y[i], x[i]))
            angles.append(angle)
        akk_priv_angles[akk_id] = angles

    color_col = pd.Series(c_disease_status)
    df_priv = pd.DataFrame.from_dict(akk_priv_angles)
    df_priv["disease"] = color_col
    # df_split = pd.DataFrame.from_dict(akk_split)
    # df_split["disease"] = color_col
    # df_total = pd.DataFrame.from_dict(akk_total)
    # df_total["disease"] = color_col

    # print(df_priv)
    # parallel_coordinates(df_priv, class_column="disease")
    # plt.show()

    # parallel_coordinates(df_split, class_column="disease")
    # plt.show()
    #
    # parallel_coordinates(df_total, class_column="disease")
    # plt.show()

    clusters = []
    for sample_id in range(len(akk_priv[akk_ids[0]])):
        akk_priv_counts = []
        for akk_id in akk_ids:
            akk_priv_counts.append(akk_priv[akk_id][sample_id])
        max_priv = max(akk_priv_counts)
        max_index = akk_priv_counts.index(max_priv)
        clusters.append(mpl_colors[max_index])

    mat_akk_overlap = build_overlap_mat_powerset(akk_ids, merged_gotu_data.gotu_to_overlap, merged_gotu_data.gotu_to_length)

    print(mat_akk_overlap)
    # akk_lines = [calc_distribution_2D(i, mat_akk_overlap) for i in range(len(akk_ids))]
    akk_lines = [calc_distribution_powerset(i, mat_akk_overlap) for i in range(len(akk_ids))]
    for line in akk_lines:
        print(line)
        print(sum(line))

    mom = np.array(akk_lines)
    print(mom)

    try:
        mom_inv = np.linalg.inv(mom)
        print("Inverted!")
        print(mom_inv)
    except LinAlgError as e:
        print("Couldn't invert, oops")
        mom_inv = None

    vectors = [np.array(line) for line in akk_lines]
    vectors = [v / np.linalg.norm(v) for v in vectors]

    residuals = []
    corrected_points = []
    scaled_residuals = []
    all_scaled_residuals = []
    for sample_id in range(len(akk_priv[akk_ids[0]])):
        sample_point = [akk_total[akk_ids[i]][sample_id] for i in range(len(akk_ids))]
        sample_point = np.array(sample_point)
        if mom_inv is not None:
            corrected_point = np.matmul(mom_inv, sample_point)
            corrected_points.append(corrected_point)
        min_index = 0
        min_dist = 999999999
        for index in range(len(vectors)):
            v = vectors[index]
            # compute distance between sample point and presumed regression line if sample was all
            # all from one akk id
            dot_prod = np.dot(sample_point, v)
            pt_on_line = dot_prod * v
            dist = np.linalg.norm(pt_on_line - sample_point)
            if dist < min_dist:
                min_dist = dist
                min_index = index
        c_fit_line.append(mpl_colors[min_index])

        num_akkermansia_in_sample = np.linalg.norm(sample_point)
        if num_akkermansia_in_sample > 500:
            residuals.append(min_dist)
            scaled_residuals.append(min_dist / num_akkermansia_in_sample)
            all_scaled_residuals.append(min_dist / num_akkermansia_in_sample)
        else:
            all_scaled_residuals.append(0)

    plt.hist(residuals)
    plt.title("Sample residuals- Only 1 Akkermansia")
    plt.show()

    plt.hist(scaled_residuals)
    plt.title("Sample residuals scaled - Only 1 Akkermansia")
    plt.show()

    corrected_points = np.array(corrected_points)
    for i in range(len(akk_ids)):
        for j in range(i + 1, len(akk_ids)):
            plt.subplot(2, 2, 1)
            plt.scatter(akk_priv[akk_ids[i]], akk_priv[akk_ids[j]], c=clusters)
            plt.axis('equal')
            plt.title('private reads')
            plt.xlabel(akk_ids[i])
            plt.ylabel(akk_ids[j])

            plt.subplot(2, 2, 2)
            plt.scatter(akk_split[akk_ids[i]], akk_split[akk_ids[j]], c=clusters)
            plt.axis('equal')
            plt.title('split reads')
            plt.xlabel(akk_ids[i])
            plt.ylabel(akk_ids[j])

            plt.subplot(2, 2, 3)
            plt.scatter(akk_total[akk_ids[i]], akk_total[akk_ids[j]], c=clusters)
            maxx = max(akk_total[akk_ids[i]])
            maxy = max(akk_total[akk_ids[j]])
            for line_index in range(len(akk_lines)):
                line = akk_lines[line_index]
                color = mpl_colors[line_index]
                x = maxx
                y = line[j] * maxx / line[i]
                if y > maxy:
                    x = line[i] * maxy / line[j]
                    y = maxy

                plt.plot([0, x], [0, y], c=color, label=akk_ids[line_index])

            plt.axis('equal')
            plt.title('total reads')
            plt.xlabel(akk_ids[i])
            plt.ylabel(akk_ids[j])
            plt.legend()

            plt.subplot(2, 2, 4)
            plt.scatter(corrected_points[:,i], corrected_points[:,j], c=clusters)
            plt.axis('equal')
            plt.title('corrected')
            plt.xlabel(akk_ids[i])
            plt.ylabel(akk_ids[j])
            plt.legend()
            plt.show()


def load_gotu_to_species():
    gotu_to_species = {}
    conn = sqlite3.connect("mimics.db")
    c = conn.cursor()
    c.execute("select genome_id, species from genome")
    for r in c.fetchall():
        gotu_to_species[r[0]] = r[1]
    conn.close()
    return gotu_to_species


def load_gotu_to_coverage():
    gotu_to_coverage = {}
    conn = sqlite3.connect("mimics.db")
    c = conn.cursor()
    c.execute("select genome_id, coverage_ratio from zebra")
    for r in c.fetchall():
        gotu_to_coverage[r[0]] = r[1]
    conn.close()
    return gotu_to_coverage


def load_gotu_to_length():
    gotu_to_length = {}
    conn = sqlite3.connect("mimics.db")
    c = conn.cursor()
    c.execute("select genome_id, total_length from zebra")
    for r in c.fetchall():
        gotu_to_length[r[0]] = r[1]
    conn.close()
    return gotu_to_length


def load_gotu_data(fname):
    with open(fname) as f:
        gotu_data = GOTUData()
        part = 0

        p3_range_state = ["","",0]
        for line in f:
            if line.startswith("---"):
                part+=1
                continue

            if part == 0:
                ss = line.strip().split(",")
                gotu, private_count = ss[0], int(ss[1])
                gotu_data.gotu_to_priv[gotu] = private_count
                gotu_data.all_gotus_set.add(gotu)
            elif part == 1:
                ss = line.strip().split(",")
                gotu, split_count = ss[0], float(ss[1])
                gotu_data.gotu_to_split[gotu] = split_count
                gotu_data.all_gotus_set.add(gotu)
            elif part == 2:
                ss = line.strip().split(",")
                gotu1,gotu2,count = ss[0], ss[1], float(ss[2])
                gotu_data.all_gotus_set.add(gotu1)
                gotu_data.all_gotus_set.add(gotu2)
                gotu_data.gotu_to_shared[(gotu1, gotu2)] = count
            elif part == 3:
                if p3_range_state[2] == 0:
                    ss = line.strip().split(",")
                    gotu1,gotu2_powerset,num_ranges = ss[0],ss[1], int(ss[2])
                    p3_range_state = [gotu1,gotu2_powerset,num_ranges]
                    gotu_data.all_gotus_set.add(gotu1)
                    for gotu in gotu2_powerset.split(";"):
                        gotu_data.all_gotus_set.add(gotu)
                else:
                    ss = line.strip().split("-")
                    gotu1_start, gotu1_end = int(ss[0]), int(ss[1])
                    gotu_data.gotu_to_overlap[(p3_range_state[0], p3_range_state[1])].append((gotu1_start, gotu1_end))
                    p3_range_state[2] -= 1

    if p3_range_state[2] != 0:
        print("OOPSIE, BAD OVERLAP DATA!?")
        exit(-1)

    return gotu_data


def write_gotu_data(gotu_data, fname):
    with open(fname, "w") as f:
        for gotu in gotu_data.gotu_to_priv:
            f.write(gotu + "," + str(gotu_data.gotu_to_priv[gotu]) + "\n")
        f.write("---\n")
        for gotu in gotu_data.gotu_to_split:
            f.write(gotu + "," + str(gotu_data.gotu_to_split[gotu]) + "\n")
        f.write("---\n")
        for tup in gotu_data.gotu_to_shared:
            gotu1 = tup[0]
            gotu2 = tup[1]
            f.write(gotu1 + "," + gotu2 + "," + str(gotu_data.gotu_to_shared[tup]) + "\n")
        f.write("---\n")
        for tup in gotu_data.gotu_to_overlap:
            gotu1 = tup[0]
            gotu2 = tup[1]
            f.write(gotu1 + "," + gotu2 + "," + str(len(gotu_data.gotu_to_overlap[tup])) + "\n")
            for start,end in gotu_data.gotu_to_overlap[tup]:
                f.write(str(start) + "-" + str(end) + "\n")


def merge_gotu_data(gotu_datas, accumulator=None):
    merged = accumulator
    if merged is None:
        merged = GOTUData()

    def merge_map(gotu_data_map, merged_map):
        for gotu in gotu_data_map:
            if gotu in merged_map:
                merged_map[gotu] += gotu_data_map[gotu]
            else:
                merged_map[gotu] = gotu_data_map[gotu]

    for gotu_data in gotu_datas:
        for gotu in gotu_data.all_gotus_set:
            merged.all_gotus_set.add(gotu)

        merge_map(gotu_data.gotu_to_priv, merged.gotu_to_priv)
        merge_map(gotu_data.gotu_to_split, merged.gotu_to_split)
        merge_map(gotu_data.gotu_to_shared, merged.gotu_to_shared)
        merge_map(gotu_data.gotu_to_overlap, merged.gotu_to_overlap)

    return merged


def main(fnames):
    print("Number of files: ", len(fnames))
    print("Loading Data")
    gotu_to_species = load_gotu_to_species()
    gotu_to_coverage = load_gotu_to_coverage()
    gotu_to_length = load_gotu_to_length()
    i = 0

    file_stats = []
    gotu_data = GOTUData()
    gotu_data_by_file = {}
    for fname in fnames:
        i += 1
        if i % 10 == 0:
            print(i)
        gotu_data_from_file = load_gotu_data(fname)
        gotu_data_by_file[fname] = gotu_data_from_file

        total_priv = sum(gotu_data_from_file.gotu_to_priv.values())
        total_split = sum(gotu_data_from_file.gotu_to_split.values())
        total_reads = total_priv + total_split
        if total_reads == 0:
            print(fname, "has no reads?")
            continue
        percent_split = (total_split / (total_priv + total_split)) * 100
        file_stats.append((total_priv, total_split, total_reads, percent_split))
        gotu_data = merge_gotu_data([gotu_data_from_file], accumulator=gotu_data)

    print("Data Loaded!")

    if len(fnames) > 1:
        pct_split = [x[3] for x in file_stats]
        plt.hist(pct_split)
        plt.title("Read Split Histogram")
        plt.xlabel("Percent Split Reads")
        plt.ylabel("# .sam files")
        plt.show()

        gotu_data.set_coverage(gotu_to_coverage)
        gotu_data.set_length(gotu_to_length)
        gotu_data.compress_ranges()
        print("Writing merged file to merged.outsam")
        write_gotu_data(gotu_data, "merged.outsam")

        akkermansia_parallel_plots(gotu_data_by_file, gotu_data)

    total_priv = sum(gotu_data.gotu_to_priv.values())
    total_split = sum(gotu_data.gotu_to_split.values())
    print("Total Private Reads", total_priv)
    print("Total Split Reads", total_split)
    print("Total Reads", total_priv + total_split)
    print("Percent Reads Split: ", total_split / (total_priv + total_split))

    all_gotus = sorted(list(gotu_data.all_gotus_set), key=lambda x: gotu_to_species[x])
    gotu_names = [gotu_to_species[x] for x in all_gotus]
    for i in range(len(gotu_names)):
        print(i, gotu_names[i])

    mat_all_priv = build_priv_mat(all_gotus, gotu_data.gotu_to_priv)
    mat_all_split = build_priv_mat(all_gotus, gotu_data.gotu_to_split)
    mat_all = build_shared_mat(all_gotus, gotu_data.gotu_to_shared)

    gotu_split_data = []
    for gotu in range(len(all_gotus)):
        priv = mat_all_priv[gotu][0]
        split = mat_all_split[gotu][0]
        fraction_split = split / (priv + split)
        gotu_split_data.append(fraction_split)

    plt.hist(gotu_split_data)
    plt.xlabel("Fraction Of Reads Split Across Multiple Genomes")
    plt.ylabel("# Reference Genomes")
    plt.title("Split Ratio")
    plt.show()

    for precision in [1000, 10000, 100000, 1000000, 10000000]:
        plt.spy(mat_all, precision=precision, markersize=4)
        # plt.colorbar()
        plt.title("Share Matrix Val > " + str(precision))
        plt.show()

    mat_1_in_200000 = []
    mat_1_in_10000 = []
    mat_1_in_1000 = []
    mat_1_in_100 = []
    mat_1_in_10 = []
    mat_all_row_normalized = []
    for y in range(len(all_gotus)):
        norm_row = []
        row_total = mat_all_priv[y][0] + mat_all_split[y][0]
        num_shared_1_in_200000 = 0
        num_shared_1_in_10000 = 0
        num_shared_1_in_1000 = 0
        num_shared_1_in_100 = 0
        num_shared_1_in_10 = 0
        for x in range(len(all_gotus)):
            if (all_gotus[x],all_gotus[y]) in gotu_data.gotu_to_shared:
                val = gotu_data.gotu_to_shared[(all_gotus[x],all_gotus[y])]
            else:
                val = 0
            norm_val = val / row_total
            if norm_val > 0.000005:
                num_shared_1_in_200000 += 1
            if norm_val > 0.0001:
                num_shared_1_in_10000 += 1
            if norm_val > 0.001:
                num_shared_1_in_1000 += 1
            if norm_val > 0.01:
                num_shared_1_in_100 += 1
            if norm_val > 0.1:
                num_shared_1_in_10 += 1
            norm_row.append(val / row_total)

        mat_all_row_normalized.append(norm_row)
        mat_1_in_10.append(num_shared_1_in_10)
        mat_1_in_100.append(num_shared_1_in_100)
        mat_1_in_1000.append(num_shared_1_in_1000)
        mat_1_in_10000.append(num_shared_1_in_10000)
        mat_1_in_200000.append(num_shared_1_in_200000)

    plt.hist(mat_1_in_200000)
    plt.xlabel("Number of genomes sharing 1/200000 of reads")
    plt.ylabel("# Reference Genomes")
    plt.title("Identical Neighbors or Death By 200000 Cuts Read Bleed")
    plt.show()

    plt.hist(mat_1_in_10000)
    plt.xlabel("Number of genomes sharing 1/10000 of reads")
    plt.ylabel("# Reference Genomes")
    plt.title("Identical Neighbors or Death By 10000 Cuts Read Bleed")
    plt.show()

    plt.hist(mat_1_in_1000)
    plt.xlabel("Number of genomes sharing 1/1000 of reads")
    plt.ylabel("# Reference Genomes")
    plt.title("Identical Neighbors or Death By 1000 Cuts Read Bleed")
    plt.show()

    plt.hist(mat_1_in_100)
    plt.xlabel("Number of genomes sharing 1/100 of reads")
    plt.ylabel("# Reference Genomes")
    plt.title("Identical Neighbors or Death By 100 Cuts Read Bleed")
    plt.show()

    plt.hist(mat_1_in_10)
    plt.xlabel("Number of genomes sharing 1/10 of reads")
    plt.ylabel("# Reference Genomes")
    plt.title("Identical Neighbors or Death By 10 Cuts Read Bleed")
    plt.show()

    # mat_all_row_normalized = np.array(mat_all_row_normalized)
    #
    # plt.imshow(mat_all_row_normalized, cmap='hot', interpolation='nearest', vmin=0, vmax=0.0001)
    # plt.title("Read Fraction Clamped to 1/10000")
    # plt.show()
    # for precision in [0.0001, 0.001, 0.01, 0.1, 0.25]:
    #     plt.spy(mat_all_row_normalized, precision=precision, markersize=4)
    #     plt.title("Read Fraction > " + str(precision))
    #     plt.show()

    # bacteroides: 805-855
    # clostridium: 1779-1897
    # Faecalibacterium 2569-2572
    # Cronobacter sakazakii 2107-2108
    range_min = 805
    range_max = 855
    precision = 100000
    mat_all = build_shared_mat(all_gotus, gotu_data.gotu_to_shared, gotus_y=all_gotus[range_min:range_max])
    plt.spy(mat_all, precision=precision, markersize=4, aspect='auto')
    plt.yticks(range(len(gotu_names[range_min:range_max])), gotu_names[range_min:range_max])
    print(mat_all.shape)
    active_genus = ""
    for r in range(mat_all.shape[0]):
        offset_index = 0
        for c in range(mat_all.shape[1]):
            if range_min <= c < range_max:
                continue
            if mat_all[r,c] > precision:
                gname = gotu_names[c]
                if gname.split()[0] != active_genus:
                    rname = gotu_names[r + range_min]
                    cname = gotu_names[c]
                    rpriv = mat_all_priv[r + range_min][0]
                    cpriv = mat_all_priv[c][0]
                    rsplit = mat_all_split[r + range_min][0]
                    csplit = mat_all_split[c][0]
                    val = mat_all[r,c]

                    pct_r_reads = val / (rpriv + rsplit)
                    pct_c_reads = val / (cpriv + csplit)
                    print(rname, cname, pct_r_reads, pct_c_reads)
                    plt.annotate(gotu_names[c] + "\n" + "{:2.0f}".format(val) + "," + "{:2.0f}%".format(pct_r_reads*100) + "," + "{:2.0f}%".format(pct_c_reads*100),
                        xy=(c,r), xycoords='data',
                        xytext=(0, [36, 24, 12][offset_index]), textcoords='offset points',
                        arrowprops=dict(arrowstyle="-", relpos=(0, 0)),
                        horizontalalignment='left',
                        verticalalignment='bottom',
                    )
                    active_genus = gname.split()[0]
                    offset_index = (offset_index + 1) % 3

    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)
    plt.title("Share Matrix Val > " + str(precision))

    plt.show()

    conn = sqlite3.connect("mimics.db")
    c = conn.cursor()
    c.execute("select distinct genus from genome order by genus")
    genera = [r[0] for r in c.fetchall()]
    for genus in genera:
        if genus not in ["Bacteroides"]:
            continue
        c.execute("select genome_id, species from genome where genus = ? order by species", (genus,))
        rows = c.fetchall()
        genus_gotus = [r[0] for r in rows if r[0] in gotu_data.all_gotus_set]

        if len(genus_gotus) <= 1:
            continue

        if genus == "":
            title = "Unknown Genus"
        else:
            title = genus
        plot(gotu_data, title, genus_gotus, names=[gotu_to_species[g] for g in genus_gotus])
        # plot_MOM(gotu_data, title + "- Overlap", genus_gotus, names=[gotu_to_species[g] for g in genus_gotus])

    conn.close()

    print("-------------------")
    gotu_ratios = []
    for gotu in all_gotus:
        priv = gotu_data.gotu_to_priv.get(gotu, 0)
        split = gotu_data.gotu_to_split.get(gotu, 0)

        denominator = priv+split
        if denominator == 0:
            ratio = 0
        else:
            ratio = priv/denominator

        name = gotu_to_species[gotu]
        gotu_ratios.append((name, ratio, gotu, priv, split))

    gotu_ratios.sort(key=lambda x: x[3] + x[4], reverse=True)
    ind = min(25, len(gotu_ratios))
    thresh_total = gotu_ratios[ind][3] + gotu_ratios[ind][4]

    # Show top 10 by total reads assigned
    most_confused = [x for x in gotu_ratios if x[3] + x[4] >= thresh_total]
    most_confused.sort(key=lambda x: x[0])
    names = []
    for i in range(len(most_confused)):
        names.append(most_confused[i][0] + ": " + str(i))
    plot(gotu_data, fname + "-Top GOTUs by Total Reads Assigned", [x[2] for x in most_confused], names=names)


    print("Top 25 Species By Total Reads")
    for x in gotu_ratios[:10]:
        print(x[0], x[3]+x[4])

    gotu_ratios.sort(key=lambda x: x[1])
    zero_priv_count = len([x for x in gotu_ratios if x[3] == 0])
    print(zero_priv_count, " species have only shared reads")
    ind = min(25 + zero_priv_count, len(gotu_ratios))
    # Pick bottom 25 by private reads/total reads
    thresh_ratio = gotu_ratios[ind][1]

    # Show bottom 25 gotus by ratio
    most_confused = [x for x in gotu_ratios if x[1] <= thresh_ratio and x[3] > 0]
    most_confused.sort(key=lambda x: x[0])
    names = []
    for i in range(len(most_confused)):
        names.append(most_confused[i][0] + ": " + str(i))
    plot(gotu_data, fname + "-Most Confused By Ratio (and >1 Private Reads)", [x[2] for x in most_confused], names=names)

    # Show gotus with ratio 0
    most_confused = [x for x in gotu_ratios if x[3] == 0]
    most_confused.sort(key=lambda x: x[0])
    names = []
    for i in range(len(most_confused)):
        names.append(most_confused[i][0] + ": " + str(i))
    plot(gotu_data, fname + "-Most Confused (No private reads)", [x[2] for x in most_confused], names=names)

    gotu_ratios.sort(key=lambda x: x[4], reverse=True)
    ind = min(25, len(gotu_ratios))
    # Pick top 25 by split count
    thresh_split = gotu_ratios[ind][4]

    # Show top 25 gotus by split count
    most_confused = [x for x in gotu_ratios if x[4] >= thresh_split]
    most_confused.sort(key=lambda x: x[0])
    names = []
    for i in range(len(most_confused)):
        names.append(most_confused[i][0] + ": " + str(i))
    plot(gotu_data, fname + "-Most Confused By Split Count", [x[2] for x in most_confused], names=names)

    # Show things that get confused with Yersinia Pestis
    yersinia = [x for x in gotu_ratios if x[0] == "Yersinia pestis"]
    # yersinia = [x for x in gotu_ratios if x[0] == "Staphylococcus aureus"]
    # yersinia = [x for x in gotu_ratios if x[0] == "Bacteroidales bacterium K10"]
    yersinia_gotus = set([y[2] for y in yersinia])
    yersinia_counts = []
    for key in gotu_data.gotu_to_shared:
        if key[0] in yersinia_gotus or key[1] in yersinia_gotus:
            yersinia_counts.append(gotu_data.gotu_to_shared[key])

    yersinia_counts.sort(reverse=True)
    ind = min(50, len(yersinia_counts) - 1)
    yersinia_share_thresh = yersinia_counts[ind]

    yersinia_shared = set([])
    for key in gotu_data.gotu_to_shared:
        if key[0] in yersinia_gotus or key[1] in yersinia_gotus:
            if gotu_data.gotu_to_shared[key] >= yersinia_share_thresh:
                yersinia_shared.add(key[0])
                yersinia_shared.add(key[1])

    # Show Confused with Yersinia pestis
    most_confused = [x for x in gotu_ratios if x[2] in yersinia_shared]
    most_confused.sort(key=lambda x: x[0])
    names = []
    for i in range(len(most_confused)):
        names.append(most_confused[i][0] + ": " + str(i))
    plot(gotu_data, "Confused With Yersinia pestis", [x[2] for x in most_confused], names=names, vmin=0, vmax=50000)
    # plot(gotu_data, "Confused With Staphylococcus aureus", [x[2] for x in most_confused], names=names)
    # plot(gotu_data, "Confused With Bacteroidales bacterium K10", [x[2] for x in most_confused], names=names)


if __name__ == "__main__":
    main(glob.glob("./merged_imsms.outsam"))
    # main(glob.glob("./outsams_akkermansia_coverage/*.outsam")[:100])
    # main(glob.glob("./outsams_akkermansia_powerset_coverage/*.outsam")[:100])
    # main(glob.glob("./187samples_qiita11919.outsam"))
    # main(glob.glob("./outsams_imsms/*.outsam"))
    # main(glob.glob("./outsams_staph_aureus/*.outsam"))
