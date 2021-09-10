import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import sqlite3
import glob


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


def build_shared_mat(gotus, counts_map):
    mat_all = []
    for gotu1 in gotus:
        row = []
        for gotu2 in gotus:
            count = 0
            if (gotu1, gotu2) in counts_map:
                count = counts_map[(gotu1, gotu2)]
            row.append(count)
        mat_all.append(row)
    mat = np.array(mat_all)

    return mat


class GOTUData:
    def __init__(self):
        self.all_gotus_set = set()
        self.gotu_to_priv = {}
        self.gotu_to_split = {}
        self.gotu_to_shared = {}
        self.gotu_to_coverage = None

    def set_coverage(self, gotu_to_coverage):
        self.gotu_to_coverage = gotu_to_coverage


def plot(gotu_data, title, gotus, names=None):
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
    plt.imshow(mat_akk_priv, cmap='hot', interpolation='nearest')
    plt.title("Private Reads")

    plt.colorbar()
    plt.subplot(1, numplots, 3)
    plt.xticks([])
    plt.imshow(mat_akk_split, cmap='hot', interpolation='nearest')
    plt.title("Split Reads")

    if mat_akk_cover is not None:
        plt.colorbar()
        plt.subplot(1, numplots, 4)
        plt.xticks([])
        plt.imshow(mat_akk_cover, cmap='hot', interpolation='nearest')
        plt.title("Coverage")

    plt.colorbar()
    plt.subplot(1, numplots, 1)
    plt.imshow(mat_akk, cmap='hot', interpolation='nearest')
    if names is not None:
        plt.yticks(range(len(names)), names)
    plt.title("Confusion Matrix")
    plt.colorbar()
    plt.suptitle(title)
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


def load_gotu_data(fname):
    with open(fname) as f:
        gotu_data = GOTUData()
        part = 0
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
            else:
                ss = line.strip().split(",")
                gotu1,gotu2,count = ss[0], ss[1], float(ss[2])
                gotu_data.all_gotus_set.add(gotu1)
                gotu_data.all_gotus_set.add(gotu2)
                gotu_data.gotu_to_shared[(gotu1, gotu2)] = count

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

    return merged


def main(fnames):
    print("Number of files: ", len(fnames))
    print("Loading Data")
    gotu_to_species = load_gotu_to_species()
    gotu_to_coverage = load_gotu_to_coverage()
    i = 0

    gotu_data = GOTUData()
    for fname in fnames:
        i += 1
        if i % 10 == 0:
            print(i)
        gotu_data_from_file = load_gotu_data(fname)
        gotu_data = merge_gotu_data([gotu_data_from_file], accumulator=gotu_data)

    print("Data Loaded!")

    if len(fnames) > 1:
        print("Writing merged file to merged.outsam")
        write_gotu_data(gotu_data, "merged.outsam")

    gotu_data.set_coverage(gotu_to_coverage)

    all_gotus = sorted(list(gotu_data.all_gotus_set), key=lambda x: gotu_to_species[x])
    mat_all = build_shared_mat(all_gotus, gotu_data.gotu_to_shared)

    for precision in [100, 1000, 10000, 100000, 1000000, 10000000, 100000000]:
        plt.spy(mat_all, precision=precision, markersize=4)
        # plt.colorbar()
        plt.title("Confused > " + str(precision))
        plt.show()

    conn = sqlite3.connect("mimics.db")
    c = conn.cursor()
    c.execute("select distinct genus from genome order by genus")
    genera = [r[0] for r in c.fetchall()]
    for genus in genera:
        if genus not in ["Akkermansia"]:
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
    ind = min(10, len(gotu_ratios))
    thresh_total = gotu_ratios[ind][3] + gotu_ratios[ind][4]

    # Show top 10 by total reads assigned
    most_confused = [x for x in gotu_ratios if x[3] + x[4] >= thresh_total]
    most_confused.sort(key=lambda x: x[0])
    names = []
    for i in range(len(most_confused)):
        names.append(most_confused[i][0] + ": " + str(i))
    plot(gotu_data, fname + "-Top GOTUs by Total Reads Assigned", [x[2] for x in most_confused], names=names)


    print("Top 10 Species By Total Reads")
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
    plot(gotu_data, fname + "-Most Confused (No private reads)", [x[2] for x in most_confused])

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
    yersinia_gotus = set([y[2] for y in yersinia])
    yersinia_counts = []
    for key in gotu_data.gotu_to_shared:
        if key[0] in yersinia_gotus or key[1] in yersinia_gotus:
            yersinia_counts.append(gotu_data.gotu_to_shared[key])

    yersinia_counts.sort(reverse=True)
    ind = min(200, len(yersinia_counts))
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
    plot(gotu_data, "Confused With Yersinia pestis", [x[2] for x in most_confused], names=names)


if __name__ == "__main__":
    main(glob.glob("./imsms_all.outsam"))
    # main(glob.glob("./outsams/*.outsam"))
