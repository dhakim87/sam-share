import numpy as np
import matplotlib.pyplot as plt
import sqlite3


def build_priv_mat(gotus, counts_map):
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


def plot(title, gotus, names=None):
    mat_akk_priv = build_priv_mat(gotus, gotu_to_priv)
    mat_akk_split = build_priv_mat(gotus, gotu_to_split)
    mat_akk = build_shared_mat(gotus, gotu_to_shared)

    if np.max(mat_akk) == 0:
        return

    plt.subplots(1,3, gridspec_kw={'width_ratios': [1, 1, 5]})
    plt.subplots_adjust(wspace=0.5)
    plt.subplot(1, 3, 1)
    plt.xticks([])
    plt.imshow(mat_akk_priv, cmap='hot', interpolation='nearest')
    plt.title("Private Reads")

    plt.colorbar()
    plt.subplot(1, 3, 2)
    plt.xticks([])
    plt.imshow(mat_akk_split, cmap='hot', interpolation='nearest')
    plt.title("Split Reads")

    plt.colorbar()
    plt.subplot(1, 3, 3)
    plt.imshow(mat_akk, cmap='hot', interpolation='nearest')
    if names is not None:
        plt.yticks(range(len(names)), names)
    plt.title("Shared2 Confusion Matrix")
    plt.colorbar()
    plt.suptitle(title)
    plt.tight_layout()
    plt.show()
    print(mat_akk_priv)
    print(mat_akk)


with open("S.71801.0038.outsam") as f:
    first_part = True
    all_gotus_set = set()
    gotu_to_priv = {}
    gotu_to_split = {}
    gotu_to_shared = {}
    gotu_to_species = {}
    part = 0
    for line in f:
        if line.startswith("---"):
            part+=1
            continue

        if part == 0:
            ss = line.strip().split(",")
            gotu, private_count = ss[0], int(ss[1])
            gotu_to_priv[gotu] = private_count
            all_gotus_set.add(gotu)
        elif part == 1:
            ss = line.strip().split(",")
            gotu, split_count = ss[0], float(ss[1])
            gotu_to_split[gotu] = split_count
            all_gotus_set.add(gotu)
        else:
            ss = line.strip().split(",")
            gotu1,gotu2,count = ss[0], ss[1], float(ss[2])
            all_gotus_set.add(gotu1)
            all_gotus_set.add(gotu2)
            gotu_to_shared[(gotu1, gotu2)] = count

    all_gotus = sorted(list(all_gotus_set))
    mat_all = build_shared_mat(all_gotus, gotu_to_shared)

    conn = sqlite3.connect("mimics.db")
    c = conn.cursor()
    c.execute("select distinct genus from genome order by genus")
    genomes = [x[0] for x in c.fetchall()]
    c.execute("select genome_id, species from genome")
    for r in c.fetchall():
        gotu_to_species[r[0]] = r[1]

    plt.imshow(mat_all, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title("All")
    plt.show()

    for genome in genomes:
        c.execute("select genome_id, species from genome where genus = ? order by species", (genome,))
        rows = c.fetchall()
        genome_gotus = [r[0] for r in rows if r[0] in all_gotus_set]

        if len(genome_gotus) <= 1:
            continue

        if genome == "":
            title = "Unknown Genus"
        else:
            title = genome
        plot(title, genome_gotus)

    conn.close()

    print("-------------------")
    gotu_ratios = []
    for gotu in all_gotus:
        priv = gotu_to_priv.get(gotu, 0)
        split = gotu_to_split.get(gotu, 0)

        denominator = priv+split
        if denominator == 0:
            ratio = 0
        else:
            ratio = priv/denominator

        name = gotu_to_species[gotu]
        gotu_ratios.append((name, ratio, gotu, priv, split))

    gotu_ratios.sort(key=lambda x: x[1])

    for tup in gotu_ratios:
        if tup[4] < 100:
            continue
        print(tup)

    most_confused = [x for x in gotu_ratios if x[4] >= 1000]
    most_confused.sort(key=lambda x: x[1])

    plot("Most Confused", [x[2] for x in most_confused], names=[x[0] for x in most_confused])
