import inspect
import time

import cooler
import numpy as np
import os
import sys
from DenseUtils import heatmap_with_breakpoints_and_cluster, create_path_if_not_exist, plot_distribution, perf_measure
from DonorRepository import DonorRepository
from DonorService import collect_chr_bins_map_with_resolution, bps_to_bins_with_resolution
from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density

import matplotlib.pyplot as plt

# from HiCUtils import zeros_to_nan, normalize_intra
from functions import functions
from generate_functions import generate_functions, max_range_pow

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
donorspath = script_dir + '/donors'


def analyze_donor(donor, cooler, f_id, f_proximity, rep: DonorRepository, hic_plot):
    print(donor)
    svs = rep.get_donor_sv_chr_1_21(donor)

    resolution = cooler.info['bin-size']
    chr_bins_map = collect_chr_bins_map_with_resolution(svs, resolution)

    # chr_n = 17
    # bin_pairs = chr_bins_map['17']
    for (chr_n, bin_pairs) in chr_bins_map.items():
        # for i in range(1):
        rep.insert_donorinfo(donor, chr_n, f_id)

        all_bins = set()
        number_of_edges = 0
        number_of_nodes = -1
        dist_mat = np.load(f'dists/dists_npy_chr{chr_n}.npy')
        for (x, y) in bin_pairs:
            all_bins.add(x)
            all_bins.add(y)
        all_bins = list(all_bins)

        dirpath = donorspath + f'/{donor}/{chr_n}'
        create_path_if_not_exist(dirpath)

        filepath = donorspath + f'/{donor}/{chr_n}'

        mat = cooler.matrix(balance=False).fetch(f'chr{chr_n}')
        # mat_nan = zeros_to_nan(mat) # better to precount it for every chr
        # mat_norm = normalize_intra(mat_nan)

        with open(filepath, 'w') as outfile:
            for i_idx in range(len(all_bins)):
                for j_idx in range(i_idx + 1, len(all_bins)):
                    i = all_bins[i_idx]
                    j = all_bins[j_idx]
                    number_of_nodes = max(number_of_nodes, max(i, j))
                    dist = dist_mat[i][j]
                    close_f = f_proximity(dist)
                    number_of_edges += 1
                    outfile.write(f'{i} {j} {close_f}\n')

        if number_of_edges == 0:
            continue

        some_delta_just_for_sure = 5
        cluster_bins = WFind_Densest_Subgraph(number_of_nodes + some_delta_just_for_sure, number_of_edges, filepath)
        assert cluster_bins != [], 'not enough accuracy in find densest subgraph algorithm'
        if hic_plot:
            heatmap_with_breakpoints_and_cluster(mat,
                                                 f'normed hic & breakpoints chr{chr_n}\n{inspect.getsource(f_proximity)}',
                                                 bin_pairs,
                                                 cluster_bins,
                                                 save_path=donorspath + f'/{donor}/{chr_n}.png')
        dens = WFind_Density(cluster_bins, filepath)
        print(dens)
        print(f'clusters {cluster_bins}')
        periphery = set()
        info_id = rep.get_info_id(donor, chr_n, f_id)
        rep.insert_cluster(info_id, tuple(cluster_bins))
        rep.insert_dense(info_id, dens)
        for b in bin_pairs:
            i = b[0]
            j = b[1]
            if i in cluster_bins and j not in cluster_bins:
                periphery.add(j)
            if j in cluster_bins and i not in cluster_bins:
                periphery.add(i)
        print(f'periphery: {periphery}')
        rep.insert_periphery(info_id, tuple(list(periphery)))
        print(f'all sv: {all_bins}')


def plot_donor_info(donor, f_id, rep: DonorRepository):
    chr_densities = []
    chr_cluster_size = []
    chr_periphery_size = []
    for chri in range(1, 22):
        info_id = rep.get_info_id(donor, chri, f_id)
        if info_id == -1:
            chr_densities.append(0)
            chr_cluster_size.append(0)
            chr_periphery_size.append(0)
            continue
        cluster = rep.get_cluster(info_id)
        periphery = rep.get_periphery(info_id)
        density = rep.get_density(info_id)
        chr_densities.append(density)
        print(f'{donor} {chri} info cl{cluster} per{periphery} d={density}')
        chr_cluster_size.append(len(cluster))
        chr_periphery_size.append(len(periphery))
    chr_pos = np.arange(21)
    plt.bar(chr_pos, chr_densities, align='center', label='density')
    plt.bar(chr_pos, chr_cluster_size, align='center', alpha=0.2, label='cluster_size')
    plt.bar(chr_pos, chr_periphery_size, align='center', label='periphery_size')
    plt.xlabel('chrN')
    plt.ylabel('info')
    xticks = [f'{i}' for i in range(1, 22)]
    plt.xticks(chr_pos, xticks)
    plt.title(f'{donor} f_id={f_id}')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.show()


def all_info_donors_plots(f_id_arr):
    t1 = time.time()
    with DonorRepository() as rep:
        donors = rep.unique_prostate_donors()
        for donor in donors[:5]:
            for f_id in f_id_arr:
                plot_donor_info(donor, f_id, rep)

    t2 = time.time()
    print(f'plots took {t2 - t1} sec')


# run generate_functions.py before
def main():
    t1 = time.time()

    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    with DonorRepository() as rep:
        rep.ddl()
        generate_functions()
        for f in functions:
            rep.insert_proximity(f)

        donors = rep.unique_prostate_donors()
        # donors = ['A31-0018_CRUK_PC_0018']
        f_ids = [i for i in range(1, max_range_pow + 1)]
        # f_ids = [9]
        for donor in donors:
            for f_id in f_ids:
                analyze_donor(donor=donor, cooler=cool, f_id=f_id, f_proximity=functions[f_id - 1], rep=rep,
                              hic_plot=False)
                print(f'f_id={f_id} filled')

    t2 = time.time()
    print(f'filling db took {t2 - t1} sec')


def hist_patients(f_id):
    t1 = time.time()
    denss_info = []
    clusters = []
    peripheries = []
    infoids = []
    with DonorRepository() as rep:
        donors = rep.unique_prostate_donors()
        for donor in donors:
            for i in range(1, 22):
                info_id = rep.get_info_id(donor, i, f_id)
                if info_id != -1:
                    dens = rep.get_density(info_id)
                    cluster = len(rep.get_cluster(info_id))
                    periphery = len(rep.get_periphery(info_id))
                    # print(cluster)
                    # if cluster == 52:
                    #    print(info_id)
                    # denss.append(dens)
                    denss_info.append((dens, info_id))
                    clusters.append(cluster)
                    peripheries.append(periphery)
        code = rep.get_proximity_code(f_id)

        dens_path = f'distribution/densities/{f_id}.png'
        ratio = 0.99
        denss = list(map(lambda x: x[0], denss_info))
        top_ratio_infos = sorted(denss_info, key=lambda p: p[0])
        cnt = len(denss)
        for i in range(int(cnt * ratio), cnt):
            infoids.append(top_ratio_infos[i][1])
        plot_distribution(denss, dens_path, code, 'density', 'density', ratio=ratio)

        cluster_path = f'distribution/clusters/{f_id}.png'
        plot_distribution(clusters, cluster_path, code, 'clusters', 'cluster_size', 'green', ratio=ratio)

        periphery_path = f'distribution/peripheries/{f_id}.png'
        plot_distribution(peripheries, periphery_path, code, 'peripheries', 'periphery_size', 'orange', ratio=ratio)

    t2 = time.time()
    print(f'plots took {t2 - t1} sec')
    return infoids


def find_chromos(infoids):
    chromos = []
    with DonorRepository() as rep:
        for id in infoids:
            donor, chr, _ = rep.get_by_infoid(id)
            cluster = rep.get_cluster(id)
            low = min(cluster)
            up = max(cluster)
            chromos.append((donor, chr, low, up))
    return chromos


def find_cromo_clusters(infoids):
    chromo_clusters = []
    with DonorRepository() as rep:
        for id in infoids:
            donor, chr, _ = rep.get_by_infoid(id)
            cluster = rep.get_cluster(id)
            chromo_clusters.append((donor, chr, cluster))
    return chromo_clusters


def compare_chromos(chromos):
    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    resolution = cool.info['bin-size']
    with DonorRepository() as rep:
        for donor, chr, low, up in chromos:
            donor_chr_svs_bp = rep.get_donor_chr_svs(donor, chr)
            donor_chr_svs = list(map(lambda p: bps_to_bins_with_resolution(p[0], p[1], resolution), donor_chr_svs_bp))
            print(donor_chr_svs)
            print('my   ', donor, f'chr{chr}', low, up)
            seek_donor, seek_chr, seek_low, seek_up = rep.get_chromo(donor, chr)
            if seek_low and seek_up:
                seek_low, seek_up = bps_to_bins_with_resolution(seek_low, seek_up, resolution)
            print('seek ', donor, f'chr{seek_chr}', seek_low, seek_up)
            print()


def measure_chromos(chromo_clusters):
    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    resolution = cool.info['bin-size']
    accs = []
    chromo_cnt = 0
    all_cnt = 0
    with DonorRepository() as rep:
        for donor, chr, cluster in chromo_clusters:
            all_cnt += 1

            donor_chr_svs_bp = rep.get_donor_chr_svs(donor, chr)
            donor_chr_svs = list(map(lambda p: bps_to_bins_with_resolution(p[0], p[1], resolution), donor_chr_svs_bp))
            donor_chr_svs = list(filter(lambda p: abs(p[0] - p[1]) > 1, donor_chr_svs))
            donor_chr_bps = set()
            for p in donor_chr_svs:
                donor_chr_bps.add(p[0])
                donor_chr_bps.add(p[1])
            donor_chr_bps = list(donor_chr_bps)

            seek_donor, seek_chr, seek_low, seek_up = rep.get_chromo(donor, chr)

            if seek_low == -2:
                print('No Chromo')
                continue

            if seek_low == -1:
                all_cnt -= 1
                print('No Data')
                continue

            seek_low, seek_up = bps_to_bins_with_resolution(seek_low, seek_up, resolution)

            chromo_cnt += 1

            seek_vector = []
            for bp in donor_chr_bps:
                if seek_low <= bp <= seek_up:
                    seek_vector.append(1)
                else:
                    seek_vector.append(0)

            own_vector = []
            cluster_set = set(cluster)
            for bp in donor_chr_bps:
                if bp in cluster_set:
                    own_vector.append(1)
                else:
                    own_vector.append(0)

            TP, FP, TN, FN = perf_measure(seek_vector, own_vector)
            acc = (TP + TN) / (TP + TN + FP + FN)
            accs.append(acc)
            print(acc)
    mean_chromo_acc = np.mean(np.array(accs))
    print(f'mean chromo acc: {mean_chromo_acc * 100}%')
    chromo_marker_acc = chromo_cnt / all_cnt
    print(f'marker acc: {chromo_marker_acc * 100}%')

if __name__ == '__main__':
    # main()
    # for i in range(1, max_range_pow + 1):
    # for i in range(1, 4):

    i = 3
    infoids = hist_patients(i)
    #chromos = find_chromos(infoids)
    chromo_clusters = find_cromo_clusters(infoids)
    measure_chromos(chromo_clusters)

    # with DonorRepository() as rep:
    #     print(rep.get_chromo('CPCG0201', 16))

    # f_id_arr = [2]
    # all_info_donors_plots([1])
