import inspect
import math
import time

import cooler
import numpy as np
import os
import sys
from DenseUtils import heatmap_with_breakpoints_and_cluster, create_path_if_not_exist, plot_distribution, perf_measure, \
    plot_seek_distibution, plot_seek_compare
from DonorRepository import DonorRepository
from DonorService import collect_chr_bins_map_with_resolution, bps_to_bins_with_resolution
from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density

import matplotlib.pyplot as plt

from HiCUtils import zeros_to_nan, normalize_intra
from functions import functions
from generate_functions import generate_functions, max_range_pow

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
donorspath = script_dir + '/donors'


def analyze_donor(donor, cooler, f_id, f_proximity, rep: DonorRepository, hic_plot, from_hic=False, oe=False):
    print(donor)
    svs = rep.get_donor_sv_chr_1_21(donor)

    resolution = cooler.info['bin-size']
    chr_bins_map = collect_chr_bins_map_with_resolution(svs, resolution)

    corr_densities = []
    corr_sv_cnt = []
    corr_edges_cnt = []

    for (chr_n, bin_pairs) in chr_bins_map.items():

        all_bins = set()
        number_of_edges = 0
        number_of_nodes = -1

        mat = cooler.matrix(balance=False).fetch(f'chr{chr_n}')

        if not from_hic:
            dist_mat = np.load(f'dists/dists_npy_chr{chr_n}.npy')
        else:
            if not oe:
                dist_mat = mat
            else:
                mat_nan = zeros_to_nan(mat)  # better to precount it for every chr
                mat_norm = normalize_intra(mat_nan)
                dist_mat = mat_norm

        for (x, y) in bin_pairs:
            all_bins.add(x)
            all_bins.add(y)
        all_bins = list(all_bins)

        dirpath = donorspath + f'/{donor}/{chr_n}'
        create_path_if_not_exist(dirpath)

        filepath = donorspath + f'/{donor}/{chr_n}'

        with open(filepath, 'w') as outfile:
            for i_idx in range(len(all_bins)):
                for j_idx in range(i_idx + 1, len(all_bins)):
                    i = all_bins[i_idx]
                    j = all_bins[j_idx]
                    number_of_nodes = max(number_of_nodes, max(i, j))
                    dist = dist_mat[i][j]
                    if dist != np.nan:
                        close_f = f_proximity(dist)
                        if close_f == 0:
                            continue
                        number_of_edges += 1
                        outfile.write(f'{i} {j} {close_f}\n')

        if number_of_edges == 0:
            continue

        some_delta_just_for_sure = 5
        cluster_bins = WFind_Densest_Subgraph(number_of_nodes + some_delta_just_for_sure, number_of_edges, filepath)
        # assert cluster_bins != [], 'not enough accuracy in find densest subgraph algorithm'
        print('ZERO CLUSTER')
        if hic_plot:
            heatmap_with_breakpoints_and_cluster(mat,
                                                 f'normed hic & breakpoints chr{chr_n}\n{inspect.getsource(f_proximity)}',
                                                 bin_pairs,
                                                 cluster_bins,
                                                 save_path=donorspath + f'/{donor}/{chr_n}.png')
        rep.insert_donorinfo(donor, chr_n, f_id)
        dens = WFind_Density(cluster_bins, filepath)
        print(dens)

        corr_densities.append(dens)
        corr_sv_cnt.append(len(bin_pairs))
        corr_edges_cnt.append(number_of_edges)

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

    return corr_densities, corr_sv_cnt, corr_edges_cnt


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
        # f_ids = [i for i in range(1, max_range_pow + 1)] todo: replace here
        f_ids = [3, 4]


        for f_id in f_ids:
            densess = []
            svscnt = []
            edgescnt = []
            corr_path = f'distribution/corr/{f_id}/corr_arrays.npy'
            create_path_if_not_exist(corr_path)
            for donor in donors:
                corr_densities, corr_sv_cnt, corr_edges_cnt = analyze_donor(donor=donor, cooler=cool, f_id=f_id, f_proximity=functions[f_id - 1], rep=rep,
                              hic_plot=False)
                densess.extend(corr_densities)
                svscnt.extend(corr_sv_cnt)
                edgescnt.extend(corr_edges_cnt)
            print(f'f_id={f_id} filled')
            with open(corr_path, 'wb') as fl:
                np.save(fl, np.array(densess))
                np.save(fl, np.array(svscnt))
                np.save(fl, np.array(edgescnt))

    t2 = time.time()
    print(f'filling db took {t2 - t1} sec')


def hist_patients(f_id, ratio, dens_plot=True, cluster_plot=True, periphery_plot=True, seek_plot=True):
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
                    denss_info.append((dens, info_id))
                    clusters.append(cluster)
                    peripheries.append(periphery)
        code = rep.get_proximity_code(f_id)

        denss = list(map(lambda x: x[0], denss_info))
        top_ratio_infos_sorted_by_dens = sorted(denss_info, key=lambda p: p[0])
        infoids_sorted_by_dens = list(map(lambda p: p[1], top_ratio_infos_sorted_by_dens))
        cnt = len(denss)

        for i in range(int(cnt * ratio), cnt):
            cur_id = top_ratio_infos_sorted_by_dens[i][1]
            infoids.append(cur_id)

        ratios = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
        ratios_map = chrs_per_donor_summary(infoids_sorted_by_dens, ratios)
        text_ratios_map = ''
        for k, v in ratios_map.items():
            text_ratios_map += f'ratio={k} : ' + 'chr_per_donor=%0.2f\n' % v

        buckets_cnt = 100
        dens_path = f'distribution/densities/{f_id}.png'
        if dens_plot:
            plot_distribution(denss, dens_path, code, 'density', 'density', ratio, buckets_cnt, text=text_ratios_map)

        cluster_path = f'distribution/clusters/{f_id}.png'
        if cluster_plot:
            plot_distribution(clusters, cluster_path, code, 'clusters', 'cluster_size', ratio, buckets_cnt, 'green')

        periphery_path = f'distribution/peripheries/{f_id}.png'
        if periphery_plot:
            plot_distribution(peripheries, periphery_path, code, 'peripheries', 'periphery_size', ratio, buckets_cnt,
                              'orange')

        seek_path = f'distribution/seek/{f_id}.png'
        if seek_plot:
            plot_seek_distibution(top_ratio_infos_sorted_by_dens, buckets_cnt, rep, seek_path,
                                  f'seek distribution\n{code}')

    t2 = time.time()
    # print(f'plots took {t2 - t1} sec')
    return infoids


def chrs_per_donor_summary(infoids_sorted_by_dens, ratios):
    cnt = len(infoids_sorted_by_dens)
    ratios_map = {}
    for _ratio in ratios:
        _ratio_infoids = []
        for i in range(int(cnt * _ratio), cnt):
            cur_id = infoids_sorted_by_dens[i]
            _ratio_infoids.append(cur_id)
        chr_per_donor = relation_chromo_chrs_per_donor(_ratio_infoids)
        ratios_map[_ratio] = chr_per_donor
    return ratios_map


def relation_chromo_chrs_per_donor(chromo_infoids):
    donor_chromo_chrs_cnt = {}
    with DonorRepository() as rep:
        for infoid in chromo_infoids:
            donor, chr, f_id = rep.get_by_infoid(infoid)
            donor_chromo_chrs_cnt[donor] = donor_chromo_chrs_cnt.get(donor, 0) + 1
    avg = 0.0
    cnt = 0
    for _, v in donor_chromo_chrs_cnt.items():
        cnt += 1
        avg += v
    return avg / cnt


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
            seek_donor, seek_chr, seek_low, seek_up, label = rep.get_chromo(donor, chr)
            if seek_low and seek_up:
                seek_low, seek_up = bps_to_bins_with_resolution(seek_low, seek_up, resolution)
            print('seek ', donor, f'chr{seek_chr}', seek_low, seek_up)
            print()


def measure_chromos(chromo_clusters):
    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    resolution = cool.info['bin-size']
    accs = []
    recalls = []
    chromo_cnt = 0
    all_cnt = 0
    labels = set()
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

            seek_donor, seek_chr, seek_low, seek_up, label = rep.get_chromo(donor, chr)

            # if label != 'High confidence' and label != 'Linked to high confidence':
            if label != 'High confidence':
                # print('No high confidence chromo')
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
            recall = TP / (TP + FN)
            accs.append(acc)
            recalls.append(recall)
            # print(acc)

    mean_chromo_acc = np.mean(np.array(accs))
    mean_chromo_recall = np.mean(np.array(recalls))
    print(f'mean chromo acc: {mean_chromo_acc * 100}%')
    print(f'mean chromo recall: {mean_chromo_recall * 100}%')
    chromo_marker_acc = chromo_cnt / all_cnt
    print(f'marker acc: {chromo_marker_acc * 100}%')
    return mean_chromo_acc, mean_chromo_recall, chromo_marker_acc


def all_hists(ratio):
    iss = [1, 2, 11, 12, 3]
    # for i in range(1, 4):
    for i in iss:
        # for i in range(1, max_range_pow + 1):
        infoids = hist_patients(i, ratio=ratio, dens_plot=True, seek_plot=True, periphery_plot=False,
                                cluster_plot=False)
        # if infoids:
        #     avg = relation_chromo_chrs_per_donor(infoids)
        #     print(f'chromo chr per 1 donor = {avg} with f_id = {i}')


def seek_test():
    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    resolution = cool.info['bin-size']
    with DonorRepository() as rep:
        cnt = 0
        no_chromo_cnt = 0
        no_data_cnt = 0
        chromo_cnt = 0
        donors = rep.unique_prostate_donors()

        all_patients = 0
        ok_cnt = 0
        for donor in donors:
            all_patients += 1
            ok = True
            for i in range(1, 22):
                cnt += 1
                seek_donor, seek_chr, seek_low, seek_up, label = rep.get_chromo(donor, i)
                print(label)
                if seek_low == -2:
                    no_chromo_cnt += 1
                    continue

                if seek_low == -1:
                    no_data_cnt += 1
                    continue
                if label == 'High confidence':
                    chromo_cnt += 1
                    ok = False
            if ok:
                ok_cnt += 1

    print(f'cnt={cnt}, no_chromo={no_chromo_cnt}, no_data={no_data_cnt}, chromo={chromo_cnt}')
    print(f'all_patients={all_patients}, ok_cnt={ok_cnt}')


def shatter_seek_compare(analyze_donor_chr_pairs=False):
    with DonorRepository() as rep:
        seek_donors = set()
        seek_chr_pairs = set()
        seek_chromo_donors = set()
        seek_chromo_donor_chr_pairs = set()

        rows = rep.get_seek()
        for _donor_seek, chr_seek, label in rows:
            delim = '::'
            seek_donor = _donor_seek[_donor_seek.index(delim) + len(delim):]
            seek_donors.add(seek_donor)
            seek_chr_pairs.add((seek_donor, chr_seek))
            if label == 'High confidence':
                seek_chromo_donors.add(seek_donor)
                seek_chromo_donor_chr_pairs.add((seek_donor, int(chr_seek)))

        pcawg_donors, pcawg_pairs = rep.get_pcawg()

        print(f'SHATTER SEEK CHROMO: {len(seek_chromo_donors) / len(seek_donors)}')
        COMMON_DONORS = seek_donors.intersection(pcawg_donors)
        COMMON_DONOR_CHR_PAIRS = list(map(lambda e: (e[0], int(e[1])), seek_chr_pairs.intersection(pcawg_pairs)))

        iss = [11, 12, 1, 2, 3]
        ratios = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
        print(f'mode donor chr pairs={analyze_donor_chr_pairs}')
        for i in iss:
            print(f'f_id = {i}')
            accs = []
            recalls = []
            precisions = []
            for ratio in ratios:
                infoids = hist_patients(i, ratio=ratio, dens_plot=False, seek_plot=False, periphery_plot=False,
                                        cluster_plot=False)
                chromo_donors = set()
                chromo_donor_chr_pairs = set()
                for infoid in infoids:
                    donor, chr, f_id = rep.get_by_infoid(infoid)
                    chromo_donors.add(donor)
                    chromo_donor_chr_pairs.add((donor, int(chr)))

                seek_v = []
                own_v = []

                if analyze_donor_chr_pairs:
                    for p in COMMON_DONOR_CHR_PAIRS:
                        if p in seek_chromo_donor_chr_pairs:
                            seek_v.append(1)
                        else:
                            seek_v.append(0)

                        if p in chromo_donor_chr_pairs:
                            own_v.append(1)
                        else:
                            own_v.append(0)
                else:
                    for d in COMMON_DONORS:
                        if d in seek_chromo_donors:
                            seek_v.append(1)
                        else:
                            seek_v.append(0)

                        if d in chromo_donors:
                            own_v.append(1)
                        else:
                            own_v.append(0)

                TP, FP, TN, FN = perf_measure(seek_v, own_v)
                acc = (TP + TN) / (TP + TN + FP + FN)
                recall = TP / (TP + FN)
                precision = TP / (TP + FP)
                accs.append(acc)
                recalls.append(recall)
                precisions.append(precision)
                print(f'  percentile_healthy_threshold={ratio}, acc={acc}, recall={recall} precision={precision}')
            f_source = rep.get_proximity_code(i)
            if analyze_donor_chr_pairs:
                save_path = f'seek_compare/by_donor_chr_pairs/{i}.png'
                title = f'seek compare donor-chr pairs\n{f_source}'
            else:
                save_path = f'seek_compare/by_donors/{i}.png'
                title = f'seek compare only donors\n{f_source}'
            plot_seek_compare(accs, recalls, precisions, ratios, title=title, save_path=save_path)


def measure_test():
    # i = 3
    for i in range(1, 4):
        print(f'f_id={i}')
        for ratio in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
            infoids = hist_patients(i, ratio=ratio, dens_plot=False, cluster_plot=False, periphery_plot=False,
                                    seek_plot=False)
            chromo_clusters = find_cromo_clusters(infoids)
            acc, recall, chromo_donorchr_pairs = measure_chromos(chromo_clusters)
            print(f'ratio={ratio}, acc={acc}, recall={recall}, chromo_donorchr_pairs={chromo_donorchr_pairs}')


def identity_hic_oe(x):
    return x


def identity_hic(x):
    return x


def hic_oe_analyzer(oe=False):
    t1 = time.time()

    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    with DonorRepository() as rep:
        if oe:
            f = identity_hic_oe
        else:
            f = identity_hic

        rep.ddl()
        rep.insert_proximity(f)
        f_id = rep.get_proximity_id(inspect.getsource(f))
        donors = rep.unique_prostate_donors()

        corr_path = f'distribution/corr/{f_id}/corr_arrays.npy'
        create_path_if_not_exist(corr_path)

        densess = []
        svscnt = []
        edgescnt = []
        for donor in donors:
            corr_densities, corr_sv_cnt, corr_edges_cnt = analyze_donor(donor=donor, cooler=cool, f_id=f_id, f_proximity=f, rep=rep,
                          hic_plot=False, from_hic=True, oe=oe)
            densess.extend(corr_densities)
            svscnt.extend(corr_sv_cnt)
            edgescnt.extend(corr_edges_cnt)

        with open(corr_path, 'wb') as fl:
            np.save(fl, np.array(densess))
            np.save(fl, np.array(svscnt))
            np.save(fl, np.array(edgescnt))

    t2 = time.time()
    print(f'filling hic oe db took {t2 - t1} sec')


# Pre: analyze donors & collect data before
def corr_test(f_id):
    corr_path = f'distribution/corr/{f_id}/corr_arrays.npy'
    with open(corr_path, 'rb') as fl:
        denss = np.load(fl)
        svs = np.load(fl)
        edgess = np.load(fl)
    print(denss)
    print(svs)
    print(edgess)
    print(len(denss), len(svs), len(edgess))
    print(np.corrcoef(denss, svs)[0, 1])
    print(np.corrcoef(denss, edgess)[0, 1])
    print(np.corrcoef(denss, list(map(math.sqrt, edgess)))[0, 1])


def cluster_sustainability_percentile_test():
    t1 = time.time()
    fid_cluster_map = {}
    ratio_healthy = 0.8
    adj = np.zeros((13, 13))
    with DonorRepository() as rep:
        for f_id in range(1, 13):
            infoids = hist_patients(f_id, ratio=ratio_healthy, dens_plot=False, cluster_plot=False, periphery_plot=False,
                                    seek_plot=False)
            for info_id in infoids:
                cluster_set = set(rep.get_cluster(info_id))
                cur_set: set = fid_cluster_map.get(f_id, set())
                fid_cluster_map[f_id] = cur_set.union(cluster_set)
    for i in range(1, 13):
        for j in range(i, 13):
            adj[i][j] = len(fid_cluster_map[i].intersection(fid_cluster_map[j]))

    print(adj)

    plt.title(f'percentile healthy = {ratio_healthy}')
    plt.imshow(adj, interpolation='none')

    for (j, i), label in np.ndenumerate(adj):
        plt.text(i, j, int(label), ha='center', va='center')

    plt.show()

    t2 = time.time()
    print(f'sustainability test took {t2 - t1} sec')


if __name__ == '__main__':
    # main()

    # all_hists(ratio=0.8)  # 0.9965: -1,  0.975: -2

    # seek_test()

    # measure_test()

    shatter_seek_compare(analyze_donor_chr_pairs=False)
    # shatter_seek_compare(analyze_donor_chr_pairs=True)


    # hic_oe_analyzer(oe=True)
    #corr_test(4)
    # cluster_sustainability_percentile_test()


