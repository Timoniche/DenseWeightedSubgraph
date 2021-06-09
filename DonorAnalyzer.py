import inspect
import math
import time
from enum import Enum

import cooler
import numpy as np
import os
import sys

from DenseUtils import heatmap_with_breakpoints_and_cluster, create_path_if_not_exist, plot_distribution, perf_measure, \
    plot_seek_distibution, plot_seek_compare, plot_sustainability, div_with_none, filter_nones, plot_mixed_denss, \
    find_first_predicate_index
from DonorRepository import DonorRepository
from DonorService import collect_chr_bins_map_with_resolution, bps_to_bins_with_resolution, filter_svs_with_resolution, \
    filter_chr_svs_with_resolution
from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density

import matplotlib.pyplot as plt

from HiCUtils import zeros_to_nan, normalize_intra
from functions import functions
from generate_functions import generate_functions, max_range_pow
from tqdm import tqdm

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
donorspath = script_dir + '/donors'


def analyze_donor(donor, cooler, f_id, f_proximity, rep: DonorRepository, hic_plot, cluster_threshold_size,
                  from_hic=False, oe=False,
                  in_bp=False, _table='sv_intra'):
    print(donor)
    svs = rep.get_donor_sv_chr_1_21(donor, _table)

    resolution = cooler.info['bin-size']

    if not in_bp:
        chr_sv_map = collect_chr_bins_map_with_resolution(svs, resolution)
    else:
        chr_sv_map = filter_svs_with_resolution(svs, resolution)

    corr_densities = []
    corr_sv_cnt = []
    corr_edges_cnt = []


    # filtered svs should present in the database
    ALL_DONOR_CHRS = rep.get_donor_chrs_1_21(donor, _table)
    ALL_DONOR_CHRS = list(map(lambda row: row[0], ALL_DONOR_CHRS))
    for _chr in ALL_DONOR_CHRS:
        if not chr_sv_map.get(_chr, []):
            rep.insert_donorinfo(donor, _chr, f_id)
            info_id = rep.get_info_id(donor, _chr, f_id)
            rep.insert_dense(info_id, 0)

    for (chr_n, coord_pairs) in chr_sv_map.items():

        all_coords = set()
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

        for (x, y) in coord_pairs:
            all_coords.add(x)
            all_coords.add(y)
        all_coords = list(all_coords)

        if not in_bp:
            all_coords_bin = all_coords.copy()
            coord_pairs_bin = coord_pairs.copy()
        else:
            all_coords_bin = list(set(map(lambda e: int(e / resolution), all_coords)))
            coord_pairs_bin = list(set(map(lambda p: bps_to_bins_with_resolution(p[0], p[1], resolution), coord_pairs)))

        dirpath = donorspath + f'/{donor}/{chr_n}'
        create_path_if_not_exist(dirpath)

        filepath = donorspath + f'/{donor}/{chr_n}'

        with open(filepath, 'w') as outfile:
            for i_idx in range(len(all_coords_bin)):
                for j_idx in range(i_idx + 1, len(all_coords_bin)):
                    i = all_coords_bin[i_idx]
                    j = all_coords_bin[j_idx]
                    number_of_nodes = max(number_of_nodes, max(i, j))
                    dist = dist_mat[i][j]
                    if dist != np.nan:
                        close_f = f_proximity(dist)
                        if close_f == 0:
                            continue
                        number_of_edges += 1
                        outfile.write(f'{i} {j} {close_f}\n')

        if number_of_edges == 0:
            rep.insert_donorinfo(donor, chr_n, f_id)
            info_id = rep.get_info_id(donor, chr_n, f_id)
            rep.insert_dense(info_id, 0)
            continue

        some_delta_just_for_sure = 5
        cluster_bins = WFind_Densest_Subgraph(number_of_nodes + some_delta_just_for_sure, number_of_edges, filepath,
                                              cluster_threshold_size)
        # assert cluster_bins != [], 'not enough accuracy in find densest subgraph algorithm'
        if not cluster_bins:
            print('ZERO CLUSTER')
        if hic_plot:
            heatmap_with_breakpoints_and_cluster(mat,
                                                 f'normed hic & breakpoints chr{chr_n}\n{inspect.getsource(f_proximity)}',
                                                 coord_pairs_bin,
                                                 cluster_bins,
                                                 save_path=donorspath + f'/{donor}/{chr_n}.png')
        rep.insert_donorinfo(donor, chr_n, f_id)
        dens = WFind_Density(cluster_bins, filepath)
        print(dens)

        corr_densities.append(dens)
        corr_sv_cnt.append(len(coord_pairs_bin))
        corr_edges_cnt.append(number_of_edges)

        info_id = rep.get_info_id(donor, chr_n, f_id)
        rep.insert_dense(info_id, dens)
        collect_periphery_and_cluster(rep, info_id, in_bp, cluster_bins, all_coords, resolution, coord_pairs)

    return corr_densities, corr_sv_cnt, corr_edges_cnt


def collect_periphery_and_cluster(rep, info_id, in_bp, cluster_bins, all_coords, resolution, coord_pairs):
    periphery = set()
    cluster_bins_set = set(cluster_bins)
    if not in_bp:
        rep.insert_cluster(info_id, tuple(cluster_bins))
    else:
        coords_bp_to_insert = []
        for coord in all_coords:
            coord_bin = int(coord / resolution)
            if coord_bin in cluster_bins_set:
                coords_bp_to_insert.append(coord)
        rep.insert_cluster(info_id, tuple(coords_bp_to_insert))

    if not in_bp:
        for b in coord_pairs:
            i = b[0]
            j = b[1]
            if i in cluster_bins and j not in cluster_bins:
                periphery.add(j)
            if j in cluster_bins and i not in cluster_bins:
                periphery.add(i)
    else:
        for b in coord_pairs:
            i = int(b[0] / resolution)
            j = int(b[1] / resolution)
            if i in cluster_bins and j not in cluster_bins:
                periphery.add(b[1])
            if j in cluster_bins and i not in cluster_bins:
                periphery.add(b[0])
    rep.insert_periphery(info_id, tuple(list(periphery)))

    print(f'clusters {cluster_bins}')
    print(f'periphery: {periphery}')
    print(f'all sv: {all_coords}')


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

    cool = cooler.Cooler('healthy_hics/new_cool_480kb.cool')
    with DonorRepository() as rep:
        rep.ddl()
        generate_functions()
        for f in functions:
            rep.insert_proximity(f)

        donors = rep.unique_prostate_donors()
        f_ids = [i for i in range(1, max_range_pow + 1)]

        for f_id in f_ids:
            densess = []
            svscnt = []
            edgescnt = []
            corr_path = f'distribution/corr/{f_id}/corr_arrays.npy'
            create_path_if_not_exist(corr_path)
            for donor in tqdm(donors):
                corr_densities, corr_sv_cnt, corr_edges_cnt = analyze_donor(donor=donor, cooler=cool, f_id=f_id,
                                                                            f_proximity=functions[f_id - 1], rep=rep,
                                                                            hic_plot=False, in_bp=True)
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


def seek_hills_distrib(ratio, fid=16):
    infoids_sorted_by_dens, infoids, thresholds = hist_patients(f_id=fid, ratio=0, dens_plot=False, seek_plot=False,
                                                                periphery_plot=False,
                                                                cluster_plot=False)
    cnt_left = 0
    cnt_right = 0
    sz = len(infoids)
    sz_left = sz * ratio
    sz_right = sz * (1 - ratio)
    print(f'sz {sz} sz_left {sz_left} sz_right {sz_right}')
    with DonorRepository() as rep:
        i = 0
        for info_id in infoids:
            donor, chrN, _ = rep.get_by_infoid(info_id)
            _, _, _, _, label = rep.get_chromo(donor, chrN)
            if label == 'High confidence':
                if i / sz >= ratio:
                    cnt_right += 1
                else:
                    cnt_left += 1
            i += 1
    print(f'concentration: left {cnt_left / sz_left} right {cnt_right / sz_right}')
    return cnt_left, cnt_right


def hist_patients(f_id, ratio, dens_plot=True, cluster_plot=True, periphery_plot=True, seek_plot=True, chr_num=-1,
                  weights=None, cluster_threshold_sz=0, dens_threshold=0.0, types5plot=False, _table='sv_intra'):
    if weights is None:
        weights = []
    t1 = time.time()
    denss_info = []
    clusters = []
    peripheries = []
    infoids_after_ratio = []
    with DonorRepository() as rep:
        donors = rep.unique_prostate_donors(_table)
        for donor in donors:
            if chr_num == -1:
                rng = range(1, 22)
            else:
                rng = range(chr_num, chr_num + 1)
            for i in rng:
                info_id = rep.get_info_id(donor, i, f_id)
                if info_id != -1:
                    dens = rep.get_density(info_id)
                    if weights:
                        dens *= weights[i - 1]
                    cluster = len(rep.get_cluster(info_id))
                    periphery = len(rep.get_periphery(info_id))

                    if cluster >= cluster_threshold_sz and dens >= dens_threshold:
                        denss_info.append((dens, info_id))
                        clusters.append(cluster)
                        peripheries.append(periphery)

        code = rep.get_proximity_code(f_id)

        denss = list(map(lambda x: x[0], denss_info))
        print(f'f_id={f_id}')
        print(f'avg density: {np.average(denss)}')
        print(f'high quantile density: {np.quantile(denss, 0.75)}')
        print(f'avg cluster: {np.average(clusters)}')
        print(f'high quantile cluster: {np.quantile(clusters, 0.75)}')
        print(f'avg periphery: {np.average(peripheries)}')
        print(f'high quantile periphery: {np.quantile(peripheries, 0.75)}')
        print(f'cnt: {len(denss)}')

        top_ratio_infos_sorted_by_dens = sorted(denss_info, key=lambda p: p[0])
        infoids_sorted_by_dens = list(map(lambda p: p[1], top_ratio_infos_sorted_by_dens))
        cnt = len(denss)

        for i in range(int(cnt * ratio), cnt):
            cur_id = top_ratio_infos_sorted_by_dens[i][1]
            infoids_after_ratio.append(cur_id)

        ratios = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
        ratios_map = chrs_per_donor_summary(infoids_sorted_by_dens, ratios)
        text_ratios_map = ''
        for k, v in ratios_map.items():
            text_ratios_map += f'ratio={k} : ' + 'chr_per_donor=%0.2f\n' % v

        buckets_cnt = 100
        dens_path = f'distribution/densities/{f_id}.png'
        if dens_plot:
            plot_distribution(denss, dens_path, code, 'density', 'density', ratio, buckets_cnt, text=text_ratios_map)
        thresholds_paths = f'distribution/thresholds/{f_id}.png'
        thresholds = plot_mixed_denss(sorted(denss), to_plot=types5plot, save_path=thresholds_paths)

        cluster_path = f'distribution/clusters/{f_id}.png'
        if cluster_plot:
            plot_distribution(clusters, cluster_path, code, 'clusters', 'cluster_size', ratio, buckets_cnt,
                              color='green')

        periphery_path = f'distribution/peripheries/{f_id}.png'
        if periphery_plot:
            plot_distribution(peripheries, periphery_path, code, 'peripheries', 'periphery_size', ratio, buckets_cnt,
                              color='orange')

        seek_path = f'distribution/seek/{f_id}.png'
        if seek_plot:
            plot_seek_distibution(top_ratio_infos_sorted_by_dens, buckets_cnt, rep, seek_path,
                                  f'seek distribution\n{code}')

    t2 = time.time()
    # print(f'plots took {t2 - t1} sec')
    return infoids_sorted_by_dens, sorted(denss), infoids_after_ratio, thresholds


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

            if label != 'High confidence':
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
    iss = [3, 4]
    # for i in range(1, 4):
    threshold_hic_dens = 200
    # threshold_hic_dens = 2.7495356
    for i in iss:
        all_infoids_by_denss, sorted_denss, infoids_after_ratios, thresholds = hist_patients(i, ratio=ratio,
                                                                                             dens_plot=True,
                                                                                             seek_plot=False,
                                                                                             periphery_plot=False,
                                                                                             cluster_plot=False,
                                                                                             cluster_threshold_sz=0,
                                                                                             dens_threshold=0.000001,
                                                                                             types5plot=True)
        densss = sorted_denss[int(ratio * len(sorted_denss))]
        print(densss)
        index_ = find_first_predicate_index(sorted_denss, lambda d: d > threshold_hic_dens)
        print(f'cnt chromo ids: {len(all_infoids_by_denss) - index_}')


def cnt_chr_weights(fid):
    weights = []
    with DonorRepository() as rep:
        for _chr in range(1, 22):
            infoids = hist_patients(f_id=fid, ratio=0, dens_plot=False, seek_plot=False, periphery_plot=False,
                                    cluster_plot=False, chr_num=_chr)
            denss = []
            for info_id in infoids:
                denss.append(rep.get_density(info_id))
            avg_chr_dense = np.average(denss)
            weights.append(1.0 / avg_chr_dense)
    return weights


def weighted_hist(fid, ratio):
    # weights = cnt_chr_weights(fid)
    infoids = hist_patients(f_id=fid, ratio=ratio, dens_plot=True, seek_plot=False, periphery_plot=False,
                            cluster_plot=False)


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


class SeekClassification(str, Enum):
    NO = 'No'
    LOW = 'Low confidence'
    LINKED_TO_LOW = 'Linked to low confidence'
    LINKED_TO_HIGH = 'Linked to high confidence'
    HIGH = 'High confidence'

    @staticmethod
    def get_range(label, thresholds):
        if label == SeekClassification.NO:
            return range(thresholds[0], thresholds[1])
        elif label == SeekClassification.LOW:
            return range(thresholds[1], thresholds[2])
        elif label == SeekClassification.LINKED_TO_LOW:
            return range(thresholds[2], thresholds[3])
        elif label == SeekClassification.LINKED_TO_HIGH:
            return range(thresholds[3], thresholds[4])
        elif label == SeekClassification.HIGH:
            return range(thresholds[4], thresholds[5])


def gmm_seek_5types_compare(analyze_donor_chr_pairs, seekLabel):
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
            if label == seekLabel:
                seek_chromo_donors.add(seek_donor)
                seek_chromo_donor_chr_pairs.add((seek_donor, int(chr_seek)))

        pcawg_donors, pcawg_pairs = rep.get_pcawg()

        print(f'SHATTER SEEK CHROMO: {len(seek_chromo_donors) / len(seek_donors)}')
        COMMON_DONORS = seek_donors.intersection(pcawg_donors)
        COMMON_DONOR_CHR_PAIRS = list(map(lambda e: (e[0], int(e[1])), seek_chr_pairs.intersection(pcawg_pairs)))

        iss = [11]
        print(f'mode donor chr pairs={analyze_donor_chr_pairs}')
        for i in iss:
            print(f'f_id = {i}')
            accs = []
            recalls = []
            precisions = []

            common_chromo_donor_pairs = []

            all_info_ids_sorted_by_denss, infoids, thresholds = hist_patients(i, ratio=1, dens_plot=False,
                                                                              seek_plot=False,
                                                                              periphery_plot=False,
                                                                              cluster_plot=False,
                                                                              cluster_threshold_sz=3)
            chromo_donors = set()
            chromo_donor_chr_pairs = set()

            for i in SeekClassification.get_range(label, thresholds):
                donor, chr, f_id = rep.get_by_infoid(all_info_ids_sorted_by_denss[i])
                chromo_donors.add(donor)
                chromo_donor_chr_pairs.add((donor, int(chr)))
            seek_v = []
            own_v = []

            if analyze_donor_chr_pairs:
                donors_pairs_selection = chromo_donor_chr_pairs.intersection(COMMON_DONOR_CHR_PAIRS)
                print(f'PAIRS_SELECTION_SIZE = {len(donors_pairs_selection)}')
                for p in chromo_donor_chr_pairs.intersection(COMMON_DONOR_CHR_PAIRS):
                    if p in seek_chromo_donor_chr_pairs:
                        seek_v.append(1)
                    else:
                        seek_v.append(0)

                    if p in chromo_donor_chr_pairs:
                        own_v.append(1)
                    else:
                        own_v.append(0)
                    if seek_v[-1] == 1 and own_v[-1] == 1:
                        common_chromo_donor_pairs.append(p)
            else:
                donors_selection = chromo_donors.intersection(COMMON_DONORS)
                print(f'DONOR_SELECTION_SIZE = {len(donors_selection)}')
                for d in donors_selection:
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
            print(f'  label: {seekLabel}, acc={acc}, recall={recall} precision={precision}')
            print(f'  common chromo donor pairs:\n{common_chromo_donor_pairs}')


def shatter_seek_compare(analyze_donor_chr_pairs, seekLabel):
    cl_szs_seek = []
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
            if label == seekLabel:
                seek_chromo_donors.add(seek_donor)
                seek_chromo_donor_chr_pairs.add((seek_donor, int(chr_seek)))

                iii = rep.get_info_id(seek_donor, chr_seek, f_id=13)
                cl_sz = len(rep.get_cluster(iii))
                cl_szs_seek.append(cl_sz)

        pcawg_donors, pcawg_pairs = rep.get_pcawg()

        print(f'SHATTER SEEK CHROMO: {len(seek_chromo_donors) / len(seek_donors)}')
        COMMON_DONORS = seek_donors.intersection(pcawg_donors)
        COMMON_DONOR_CHR_PAIRS = list(map(lambda e: (e[0], int(e[1])), seek_chr_pairs.intersection(pcawg_pairs)))

        # iss = [11, 12, 1, 2, 3]
        iss = [11]
        ratios = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
        # ratios = [0.75]
        print(f'mode donor chr pairs={analyze_donor_chr_pairs}')
        for i in iss:
            print(f'f_id = {i}')
            accs = []
            recalls = []
            precisions = []
            for ratio in ratios:
                common_chromo_donor_pairs = []

                sorted_denss, infoids, thresholds = hist_patients(i, ratio=ratio, dens_plot=False, seek_plot=False,
                                                                  periphery_plot=False,
                                                                  cluster_plot=False, cluster_threshold_sz=3)
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
                        if seek_v[-1] == 1 and own_v[-1] == 1:
                            common_chromo_donor_pairs.append(p)
                    common_chromo_pairs_path = f'seek_compare/{f_id}/{ratio}/common_pairs.npy'
                    create_path_if_not_exist(common_chromo_pairs_path)
                    with open(common_chromo_pairs_path, 'wb') as fl:
                        np.save(fl, np.array(common_chromo_donor_pairs))
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
                print(f'  common chromo donor pairs:\n{common_chromo_donor_pairs}')
            f_source = rep.get_proximity_code(i)
            if analyze_donor_chr_pairs:
                save_path = f'seek_compare/by_donor_chr_pairs/{i}.png'
                title = f'seek compare donor-chr pairs\n{f_source}'
            else:
                save_path = f'seek_compare/by_donors/{i}.png'
                title = f'seek compare only donors\n{f_source}'
            plot_seek_compare(accs, recalls, precisions, ratios, title=title, save_path=save_path)
    med = np.quantile(cl_szs_seek, 0.5)
    q75 = np.quantile(cl_szs_seek, 0.75)
    avg = np.average(cl_szs_seek)
    print(f'CHROMO SHATTER SEEK median: {med}, high quantile: {q75}, avg: {avg}')


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


def identity_hic_40kb(x):
    return x


def identity_hic(x):
    return x


def frankenstein_hic(x):
    return x


def frankenstein_hic_oe_480kb(x):
    return x


def frankenstein_hic_40kb(x):
    return x


def filter_identity_hic_480kb(x):
    return x


def filter_identity_hic_40kb(x):
    return x


def hic_oe_analyzer(cool, f, _table, oe=False):
    t1 = time.time()
    cluster_threshold_size = 0
    with DonorRepository() as rep:
        rep.ddl()
        rep.insert_proximity(f)
        f_id = rep.get_proximity_id(inspect.getsource(f))
        donors = rep.unique_prostate_donors(_table)

        corr_path = f'distribution/corr/{f_id}/corr_arrays.npy'
        create_path_if_not_exist(corr_path)

        densess = []
        svscnt = []
        edgescnt = []
        for donor in tqdm(donors):
            corr_densities, corr_sv_cnt, corr_edges_cnt = analyze_donor(donor=donor, cooler=cool, f_id=f_id,
                                                                        f_proximity=f, rep=rep,
                                                                        hic_plot=False, from_hic=True, oe=oe,
                                                                        in_bp=True,
                                                                        cluster_threshold_size=cluster_threshold_size,
                                                                        _table=_table)
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


def diff_example_own_seek(f_id, ratio):
    print(f'f_id = {f_id}, ratio = {ratio}')
    common_chromo_pairs_path = f'seek_compare/{f_id}/{ratio}/common_pairs.npy'
    with open(common_chromo_pairs_path, 'rb') as fl:
        donor_pairs: np.ndarray = np.load(fl)
        print(f'donor_pairs:\n{donor_pairs}')
        donors = list(set(map(lambda p: p[0], donor_pairs.tolist())))
        print(f'donors:\n{donors}')

        cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
        with DonorRepository() as rep:
            for donor in tqdm(donors):
                _, _, _ = analyze_donor(donor=donor, cooler=cool, f_id=f_id,
                                        f_proximity=identity_hic, rep=rep,
                                        hic_plot=True, from_hic=True, oe=False,
                                        in_bp=True)


def find_the_most_diff_seek_pair(f_id, ratio):
    print(f'f_id = {f_id}, ratio = {ratio}')
    common_chromo_pairs_path = f'seek_compare/{f_id}/{ratio}/common_pairs.npy'
    with open(common_chromo_pairs_path, 'rb') as fl:
        donor_pairs: np.ndarray = np.load(fl)
        print(f'donor_pairs:\n{donor_pairs}')
        for (donor, chr) in donor_pairs:
            min_own, max_own, min_seek, max_seek = donor_chr_pair_analyzer(donor, chr, f_id)
            if min_own < min_seek or max_own > max_seek:
                print(f'donor: {donor}, chr: {chr}')
                print(min_own, max_own, min_seek, max_seek)


def donor_chr_pair_analyzer(donor, chr, f_id):
    print(f'donor: {donor}, chr: {chr}')
    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    resolution = cool.info['bin-size']
    with DonorRepository() as rep:
        info_id = rep.get_info_id(donor, chr, f_id)
        # print(f'info_id: {info_id}')
        _, _, start_bp, end_bp, _ = rep.get_seek_markup_by_donor_chr(donor, chr)
        # print(f'start_bp: {start_bp}, end_bp: {end_bp}')
        _svs = rep.get_donor_chr_svs(donor, chr)
        filtered_svs = filter_chr_svs_with_resolution(_svs, resolution)
        # print(f'filtered_svs: {filtered_svs}')
        chromo_seek_bps = seek_chromo_svs_by_range(start_bp, end_bp, filtered_svs)
        chromo_seek_bins = list(set(map(lambda x: int(x / resolution), chromo_seek_bps)))

        print(f'cnt bp seek: {len(chromo_seek_bps)}')
        print(f'chromo seek bps: {chromo_seek_bps}')
        print(f'cnt bin seek: {len(chromo_seek_bins)}')
        print(f'chromo seek bins: {chromo_seek_bins}')

        chromo_own_bps = rep.get_cluster(info_id)
        chromo_own_bins = list(set(map(lambda x: int(x / resolution), chromo_own_bps)))

        print(f'cnt bp own: {len(chromo_own_bps)}')
        print(f'chromo own bps: {chromo_own_bps}')
        print(f'cnt bin own: {len(chromo_own_bins)}')
        print(f'chromo own bins: {chromo_own_bins}')

        set_seek_bins = set(chromo_seek_bins)
        set_own_bins = set(chromo_own_bins)

        common_bins = set_seek_bins.intersection(set_own_bins)
        extra_seek_bins = set_seek_bins.difference(set_own_bins)
        extra_own_bins = set_own_bins.difference(set_seek_bins)
        print(f'common_bins: {common_bins}')
        print(f'extra seek bins: {extra_seek_bins}')
        print(f'extra own bins: {extra_own_bins}')

        return min(chromo_own_bps), max(chromo_own_bps), min(chromo_seek_bps), max(chromo_seek_bps),


def seek_chromo_svs_by_range(_start, _end, svs):
    chromo_bps = set()
    for bp1, bp2 in svs:
        if _start <= bp1 <= _end:
            chromo_bps.add(bp1)
        if _start <= bp2 <= _end:
            chromo_bps.add(bp2)
    return list(chromo_bps)


def cluster_sustainability_percentile_test():
    t1 = time.time()
    fid_cluster_map = {}
    fid_donor_chr_pairs = {}
    fid_chromo_donors = {}
    ratio_healthy = 0.75
    functions_cnt = 12
    adj_cluster = np.zeros((functions_cnt, functions_cnt))
    adj_donor_chr = adj_cluster.copy()
    adj_donor = adj_cluster.copy()
    with DonorRepository() as rep:
        for f_id in range(1, functions_cnt + 1):
            infoids = hist_patients(f_id, ratio=ratio_healthy, dens_plot=False, cluster_plot=False,
                                    periphery_plot=False,
                                    seek_plot=False)
            for info_id in infoids:
                cluster_set = set(rep.get_cluster(info_id))
                cur_set: set = fid_cluster_map.get(f_id - 1, set())
                fid_cluster_map[f_id - 1] = cur_set.union(cluster_set)

                _donor, _chr, _ = rep.get_by_infoid(info_id)
                cur_donor_chr_set: set = fid_donor_chr_pairs.get(f_id - 1, set())
                cur_donor_chr_set.add((_donor, _chr))
                fid_donor_chr_pairs[f_id - 1] = cur_donor_chr_set

                cur_donor_set: set = fid_chromo_donors.get(f_id - 1, set())
                cur_donor_set.add(_donor)
                fid_chromo_donors[f_id - 1] = cur_donor_set

    for i in range(functions_cnt):
        for j in range(i, functions_cnt):
            adj_cluster[i][j] = len(fid_cluster_map[i].intersection(fid_cluster_map[j]))
            adj_donor_chr[i][j] = len(fid_donor_chr_pairs[i].intersection(fid_donor_chr_pairs[j]))
            adj_donor[i][j] = len(fid_chromo_donors[i].intersection(fid_chromo_donors[j]))

    cluster_title = f'common cluster breakpoints\npercentile healthy = {ratio_healthy}'
    donor_chr_title = f'common chromo donor-chr pairs\npercentile healthy = {ratio_healthy}'
    donor_title = f'common chromo donors\npercentile healthy = {ratio_healthy}'
    cluster_save_path = 'seek_compare/common_cluster_bps.png'
    donor_chr_save_path = 'seek_compare/common_donor_chr_pairs.png'
    donor_save_path = 'seek_compare/common_donors.png'

    plot_sustainability(adj=adj_cluster, functions_cnt=functions_cnt, title=cluster_title, save_path=cluster_save_path)
    plot_sustainability(adj=adj_donor_chr, functions_cnt=functions_cnt, title=donor_chr_title,
                        save_path=donor_chr_save_path)
    plot_sustainability(adj=adj_donor, functions_cnt=functions_cnt, title=donor_title, save_path=donor_save_path)

    t2 = time.time()
    print(f'sustainability test took {t2 - t1} sec')


def inner_chromo_donor_metrics(fid=13, ratio_healthy=0.7):
    with DonorRepository() as rep:
        infoids = hist_patients(f_id=fid, ratio=ratio_healthy, dens_plot=False, cluster_plot=False,
                                periphery_plot=False,
                                seek_plot=False)
        chromo_donors = set()
        for infoid in infoids:
            _donor, _chr, _ = rep.get_by_infoid(infoid)
            chromo_donors.add(_donor)

        count_avg_hill_metrics(chromo_donors, fid)

        all_donors = set(rep.unique_prostate_donors())
        non_chromo_donors = all_donors.difference(chromo_donors)

        count_avg_hill_metrics(non_chromo_donors, fid)


def count_avg_hill_metrics(donors, fid):
    tops_to_next, tops_to_avg, tops_to_quantile, tops_to_med = donor_hill_metrics(donors, fid)
    avg_tops_to_next = np.average(filter_nones(tops_to_next))
    avg_tops_to_avg = np.average(filter_nones(tops_to_avg))
    avg_tops_to_quantile = np.average(filter_nones(tops_to_quantile))
    avg_tops_to_med = np.average(filter_nones(tops_to_med))
    print(f'avg top to: next {avg_tops_to_next}, avg {avg_tops_to_avg}, ' +
          f'high quantile {avg_tops_to_quantile}, median {avg_tops_to_med}')


def donor_hill_metrics(donors, fid):
    tops_to_next = []
    tops_to_avg = []
    tops_to_quantile = []
    tops_to_med = []
    with DonorRepository() as rep:
        for chromo_donor in donors:
            # print(chromo_donor)

            donor_all_chrs_infoids = rep.get_donor_fid_infoids(chromo_donor, fid)
            donor_all_chrs_denss = list(map(rep.get_density, donor_all_chrs_infoids))
            sz = len(donor_all_chrs_denss)
            all_chrs = 21
            for i in range(all_chrs - sz):
                donor_all_chrs_denss.append(0)
            denss_sorted = sorted(donor_all_chrs_denss)
            # print(denss_sorted)

            top_hill = denss_sorted[-1]
            next_hill = denss_sorted[-2]
            avg_dens = np.average(denss_sorted)
            quantile_dens = np.quantile(denss_sorted, 0.75)
            med_dens = np.quantile(denss_sorted, 0.5)

            top_to_next = div_with_none(top_hill, next_hill)
            top_to_avg = div_with_none(top_hill, avg_dens)
            top_to_quantile = div_with_none(top_hill, quantile_dens)
            top_to_med = div_with_none(top_hill, med_dens)

            tops_to_next.append(top_to_next)
            tops_to_avg.append(top_to_avg)
            tops_to_quantile.append(top_to_quantile)
            tops_to_med.append(top_to_med)
            # print(f'top to: next {top_to_next}, avg {top_to_avg}, high quantile {top_to_quantile}, median {top_to_med}')
    cnts_in_bucket, bins = np.histogram(filter_nones(np.array(tops_to_next)), bins=100)
    plt.hist(filter_nones(np.array(tops_to_next)), bins)
    plt.show()
    plt.hist(filter_nones(np.array(tops_to_avg)), bins)
    plt.show()
    return tops_to_next, tops_to_avg, tops_to_quantile, tops_to_med


def cluster_size_test(ratio):
    fid = 13
    infoids = hist_patients(fid, ratio=ratio, dens_plot=False, seek_plot=False, periphery_plot=False,
                            cluster_plot=False)
    infoids_all = hist_patients(fid, ratio=0, dens_plot=False, seek_plot=False, periphery_plot=False,
                                cluster_plot=False)

    cluster_szs_chromo = []
    szs = []
    with DonorRepository() as rep:
        for infoid in infoids:
            cluster_szs_chromo.append(len(rep.get_cluster(infoid)))
        for infoid in infoids_all:
            if infoid not in infoids:
                szs.append(len(rep.get_cluster(infoid)))
    med = np.quantile(cluster_szs_chromo, 0.5)
    q75 = np.quantile(cluster_szs_chromo, 0.75)
    avg = np.average(cluster_szs_chromo)
    print(f'CHROMO median: {med}, high quantile: {q75}, avg: {avg}')
    med = np.quantile(szs, 0.5)
    q75 = np.quantile(szs, 0.75)
    avg = np.average(szs)
    print(f'NO chromo median: {med}, high quantile: {q75}, avg: {avg}')


def test_cnt_ids():
    #i13 = 13
    #i17f = 17
    i13 = 11
    i17f = 16
    all_infoids_by_denss13, sorted_denss13, infoids_after_ratios13, thresholds13 = hist_patients(i13, ratio=1,
                                                                                         dens_plot=True,
                                                                                         seek_plot=False,
                                                                                         periphery_plot=False,
                                                                                         cluster_plot=False,
                                                                                         cluster_threshold_sz=0,
                                                                                         dens_threshold=0,
                                                                                         types5plot=False)
    all_infoids_by_denss17, sorted_denss17, infoids_after_ratios17, thresholds17 = hist_patients(i17f, ratio=1,
                                                                                         dens_plot=True,
                                                                                         seek_plot=False,
                                                                                         periphery_plot=False,
                                                                                         cluster_plot=False,
                                                                                         cluster_threshold_sz=0,
                                                                                         dens_threshold=0,
                                                                                         types5plot=False)
    with DonorRepository() as rep:
        set13 = set()
        set17f = set()
        for id13 in all_infoids_by_denss13:
            donor, chrN, fid = rep.get_by_infoid(id13)
            set13.add((donor, chrN))
        for id17 in all_infoids_by_denss17:
            donor, chrN, fid = rep.get_by_infoid(id17)
            set17f.add((donor, chrN))

    set_diff13 = set13.difference(set17f)
    set_diff17 = set17f.difference(set13)

    print('diff 13')
    print(set_diff13)

    print('diff 17f')
    print(set_diff17)

    with DonorRepository() as rep:
        id_1 = rep.get_info_id('0082_CRUK_PC_0082', 13, 13)
        id_2 = rep.get_info_id('0082_CRUK_PC_0082', 13, 17)
        print(rep.get_cluster(id_1))
        print(id_1, id_2)


class Tables(str, Enum):
    SIMULATED = 'frankenstein',
    REAL = 'sv_intra'


def frankenstein_hic_480kb_TEST2(x):
    return x

if __name__ == '__main__':
    # main()
    #test_cnt_ids()
    all_hists(ratio=0.7)
    # cluster_size_test(ratio=0.71)

    # seek_test()

    # measure_test()

    # shatter_seek_compare(analyze_donor_chr_pairs=False, seekLabel=SeekClassification.HIGH)

    # shatter_seek_compare(analyze_donor_chr_pairs=True)

    # gmm_seek_5types_compare(analyze_donor_chr_pairs=False, seekLabel=SeekClassification.NO)
    # diff_example_own_seek(11, 0.75)

    # donor_chr_pair_analyzer('0091_CRUK_PC_0091', 13, 11)

    # inner_chromo_donor_metrics()

    # weighted_hist(11, 0.7)

    # l, r = seek_hills_distrib(0.71, fid=11)
    # print(l, r)

    # find_the_most_diff_seek_pair(11, 0.75)

    #cool = cooler.Cooler('healthy_hics/new_cool_480kb.cool')
    #cool = cooler.Cooler('healthy_hics/GSM3564252_RWPE1_HiC_40k.normalized.matrix.cool')
    #hic_oe_analyzer(cool=cool, f=identity_hic_40kb, _table=Tables.REAL, oe=False)
    #hic_oe_analyzer(cool=cool, f=frankenstein_hic_480kb_TEST2, _table=Tables.SIMULATED, oe=False)

    # corr_test(4)

    # cluster_sustainability_percentile_test()
