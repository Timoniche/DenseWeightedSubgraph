import errno
import os
import statistics

import matplotlib.pyplot as plt
import cooler
import pandas as pd
import numpy as np
import inspect

from matplotlib import rcParams

from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density
from HiCUtils import zeros_to_nan, normalize_intra


def all_oe(c):
    all_df = []
    for current_chr in c.chromnames:
        mat = c.matrix(balance=False).fetch(current_chr)
        mat_nan = zeros_to_nan(mat)
        mat_norm = normalize_intra(mat_nan)
        heatmap(mat_norm, f'normed_hic {current_chr}')
        df = pd.DataFrame(mat_norm)
        df = df.stack(dropna=True)
        df = df.to_frame()
        df.reset_index(inplace=True)
        df.columns = ['bin1', 'bin2', 'HiC_score']
        df['chrom'] = current_chr

        all_df.append(df)

    return pd.concat(all_df, ignore_index=True)


def find_first_predicate_index(arr, predicate):
    for i in range(len(arr)):
        if predicate(arr[i]):
            return i
    return -1


def plot_distribution(arr, store_path, function_code, label, xlabel, color='null', ratio=0.975):
    rcParams.update({'figure.autolayout': True})
    buckets_cnt = 100
    if color != 'null':
        cnts_in_bucket, bins, patches = plt.hist(arr, bins=buckets_cnt, color=color)
    else:
        cnts_in_bucket, bins, patches = plt.hist(arr, bins=buckets_cnt)
    cnts_sum = np.cumsum(cnts_in_bucket)
    index_ratio_from = find_first_predicate_index(cnts_sum, lambda x: x >= len(arr) * ratio)
    assert index_ratio_from != -1, 'density in the ratio percentile must exist'
    bars = len(patches)
    for i in range(index_ratio_from, bars):
        patches[i].set_facecolor('r')
    dens_med = statistics.median(arr)
    dens_quantile = np.quantile(arr, 0.75)
    plt.title(f'{label}, {function_code}' +
              'median=%0.2f\n' % dens_med +
              'quantile=%0.2f' % dens_quantile)
    create_path_if_not_exist(store_path)
    plt.xlabel(xlabel)
    plt.ylabel('bucket count')
    plt.savefig(store_path)
    plt.show()


def heatmap(arr, plot_title):
    plt.title(plot_title)
    plt.imshow(arr, cmap='hot', interpolation='nearest')
    plt.show()


'''
cluster - both i,j in chromothripsis
cluster periphery - i or j in chromothripsis
sv - neither i or j in chromothripsis (simple sv)
'''


def heatmap_with_breakpoints_and_cluster(arr, plot_title, breakpoint_edges, cluster, save_path):
    plt.title(plot_title)
    plt.imshow(arr, cmap='hot', interpolation='nearest')
    labeled_sv = False
    for b in breakpoint_edges:
        if not labeled_sv:
            plt.scatter(b[0], b[1], s=150, c='blue', marker='o', label='sv')
            labeled_sv = True
        else:
            plt.scatter(b[0], b[1], s=150, c='blue', marker='o')
    labeled_periphery = False
    for b in breakpoint_edges:
        if b[0] in cluster or b[1] in cluster:
            if not labeled_periphery:
                plt.scatter(b[0], b[1], s=150, c='orange', marker='o', label='periphery cluster')
                labeled_periphery = True
            else:
                plt.scatter(b[0], b[1], s=150, c='orange', marker='o')
    labeled_cluster = False
    for b in breakpoint_edges:
        if b[0] in cluster and b[1] in cluster:
            if not labeled_cluster:
                plt.scatter(b[0], b[1], s=150, c='green', marker='o', label='cluster')
                labeled_cluster = True
            else:
                plt.scatter(b[0], b[1], s=150, c='green', marker='o')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(save_path)
    plt.show()


def perf_measure(y_actual, y_hat):
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    for i in range(len(y_hat)):
        if y_actual[i] == y_hat[i] == 1:
            TP += 1
        if y_hat[i] == 1 and y_actual[i] != y_hat[i]:
            FP += 1
        if y_actual[i] == y_hat[i] == 0:
            TN += 1
        if y_hat[i] == 0 and y_actual[i] != y_hat[i]:
            FN += 1

    return TP, FP, TN, FN


def hic(filepath, chr):
    c = cooler.Cooler(filepath)
    mat = c.matrix(balance=False).fetch(chr)
    mat_nan = zeros_to_nan(mat)
    mat_norm = normalize_intra(mat_nan)
    heatmap(mat_norm, f'normed_hic {chr}')


def bps_to_bins_with_resolution(bp1, bp2, resolution_bases):
    return int(bp1 / resolution_bases), int(bp2 / resolution_bases)


def create_path_if_not_exist(dirpath):
    if not os.path.exists(os.path.dirname(dirpath)):
        try:
            os.makedirs(os.path.dirname(dirpath))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def f_proximity(dist):
    return (1.0 / dist) ** 2


def analyze_donor(coolpath, csvpath, chr_num):
    c = cooler.Cooler(coolpath)
    resolution = c.info['bin-size']
    mat = c.matrix(balance=False).fetch(f'chr{chr_num}')
    mat_nan = zeros_to_nan(mat)
    mat_norm = normalize_intra(mat_nan)
    projection = ['chrom1', 'end1', 'end2']
    patient_csv = pd.read_csv(csvpath, usecols=projection)
    patient_csv_chr = patient_csv[patient_csv['chrom1'] == chr_num]
    bin_ij_coordinates = []
    for idx, row in patient_csv_chr.iterrows():
        bin_ij_coordinates.append(bps_to_bins_with_resolution(int(row['end1']), int(row['end2']), resolution))

    dist_mat = np.load(f'dists_npy/dists_npychr{chr_num}.npy')

    filepath = f'breakpoints/chr{chr_num}.txt'
    number_of_edges = 0
    number_of_nodes = -1
    all_bins = set()
    for b in bin_ij_coordinates:
        all_bins.add(b[0])
        all_bins.add(b[1])
    all_bins = list(all_bins)
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

    some_delta_just_for_sure = 5
    cluster_bins = WFind_Densest_Subgraph(number_of_nodes + some_delta_just_for_sure, number_of_edges, filepath)
    heatmap_with_breakpoints_and_cluster(mat_norm,
                                         f'normed hic & breakpoints chr{chr_num}\n{inspect.getsource(f_proximity)}',
                                         bin_ij_coordinates,
                                         cluster_bins,
                                         'breakpoints')
    print(WFind_Density(cluster_bins, filepath))
    print(f'clusters {cluster_bins}')
    periphery = set()
    for b in bin_ij_coordinates:
        i = b[0]
        j = b[1]
        if i in cluster_bins and j not in cluster_bins:
            periphery.add(j)
        if j in cluster_bins and i not in cluster_bins:
            periphery.add(i)
    print(f'periphery: {periphery}')
    print(f'all sv: {all_bins}')


def main():
    analyze_donor(coolpath='healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool',
                  csvpath='CGP_donor_1347756.csv',
                  chr_num=10)


if __name__ == '__main__':
    main()
