import errno
import inspect

import numpy as np
import os
import sys
from DenseUtils import bps_to_bins_with_resolution, heatmap_with_breakpoints_and_cluster
from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density
from HiCUtils import zeros_to_nan, normalize_intra


def f_proximity(dist):
    return (1.0 / dist) ** 2

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
donorspath = script_dir + '/donors'


def analyze_donor(donor, cur, cooler):
    print(donor)
    regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
    cur.execute(
        f'SELECT * FROM sv_intra \
        WHERE donor_id = \'{donor}\' \
        AND chr SIMILAR TO \'{regex_up_to_21}\'')
    rows = cur.fetchall()

    chr_bins_map = {}

    resolution = cooler.info['bin-size']

    for row in rows:
        chr_n, bp1, bp2, _ = row
        bin1, bin2 = bps_to_bins_with_resolution(bp1, bp2, resolution)
        chr_breakpoints = chr_bins_map.get(chr_n, [])
        if not chr_breakpoints:
            chr_bins_map[chr_n] = [(bin1, bin2)]
        else:
            chr_bins_map[chr_n].append((bin1, bin2))

    for (chr_n, bin_pairs) in chr_bins_map.items():
        all_bins = set()
        number_of_edges = 0
        number_of_nodes = -1
        dist_mat = np.load(f'dists/dists_npy_chr{chr_n}.npy')
        for (x, y) in bin_pairs:
            all_bins.add(x)
            all_bins.add(y)
        all_bins = list(all_bins)

        dirpath = donorspath + f'/{donor}/{chr_n}'
        if not os.path.exists(os.path.dirname(dirpath)):
            try:
                os.makedirs(os.path.dirname(dirpath))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        filepath = donorspath + f'/{donor}/{chr_n}'

        mat = cooler.matrix(balance=False).fetch(f'chr{chr_n}')
        mat_nan = zeros_to_nan(mat)
        mat_norm = normalize_intra(mat_nan)

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
        heatmap_with_breakpoints_and_cluster(mat_norm,
                                             f'normed hic & breakpoints chr{chr_n}\n{inspect.getsource(f_proximity)}',
                                             bin_pairs,
                                             cluster_bins)
        print(WFind_Density(cluster_bins, filepath))
        print(f'clusters {cluster_bins}')
        periphery = set()
        for b in bin_pairs:
            i = b[0]
            j = b[1]
            if i in cluster_bins and j not in cluster_bins:
                periphery.add(j)
            if j in cluster_bins and i not in cluster_bins:
                periphery.add(i)
        print(f'periphery: {periphery}')
        print(f'all sv: {all_bins}')
