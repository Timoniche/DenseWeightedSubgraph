import errno
import inspect
import time

import cooler
import numpy as np
import os
import sys
from DenseUtils import heatmap_with_breakpoints_and_cluster, create_path_if_not_exist
from DonorRepository import DonorRepository
from DonorService import collect_chr_bins_map_with_resolution
from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density
from HiCUtils import zeros_to_nan, normalize_intra


def f_proximity(dist):
    return (1.0 / dist) ** 2


script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
donorspath = script_dir + '/donors'


def analyze_donor(donor, cooler, f_id, rep: DonorRepository):
    print(donor)
    svs = rep.get_donor_sv_chr_1_21(donor)

    resolution = cooler.info['bin-size']
    chr_bins_map = collect_chr_bins_map_with_resolution(svs, resolution)

    for (chr_n, bin_pairs) in chr_bins_map.items():
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
                                             cluster_bins,
                                             save_path=donorspath + f'/{donor}/{chr_n}.png')
        dens = WFind_Density(cluster_bins, filepath)
        print(dens)
        print(f'clusters {cluster_bins}')
        periphery = set()
        info_id = rep.get_info_id(donor, chr_n, f_id)
        rep.insert_dense(info_id, dens)
        for b in bin_pairs:
            i = b[0]
            j = b[1]
            if i in cluster_bins and j not in cluster_bins:
                periphery.add(j)
            if j in cluster_bins and i not in cluster_bins:
                periphery.add(i)
        print(f'periphery: {periphery}')
        print(f'all sv: {all_bins}')


def main():
    t1 = time.time()

    cool = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')
    with DonorRepository() as rep:
        rep.ddl()
        rep.insert_proximity(f_proximity)

        donors = rep.unique_prostate_donors()
        f_count = 1
        for donor in donors[:5]:
            for f_id in range(1, f_count + 1):
                analyze_donor(donor=donor, cooler=cool, f_id=f_id, rep=rep)

    t2 = time.time()
    print(f'took {t2 - t1} sec')


if __name__ == '__main__':
    main()
