import matplotlib.pyplot as plt
import cooler
import pandas as pd
from IPython.core.display import display

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


def heatmap(arr, plot_title):
    plt.title(plot_title)
    plt.imshow(arr, cmap='hot', interpolation='nearest')
    plt.show()


def hic(filepath, chr):
    c = cooler.Cooler(filepath)
    mat = c.matrix(balance=False).fetch(chr)
    mat_nan = zeros_to_nan(mat)
    mat_norm = normalize_intra(mat_nan)
    heatmap(mat_norm, f'normed_hic {chr}')


def bps_to_bins_with_resolution(bp1, bp2, resolution_bases):
    return int(bp1 / resolution_bases), int(bp2 / resolution_bases)


def analyze_donor(coolpath, csvpath, chr_num):
    c = cooler.Cooler(coolpath)
    resolution = c.info['bin-size']
    projection = ['chrom1', 'end1', 'end2']
    patient_csv = pd.read_csv(csvpath, usecols=projection)
    patient_csv_chr = patient_csv[patient_csv['chrom1'] == chr_num]
    # display(patient_csv_chr)
    for idx, row in patient_csv_chr.iterrows():
        print(bps_to_bins_with_resolution(row['end1'], row['end2'], resolution))


def main():
    analyze_donor(coolpath='healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool',
                  csvpath='CGP_donor_1347756.csv',
                  chr_num=10)

if __name__ == '__main__':
    main()
