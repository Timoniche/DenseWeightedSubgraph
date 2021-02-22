import pandas as pd
from IPython.core.display import display


def main():
    f = pd.read_csv('combined_csv.csv')
    intra_csv = f[f['chrom1'] == f['chrom2']]
    intra_csv.to_csv('intra_csv.csv', index=False)


if __name__ == '__main__':
    main()