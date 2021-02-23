import pandas as pd

donor_id = 'CGP_donor_1347756'


def main():
    f = pd.read_csv('combined_csv.csv')
    intra_csv = f[f['chrom1'] == f['chrom2']]
    intra_csv_of_id = intra_csv[intra_csv['donor_id'] == donor_id]
    # intra_csv.to_csv('intra_csv.csv', index=False)
    intra_csv_of_id.to_csv('CGP_donor_1347756.csv', index=False)


if __name__ == '__main__':
    main()