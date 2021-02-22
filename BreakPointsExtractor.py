import glob
import pandas as pd
import os
import sys

from IPython.core.display import display

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')


def main():
    donors_file = pd.read_excel(script_dir + r'/pcawg-data-releases.xlsx')
    donors_file.to_csv(script_dir + r'/donors_csv.csv', index=False)
    projection = ['tumor_wgs_aliquot_id', 'submitter_donor_id']
    filehash_id = pd.read_csv(script_dir + r'/donors_csv.csv', usecols=projection)
    filehash_id.to_csv(script_dir + r'/filehash_id.csv', index=False)
    ids_by_hash = pd.read_csv(script_dir + r'/filehash_id.csv')

    filehash_donorid_map = {}
    for idx, row in ids_by_hash.iterrows():
        filehashes = row['tumor_wgs_aliquot_id'].split(sep=',')
        for hash in filehashes:
            filehash_donorid_map[hash] = row['submitter_donor_id']

    os.chdir('icgc/open')
    csvs = []
    for file in glob.glob('*.gz'):
        file_prefix = file.split(sep='.')[0]
        donor_id = filehash_donorid_map.get(file_prefix)
        cur_csv = pd.read_csv(file, sep='\t')
        cur_csv['donor_id'] = donor_id
        csvs.append(cur_csv)
    combined_csv = pd.concat(csvs)
    combined_csv.to_csv(script_dir + r'/combined_csv.csv', index=False)


if __name__ == '__main__':
    main()
