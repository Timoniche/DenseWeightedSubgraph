import glob
import pandas as pd
import os
import sys

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')


def main():
    os.chdir('icgc/open')
    csvs = []
    for file in glob.glob('*.gz'):
        csvs.append(pd.read_csv(file, sep='\t'))
    combined_csv = pd.concat(csvs)
    combined_csv.to_csv(script_dir + '/combined_csv.csv', index=False)


if __name__ == '__main__':
    main()
