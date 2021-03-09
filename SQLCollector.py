import configparser
import inspect
import os
import re
import sys
import time
import warnings

import cooler
import psycopg2

from DonorAnalyzer import analyze_donor, f_proximity

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
sqlpath = script_dir + '/sql_scripts'


def unique_donors(cur):
    donors = []
    cur.execute(open(sqlpath + '/unique_donors_list.sql', "r").read())
    rows = cur.fetchall()
    for row in rows:
        donors.append(row[0])
    return donors


def unique_prostate_donors(cur):
    donors = []
    cur.execute('SELECT DISTINCT sv_intra.donor_id FROM sv_intra ' +
                'INNER JOIN donor_tumour ON ' +
                'sv_intra.donor_id = submitted_donor_id ' +
                'WHERE histology_tier2 = \'Prostate\'')
    rows = cur.fetchall()
    for row in rows:
        donors.append(row[0])
    return donors


def ddl(con):
    cur = con.cursor()
    try:
        cur.execute(open(sqlpath + '/ddl.sql', "r").read())
    except BaseException as e:
        warnings.warn(e.__str__())
        con.rollback()


def insert_proximity(con):
    cur = con.cursor()
    code = inspect.getsource(f_proximity)
    cur.execute(
        f'INSERT INTO proximity (function_code) VALUES (\'{code}\')'
    )
    con.commit()


def find_type_by_donor(cur, donor):
    cur.execute(f'SELECT histology_tier2 FROM donor_tumour WHERE submitted_donor_id = \'{donor}\'')
    rows = cur.fetchall()
    return rows


def main():
    t1 = time.time()

    config = configparser.ConfigParser()
    config.read(script_dir + r'/application.properties')
    user = config['DEFAULT']['pcawg.user']
    password = config['DEFAULT']['pcawg.password']

    c = cooler.Cooler('healthy_hics/Rao2014-IMR90-MboI-allreps-filtered.500kb.cool')

    with psycopg2.connect(
            database="pcawg",
            user=user,
            password=password,
            host="localhost",
            port="5432"
    ) as con:
        ddl(con)
        # insert_proximity(con) todo: only for first run!

        cur = con.cursor()
        donors = unique_donors(cur)
        donors = donors[:10]
        for donor in donors:
            analyze_donor(donor, con, c)

    t2 = time.time()
    print(f'took {t2 - t1} sec')


def check_donors():
    with psycopg2.connect(
            database="pcawg",
            user='',
            password='',
            host="localhost",
            port="5432"
    ) as con:
        cur = con.cursor()
        all = unique_donors(cur)
        all = set(all)
        print(len(all))
        # print(all, '\n')
        filepath = script_dir + '/donors.txt'
        with open(filepath, 'r') as filesmall:
            data = filesmall.readlines()
            data = list(map(lambda x: re.sub(r"[\n\t\s]*", "", x), data))
            print(len(data))
            # print(data)
            # for donor in data:
            #     if donor not in all:
            #         print(f'{donor} in all')


def donors_tumour_type_test():
    with psycopg2.connect(
            database="pcawg",
            user='',
            password='',
            host="localhost",
            port="5432"
    ) as con:
        cur = con.cursor()
        filepath = script_dir + '/donors.txt'
        with open(filepath, 'r') as filesmall:
            data = filesmall.readlines()
            data = list(map(lambda x: re.sub(r"[\n\t\s]*", "", x), data))
            for donor in data:
                t_type = find_type_by_donor(cur, donor)
                if not t_type:
                    warnings.warn('There is no such donor in the donor_tumour')
                else:
                    print(t_type)


def prostate():
    with psycopg2.connect(
            database="pcawg",
            user='',
            password='',
            host="localhost",
            port="5432"
    ) as con:
        cur = con.cursor()
        donors = unique_prostate_donors(cur)
        filepath = script_dir + '/my_donors.txt'
        with open(filepath, 'w') as file_donors:
            for donor in donors:
                file_donors.write("%s\n" % donor)
        # print(donors)
        print(len(donors))


if __name__ == '__main__':
    # main()
    # check_donors()
    # donors_tumour_type_test()
    prostate()
