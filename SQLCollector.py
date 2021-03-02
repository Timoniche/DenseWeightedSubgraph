import configparser
import os
import sys
import time
import warnings

import cooler
import psycopg2

from DonorAnalyzer import analyze_donor

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
sqlpath = script_dir + '/sql_scripts'


def unique_donors(cur):
    donors = []
    cur.execute(open(sqlpath + '/unique_donors_list.sql', "r").read())
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
        cur = con.cursor()
        donors = unique_donors(cur)
        donors = donors[:10]
        for donor in donors:
            analyze_donor(donor, cur, c)

    t2 = time.time()
    print(f'took {t2 - t1} sec')


if __name__ == '__main__':
    main()
