import configparser
import os
import sys
import psycopg2

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
sqlpath = script_dir + '/sql_scripts'


def unique_donors(cur):
    donors = []
    cur.execute(open(sqlpath + '/unique_donors_list.sql', "r").read())
    rows = cur.fetchall()
    for row in rows:
        donors.append(row[0])
    return donors


def collect_donor_info(donor, cur):
    pass


def main():
    config = configparser.ConfigParser()
    config.read(script_dir + r'/application.properties')
    user = config['DEFAULT']['pcawg.user']
    password = config['DEFAULT']['pcawg.password']

    with psycopg2.connect(
            database="pcawg",
            user=user,
            password=password,
            host="localhost",
            port="5432"
    ) as con:
        cur = con.cursor()
        donors = unique_donors(cur)
        donors = donors[:10]
        for donor in donors:
            # collect_donor_info(donor, cur)
            print(donor)


if __name__ == '__main__':
    main()
