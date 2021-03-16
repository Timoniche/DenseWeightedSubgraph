import configparser
import inspect
import os
import sys
import warnings

import psycopg2


class DonorRepository:

    def __init__(self) -> None:
        script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or '.')
        self.sqlpath = script_dir + '/sql_scripts'
        config = configparser.ConfigParser()
        config.read(script_dir + r'/application.properties')
        user = config['DEFAULT']['pcawg.user']
        password = config['DEFAULT']['pcawg.password']
        self.con = psycopg2.connect(
            database="pcawg",
            user=user,
            password=password,
            host="localhost",
            port="5432"
        )
        self.cur = self.con.cursor()

    def ddl(self):
        try:
            self.cur.execute(open(self.sqlpath + '/ddl.sql', "r").read())
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def unique_prostate_donors(self) -> list:
        donors = []
        self.cur.execute('SELECT DISTINCT sv_intra.donor_id FROM sv_intra ' +
                         'INNER JOIN donor_tumour ON sv_intra.donor_id = submitted_donor_id ' +
                         'WHERE histology_tier2 = \'Prostate\'')
        rows = self.cur.fetchall()
        for row in rows:
            donors.append(row[0])
        return donors

    def get_info_id(self, donor, chr_n, f_id):
        self.cur.execute(
            ' SELECT info_id FROM donorinfo ' +
            f'WHERE donor_id = \'{donor}\' AND chr = \'{chr_n}\' AND function_id = \'{f_id}\' '
        )
        row = self.cur.fetchone()
        return row[0]

    def get_donor_sv_chr_1_21(self, donor):
        regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
        self.cur.execute(
            f'SELECT * FROM sv_intra \
            WHERE donor_id = \'{donor}\' \
            AND chr SIMILAR TO \'{regex_up_to_21}\'')
        rows = self.cur.fetchall()
        return rows

    def insert_dense(self, info_id, dense):
        try:
            self.cur.execute(
                f'INSERT INTO densities (info_id, density) VALUES (\'{info_id}\', \'{dense}\')'
            )
            self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def insert_donorinfo(self, donor, chr_n, f_id):
        try:
            self.cur.execute(
                f'INSERT INTO donorinfo (donor_id, chr, function_id) VALUES (\'{donor}\', \'{chr_n}\', \'{f_id}\')'
            )
            self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def insert_proximity(self, f):
        try:
            code = inspect.getsource(f)
            self.cur.execute(
                f'INSERT INTO proximity (function_code) VALUES (\'{code}\')'
            )
            self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.con.close()
