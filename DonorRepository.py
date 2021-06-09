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

    def random_frankenstein_ddl(self):
        try:
            self.cur.execute('SELECT ' +
                             'to_regclass(\'random_frankenstein\') ')
            rows = self.cur.fetchall()
            if None in rows[0]:
                self.cur.execute(open(self.sqlpath + '/random_frankenstein.sql', "r").read())
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def frankenstein_ddl(self):
        try:
            self.cur.execute('SELECT ' +
                             'to_regclass(\'frankenstein\') ')
            rows = self.cur.fetchall()
            if None in rows[0]:
                self.cur.execute(open(self.sqlpath + '/frankenstein_ddl.sql', "r").read())
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def ddl(self):
        try:
            self.cur.execute('SELECT ' +
                             'to_regclass(\'clusters\'), '
                             'to_regclass(\'proximity\'), '
                             'to_regclass(\'densities\'), '
                             'to_regclass(\'periphery\'), '
                             'to_regclass(\'proximity\') ')
            rows = self.cur.fetchall()
            if None in rows[0]:
                self.cur.execute(open(self.sqlpath + '/ddl.sql', "r").read())
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def get_random_frankenstein_donors(self):
        donors = []
        self.cur.execute('SELECT DISTINCT(donor_id) FROM random_frankenstein')
        rows = self.cur.fetchall()
        for row in rows:
            donors.append(row[0])
        return donors

    def unique_prostate_donors(self, _table) -> list:
        donors = []
        self.cur.execute(f'SELECT DISTINCT {_table}.donor_id FROM {_table} ' +
                         f'INNER JOIN donor_tumour ON {_table}.donor_id = submitted_donor_id ' +
                         'WHERE histology_tier2 = \'Prostate\'')
        rows = self.cur.fetchall()
        for row in rows:
            donors.append(row[0])
        return donors

    def get_pcawg(self, _table):
        donors = self.unique_prostate_donors(_table)
        donors = set(donors)
        donor_chr_pairs = set()
        for donor in donors:
            self.cur.execute(f'SELECT chr FROM {_table} WHERE donor_id = \'{donor}\'')
            rows = self.cur.fetchall()
            for row in rows:
                donor_chr_pairs.add((donor, row[0]))
        return donors, donor_chr_pairs

    def get_info_id(self, donor, chr_n, f_id):
        self.cur.execute(
            ' SELECT info_id FROM donorinfo ' +
            f'WHERE donor_id = \'{donor}\' AND chr = \'{chr_n}\' AND function_id = \'{f_id}\' '
        )
        row = self.cur.fetchone()
        if not row:
            return -1
        return row[0]

    def get_cluster(self, info_id) -> list:
        cluster = []
        self.cur.execute(
            ' SELECT bp FROM clusters ' +
            f'WHERE info_id = \'{info_id}\' '
        )
        rows = self.cur.fetchall()
        for row in rows:
            cluster.append(row[0])
        return cluster

    def get_density(self, info_id) -> float:
        self.cur.execute(
            ' SELECT density FROM densities ' +
            f'WHERE info_id = {info_id} '
        )
        row = self.cur.fetchone()
        return row[0]

    def get_periphery(self, info_id) -> list:
        periphery = []
        self.cur.execute(
            ' SELECT bp FROM periphery ' +
            f'WHERE info_id = \'{info_id}\' '
        )
        rows = self.cur.fetchall()
        for row in rows:
            periphery.append(row[0])
        return periphery

    def get_proximity_code(self, function_id):
        self.cur.execute(
            f'SELECT function_code FROM proximity WHERE function_id = {function_id}'
        )
        row = self.cur.fetchone()
        return row[0]

    def get_proximity_id(self, function_code):
        self.cur.execute(
            f'SELECT function_id FROM proximity WHERE function_code = \'{function_code}\''
        )
        row = self.cur.fetchone()
        return row[0]

    def get_by_infoid(self, info_id):
        self.cur.execute(
            ' SELECT donor_id, chr, function_id FROM donorinfo ' +
            f'WHERE info_id = {info_id} '
        )
        triplet = self.cur.fetchone()
        return triplet

    def get_donor_sv_chr_1_21(self, donor, _table):
        regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
        self.cur.execute(
            f'SELECT * FROM {_table} \
            WHERE donor_id = \'{donor}\' \
            AND chr SIMILAR TO \'{regex_up_to_21}\'')
        rows = self.cur.fetchall()
        return rows

    def get_donor_chrs_1_21(self, donor, _table):
        regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
        self.cur.execute(
            f'SELECT DISTINCT chr FROM {_table} \
            WHERE donor_id = \'{donor}\' \
            AND chr SIMILAR TO \'{regex_up_to_21}\'')
        rows = self.cur.fetchall()
        return rows

    def get_donor_chr_svs(self, donor, chr, _table):
        self.cur.execute(
            f'SELECT bp1, bp2 FROM {_table} \
                WHERE donor_id = \'{donor}\' \
                  AND chr = \'{chr}\'')
        rows = self.cur.fetchall()
        return rows

    def get_donor_fid_infoids(self, donor, fid):
        self.cur.execute(
            f'SELECT info_id FROM donorinfo \
                WHERE donor_id = \'{donor}\' \
                  AND function_id = \'{fid}\' '
        )
        rows = self.cur.fetchall()
        rows = list(map(lambda p: p[0], rows))
        return rows

    def get_seek(self):
        regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
        self.cur.execute(
            ' SELECT ' +
            '   donor_unique_id, ' +
            '   "Chr", ' +
            '   chromo_label ' +
            '  FROM chromo ' +
            f'  WHERE "Chr" SIMILAR TO \'{regex_up_to_21}\''
        )
        rows = self.cur.fetchall()
        return rows

    def get_seek_markup_by_donor_chr(self, donor, chr):
        self.cur.execute(
            ' SELECT ' +
            '   donor_unique_id, ' +
            '   "Chr", ' +
            '   "Start", ' +
            '   "End", ' +
            '   chromo_label ' +
            '  FROM chromo ' +
            f'  WHERE "Chr" SIMILAR TO \'{chr}\' ' +
            f'    AND donor_unique_id SIMILAR TO \'%{donor}%\' '
        )
        row = self.cur.fetchone()
        return row

    # -1 not found
    def get_chromo(self, donor, chr):
        self.cur.execute(
            ' SELECT ' +
            '   donor_unique_id, ' +
            '   "Chr", ' +
            '   "Start", ' +
            '   "End", ' +
            '   chromo_label ' +
            '  FROM chromo ' +
            f'WHERE donor_unique_id SIMILAR TO \'%{donor}\' ' +
            f'  AND "Chr" = \'{chr}\' '
        )
        row = self.cur.fetchone()
        if not row:
            return -1, -1, -1, -1, -1
        return row

    def insert_cluster(self, info_id, cluster: tuple):
        try:
            if len(cluster) == 1:  # tuple of size 1 is (value,) -> SQL IN throws error to pass arg
                self.cur.execute(
                    'SELECT * FROM clusters ' +
                    f'WHERE info_id = \'{info_id}\' ' +
                    f'  AND bp = {cluster[0]} '
                )
                rows = self.cur.fetchall()
            elif cluster:
                self.cur.execute(
                    'SELECT * FROM clusters ' +
                    f'WHERE info_id = \'{info_id}\' ' +
                    f'  AND bp IN {cluster} '
                )
                rows = self.cur.fetchall()
            else:
                rows = ()
            if len(rows) != len(cluster):
                for bp in cluster:
                    self.cur.execute(
                        f'INSERT INTO clusters (info_id, bp) VALUES (\'{info_id}\', \'{bp}\')'
                    )
                self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def insert_periphery(self, info_id, periphery: tuple):
        try:
            if len(periphery) == 1:  # tuple of size 1 is (value,) -> SQL IN throws error to pass arg
                self.cur.execute(
                    'SELECT * FROM periphery ' +
                    f'WHERE info_id = \'{info_id}\' ' +
                    f'  AND bp = {periphery[0]} '
                )
                rows = self.cur.fetchall()
            elif periphery:
                self.cur.execute(
                    'SELECT * FROM periphery ' +
                    f'WHERE info_id = \'{info_id}\' ' +
                    f'  AND bp IN {periphery} '
                )
                rows = self.cur.fetchall()
            else:
                rows = ()
            if len(rows) != len(periphery):
                for bp in periphery:
                    self.cur.execute(
                        f'INSERT INTO periphery (info_id, bp) VALUES (\'{info_id}\', \'{bp}\')'
                    )
                self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def insert_dense(self, info_id, dense):
        try:
            self.cur.execute(
                'SELECT * FROM densities ' +
                f'WHERE info_id = \'{info_id}\' ' +
                f'  AND density = \'{dense}\' '
            )
            rows = self.cur.fetchall()
            if not rows:
                self.cur.execute(
                    f'INSERT INTO densities (info_id, density) VALUES (\'{info_id}\', \'{dense}\')'
                )
                self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def insert_sv_frankenstein(self, donor, chr_n, bp1, bp2, _table):
        try:
            self.cur.execute(
                f'SELECT * FROM {_table} ' +
                f'WHERE donor_id = \'{donor}\' ' +
                f'  AND chr = \'{chr_n}\' ' +
                f'  AND bp1 = \'{bp1}\' ' +
                f'  AND bp2 = \'{bp2}\' '
            )
            rows = self.cur.fetchall()
            if not rows:
                self.cur.execute(
                    f'INSERT INTO {_table} (donor_id, chr, bp1, bp2) '
                    f'VALUES (\'{donor}\', \'{chr_n}\', \'{bp1}\', \'{bp2}\')'
                )
                self.con.commit()
        except BaseException as e:
            warnings.warn(e.__str__())
            self.con.rollback()

    def insert_donorinfo(self, donor, chr_n, f_id):
        try:
            self.cur.execute(
                'SELECT * FROM donorinfo ' +
                f'WHERE donor_id = \'{donor}\' ' +
                f'  AND chr = \'{chr_n}\' ' +
                f'  AND function_id = \'{f_id}\' '
            )
            rows = self.cur.fetchall()
            if not rows:
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
                ' SELECT * FROM proximity ' +
                f'WHERE function_code = \'{code}\' '
            )
            rows = self.cur.fetchall()
            if not rows:
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
