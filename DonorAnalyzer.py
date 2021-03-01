def f_proximity(dist):
    return (1.0 / dist) ** 2


def analyze_donor(donor, cur):
    cur.execute(f'SELECT * FROM sv_intra \
                WHERE donor_id = \'{donor}\'')
    rows = cur.fetchall()
    for row in rows:
        print(row)
