def f_proximity(dist):
    return (1.0 / dist) ** 2


def analyze_donor(donor, cur):
    regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
    cur.execute(
        f'SELECT * FROM sv_intra \
        WHERE donor_id = \'{donor}\' \
        AND chr SIMILAR TO \'{regex_up_to_21}\'')

    rows = cur.fetchall()
    for row in rows:
        print(row)
