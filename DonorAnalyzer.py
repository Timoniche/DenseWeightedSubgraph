def f_proximity(dist):
    return (1.0 / dist) ** 2


def analyze_donor(donor, cur):
    regex_up_to_21 = '(2[0-1]|1[0-9]|[1-9])'
    cur.execute(
        f'SELECT * FROM sv_intra \
        WHERE donor_id = \'{donor}\' \
        AND chr SIMILAR TO \'{regex_up_to_21}\'')
    rows = cur.fetchall()

    chr_bins_map = {}

    for row in rows:
        chr, bp1, bp2, _ = row
        chr_breakpoints = chr_bins_map.get(chr, [])
        if not chr_breakpoints:
            chr_bins_map[chr] = [(bp1, bp2)]
        else:
            chr_bins_map[chr].append((bp1, bp2))

    print(donor)
    print(chr_bins_map)
    print()
