def bps_to_bins_with_resolution(bp1, bp2, resolution_bases):
    return int(bp1 / resolution_bases), int(bp2 / resolution_bases)


# collect svs in bp -> map them to bins -> filter svs with predicate NOT |bin1 - bin2| <= 1
def collect_chr_bins_map_with_resolution(svs, resolution):
    chr_bins_map = {}
    for sv in svs:
        chr_n, bp1, bp2, _ = sv
        bin1, bin2 = bps_to_bins_with_resolution(bp1, bp2, resolution)
        if abs(bin1 - bin2) <= 1:
            continue
        chr_breakpoints = chr_bins_map.get(chr_n, [])
        if not chr_breakpoints:
            chr_bins_map[chr_n] = [(bin1, bin2)]
        else:
            chr_bins_map[chr_n].append((bin1, bin2))
    return chr_bins_map


# collect svs in bp -> filter svs with predicate NOT |bp1 - bp2| <= resolution
def filter_svs_with_resolution(svs, resolution):
    chr_sv_map = {}
    for sv in svs:
        chr_n, bp1, bp2, _ = sv
        # if abs(bp1 - bp2) <= resolution:
        bin1, bin2 = bps_to_bins_with_resolution(bp1, bp2, resolution)
        if abs(bin1 - bin2) <= 1:
            continue
        chr_breakpoints = chr_sv_map.get(chr_n, [])
        if not chr_breakpoints:
            chr_sv_map[chr_n] = [(bp1, bp2)]
        else:
            chr_sv_map[chr_n].append((bp1, bp2))
    return chr_sv_map
