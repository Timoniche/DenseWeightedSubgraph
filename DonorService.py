def bps_to_bins_with_resolution(bp1, bp2, resolution_bases):
    return int(bp1 / resolution_bases), int(bp2 / resolution_bases)


def collect_chr_bins_map_with_resolution(svs, resolution):
    chr_bins_map = {}
    for sv in svs:
        chr_n, bp1, bp2, _ = sv
        bin1, bin2 = bps_to_bins_with_resolution(bp1, bp2, resolution)
        chr_breakpoints = chr_bins_map.get(chr_n, [])
        if not chr_breakpoints:
            chr_bins_map[chr_n] = [(bin1, bin2)]
        else:
            chr_bins_map[chr_n].append((bin1, bin2))
    return chr_bins_map
