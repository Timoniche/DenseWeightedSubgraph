from DonorRepository import DonorRepository
import random


def collect_frankenshteins():
    CHR_SV_PROVIDER_POOL = {}
    with DonorRepository() as rep:
        rep.frankenstein_ddl()
        prostate_pcawg_donors, prostate_pcawg_pairs = rep.get_pcawg(_table='sv_intra')
        donor_chr_sv_sz_map = {}
        for donor in prostate_pcawg_donors:
            svs = rep.get_donor_sv_chr_1_21(donor, _table='sv_intra')
            for sv in svs:
                chrN, bp1, bp2, _ = sv
                if chrN != '22' and chrN != 'X' and chrN != 'Y':
                    chr_breakpoints = CHR_SV_PROVIDER_POOL.get(chrN, [])
                    if not chr_breakpoints:
                        CHR_SV_PROVIDER_POOL[chrN] = [(bp1, bp2)]
                    else:
                        CHR_SV_PROVIDER_POOL[chrN].append((bp1, bp2))

                    cur_sz = donor_chr_sv_sz_map.get((donor, chrN), 0)
                    donor_chr_sv_sz_map[(donor, chrN)] = cur_sz + 1
        for i in range(1, 22):
            random.shuffle(CHR_SV_PROVIDER_POOL[str(i)])

        for (_donor, _chr) in prostate_pcawg_pairs:
            if _chr != '22' and _chr != 'X' and _chr != 'Y':
                sz = donor_chr_sv_sz_map[(_donor, _chr)]
                for i in range(sz):
                    bp1, bp2 = CHR_SV_PROVIDER_POOL[_chr].pop()
                    rep.insert_sv_frankenstein(_donor, _chr, bp1, bp2, 'frankenstein')



def collect_random_endpoints():
    CHR_BP_POOL = {}
    CHR_SV_SIZES = {}
    with DonorRepository() as rep:
        rep.random_frankenstein_ddl()
        prostate_pcawg_donors, prostate_pcawg_pairs = rep.get_pcawg(_table='sv_intra')
        for donor in prostate_pcawg_donors:
            svs = rep.get_donor_sv_chr_1_21(donor, _table='sv_intra')
            donor_chr_sv_cnt = {}
            for sv in svs:
                chrN, bp1, bp2, _ = sv
                if chrN != '22' and chrN != 'X' and chrN != 'Y':
                    chr_breakpoints = CHR_BP_POOL.get(chrN, [])
                    if not chr_breakpoints:
                        CHR_BP_POOL[chrN] = [bp1, bp2]
                    else:
                        CHR_BP_POOL[chrN].append(bp1)
                        CHR_BP_POOL[chrN].append(bp2)

                    donor_chr_sv_cnt[chrN] = donor_chr_sv_cnt.get(chrN, 0) + 1

            for (chrN, cnt) in donor_chr_sv_cnt.items():
                cur_sz_arr = CHR_SV_SIZES.get(chrN, [])
                if not cur_sz_arr:
                    CHR_SV_SIZES[chrN] = [cnt]
                else:
                    CHR_SV_SIZES[chrN].append(cnt)

        for i in range(1, 22):
            random.shuffle(CHR_BP_POOL[str(i)])

        DONOR_CHR_COUNT = 2424
        DONORS_CNT = 187
        donor_chr_pairs = set()
        for i in range(DONOR_CHR_COUNT):
            _chr, _random_svs = generate_frankenstein_donor(CHR_BP_POOL, CHR_SV_SIZES)
            donor_i = random.choice([i for i in range(1, DONORS_CNT + 1)])
            while (donor_i, _chr) in donor_chr_pairs:
                donor_i = random.choice([i for i in range(1, DONORS_CNT + 1)])
            donor_chr_pairs.add((donor_i, _chr))
            for bp1, bp2 in _random_svs:
                rep.insert_sv_frankenstein(f'random{donor_i}', _chr, bp1, bp2, 'random_frankenstein')


def normalize_list(raw):
    return [float(i) / sum(raw) for i in raw]


def generate_frankenstein_donor(CHR_BP_POOL: dict, CHR_SV_SIZES: dict):
    chrKeys = list(CHR_SV_SIZES.keys())
    chrProbs = normalize_list(list(map(lambda x: len(x), CHR_SV_SIZES.values())))
    chrN = random.choices(chrKeys, chrProbs)[0]
    print(chrN)
    random_sv_sz = random.choice(CHR_SV_SIZES[chrN])
    print(random_sv_sz)
    random_svs = []
    for i in range(random_sv_sz):
        bp1 = random.choice(CHR_BP_POOL[chrN])
        bp2 = random.choice(CHR_BP_POOL[chrN])
        random_svs.append((bp1, bp2))
    print(random_svs)
    return chrN, random_svs

def main():
    # collect_frankenshteins()
    collect_random_endpoints()


if __name__ == '__main__':
    main()
