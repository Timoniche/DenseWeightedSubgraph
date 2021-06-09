from DonorRepository import DonorRepository
import random


def collect_frankenshteins():
    CHR_SV_PROVIDER_POOL = {}
    with DonorRepository() as rep:
        #rep.frankenstein_ddl()
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

        summm = 0
        for i in range(1, 22):
            summm += len(CHR_SV_PROVIDER_POOL[str(i)])

        for (_donor, _chr) in prostate_pcawg_pairs:
            if _chr != '22' and _chr != 'X' and _chr != 'Y':
                sz = donor_chr_sv_sz_map[(_donor, _chr)]
                for i in range(sz):
                    bp1, bp2 = CHR_SV_PROVIDER_POOL[_chr].pop()
                    rep.insert_sv_frankenstein(_donor, _chr, bp1, bp2)
        print('Hello!')


def main():
    collect_frankenshteins()


if __name__ == '__main__':
    main()
