#!/usr/bin/env python3
"""Step 3 of docs/PLAN_marker_union_pilot.md: the donor's allele (DHd) at every biallelic union site.
Per donor and union allele (chrom, pos, ref, alt): tier from the donor's own step-4 table (discovery) if the allele is there, else from
the second-pass step-4 table (no veto). State (pilot rule): ALT = tier A, REF = tier 'ref', missing = anything else (B, C, '-', absent).
Output: OUTDIR/dhd_union_<chrom>.tsv.gz  chrom pos ref alt <donor>_state <donor>_tier <donor>_src ... (state 1 = ALT, 0 = REF, NA)
        OUTDIR/dhd_summary.tsv  donor x (ALT own / ALT second pass / REF own / REF second pass / missing by tier)
Usage: dhd_union.py OUTDIR union.tsv.gz DONOR=own.sites.tsv.gz,pass2.sites.tsv.gz [DONOR=...]"""
import sys, gzip, os, collections

out, uf = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
spec = dict(a.split('=', 1) for a in sys.argv[3:])

def tiers(path):
    t = {}
    with gzip.open(path, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t'); t[(x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']])] = x[c['tier']]
    return t

tabs = {}
for d, v in spec.items():
    own, p2 = v.split(','); tabs[d] = (tiers(own), tiers(p2))
union = []
with gzip.open(uf, 'rt') as f:
    f.readline()
    for l in f:
        x = l.rstrip('\n').split('\t')
        if x[6] == '1': continue                                    # multi-allelic positions dropped
        union.append((x[0], int(x[1]), x[2], x[3]))
donors = list(spec); summ = {d: collections.Counter() for d in donors}
chrom = union[0][0] if union else 'chr'
with gzip.open(os.path.join(out, f'dhd_union_{chrom}.tsv.gz'), 'wt') as o:
    o.write('chrom\tpos\tref\talt\t' + '\t'.join(f'{d}_state\t{d}_tier\t{d}_src' for d in donors) + '\n')
    for k in union:
        cells = []
        for d in donors:
            own, p2 = tabs[d]
            if k in own: t, src = own[k], 'own'
            elif k in p2: t, src = p2[k], 'pass2'
            else: t, src = 'absent', 'none'
            st = '1' if t == 'A' else '0' if t == 'ref' else 'NA'
            lab = 'ALT' if st == '1' else 'REF' if st == '0' else 'missing'
            summ[d][f'{lab}_{src}'] += 1; summ[d][lab] += 1
            if st == 'NA': summ[d][f'missing_tier_{t}'] += 1
            cells += [st, t, src]
        o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t" + '\t'.join(cells) + '\n')
keys = ['ALT', 'ALT_own', 'ALT_pass2', 'REF', 'REF_own', 'REF_pass2', 'missing', 'missing_tier_B', 'missing_tier_C', 'missing_tier_-', 'missing_tier_absent']
with open(os.path.join(out, 'dhd_summary.tsv'), 'w') as o:
    o.write('donor\tunion_sites\t' + '\t'.join(keys) + '\n')
    for d in donors: o.write(f"{d}\t{len(union)}\t" + '\t'.join(str(summ[d][k]) for k in keys) + '\n')
print(f"[dhd_union] biallelic union sites: {len(union)}")
for d in donors:
    s = summ[d]; n = len(union)
    print(f"[dhd_union] {d}: ALT {s['ALT']} ({100*s['ALT']/n:.1f}%; own {s['ALT_own']}, pass2 {s['ALT_pass2']}) | REF {s['REF']} ({100*s['REF']/n:.1f}%; "
          f"own {s['REF_own']}, pass2 {s['REF_pass2']}) | missing {s['missing']} ({100*s['missing']/n:.1f}%; B {s['missing_tier_B']}, C {s['missing_tier_C']}, "
          f"- {s['missing_tier_-']}, absent {s['missing_tier_absent']})")
