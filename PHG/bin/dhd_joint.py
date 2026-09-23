#!/usr/bin/env python3
"""Count-once donor allele table: each donor's allele at every biallelic union site, from ONE joint step 4 on per-sample counts.
Per donor: tier, n (reads), n_pools_alt, and two states (1 = ALT, 0 = REF, NA = unknown):
  state_A  : ALT = tier A                                   REF = tier ref
  state_AB : ALT = tier A, or tier B with ALT in >1 BC1 sample   REF = tier ref   (the A+B founder rule, restricted to union sites)
Output: OUTDIR/dhd_<chrom>.tsv.gz, OUTDIR/dhd_summary.tsv
Usage: dhd_joint.py OUTDIR union.tsv.gz DONOR=step4/DONOR.sites.tsv.gz [DONOR=...]"""
import sys, gzip, os, collections

out, uf = sys.argv[1], sys.argv[2]; tabs = dict(a.split('=', 1) for a in sys.argv[3:])
def load(p):
    t = {}
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t')
            t[(x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']])] = (x[c['tier']], x[c['n']], int(x[c['n_pools_alt']]))
    return t
T = {d: load(p) for d, p in tabs.items()}; donors = list(tabs)
union = []
with gzip.open(uf, 'rt') as f:
    f.readline()
    for l in f:
        x = l.rstrip('\n').split('\t')
        if x[6] == '0': union.append((x[0], int(x[1]), x[2], x[3], x[5]))
S = {d: collections.Counter() for d in donors}
chrom = union[0][0] if union else 'chr'
with gzip.open(os.path.join(out, f'dhd_{chrom}.tsv.gz'), 'wt') as o:
    o.write('chrom\tpos\tref\talt\tdiscovered_in\t' + '\t'.join(f'{d}_tier\t{d}_n\t{d}_npa\t{d}_state_A\t{d}_state_AB' for d in donors) + '\n')
    for c, p, r, a, disc in union:
        cells = []
        for d in donors:
            tier, n, npa = T[d].get((c, p, r, a), ('absent', '0', 0))
            sA = '1' if tier == 'A' else '0' if tier == 'ref' else 'NA'
            sAB = '1' if tier == 'A' or (tier == 'B' and npa >= 2) else '0' if tier == 'ref' else 'NA'
            own = 'own' if d in disc.split(',') else 'other'
            S[d][f'tier_{tier}_{own}'] += 1; S[d][f'A_{sA}'] += 1; S[d][f'AB_{sAB}'] += 1
            cells += [tier, n, str(npa), sA, sAB]
        o.write(f"{c}\t{p}\t{r}\t{a}\t{disc}\t" + '\t'.join(cells) + '\n')
n = len(union)
with open(os.path.join(out, 'dhd_summary.tsv'), 'w') as o:
    ks = sorted({k for d in donors for k in S[d]})
    o.write('donor\t' + '\t'.join(ks) + '\n')
    for d in donors: o.write(d + '\t' + '\t'.join(str(S[d][k]) for k in ks) + '\n')
print(f"[dhd_joint] biallelic union sites: {n}")
for d in donors:
    s = S[d]; pct = lambda v: f"{v} ({100*v/n:.1f}%)"
    print(f"[dhd_joint] {d} rule A : ALT {pct(s['A_1'])} | REF {pct(s['A_0'])} | unknown {pct(s['A_NA'])}")
    print(f"[dhd_joint] {d} rule AB: ALT {pct(s['AB_1'])} | REF {pct(s['AB_0'])} | unknown {pct(s['AB_NA'])}")
    own = {t: s[f'tier_{t}_own'] for t in ('A', 'B', 'C', 'ref', '-', 'absent')}
    print(f"[dhd_joint] {d} own tier-A discovery sites re-called: " + ', '.join(f"{t} {v}" for t, v in own.items()))
