#!/usr/bin/env python3
"""QC-set benchmark of the marker-union gap filling, scoring step: every gap call of dhd_bayes.py against the assembly truth.
Truth per founder = its PERFECT founder alleles (assembly SNPs vs B73 in the lowcopy ranges; build_founder_gvcf.py .alt.tsv):
  truth ALT = (chrom, pos, ref, alt) in the founder's truth set;  truth REF = not in it (the assembly matches B73 there, within the ranges).
Reports per founder: own sites (discovered, subsampled) precision; gap calls ALT / REF / missing vs truth, for the Bayes rule and the
heuristic A+B rule; and gap ALT precision by posterior bin (the cut-off question).
Usage: union_bench_score.py dhd_bayes_chr10.tsv.gz OUT.tsv FOUNDER=truth.alt.tsv [FOUNDER=...]"""
import sys, gzip, collections
fb, fo = sys.argv[1], sys.argv[2]; truth = {}
for spec in sys.argv[3:]:
    f, p = spec.split('=', 1); s = set()
    for l in open(p):
        x = l.split('\t'); s.add((x[0], x[1], x[2], x[3].strip()))
    truth[f] = s; print(f"[union_bench_score] {f}: truth alleles {len(s)}")
bins = [(0, 0.5), (0.5, 0.9), (0.9, 0.95), (0.95, 0.99), (0.99, 0.999), (0.999, 1.01)]
R = collections.defaultdict(collections.Counter)
with gzip.open(fb, 'rt') as f:
    h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
    for l in f:
        x = l.rstrip('\n').split('\t'); k = tuple(x[:4])
        for F in truth:
            t = 'ALT' if k in truth[F] else 'REF'; src = x[c[f'{F}_src']]
            if src == 'own': R[F][('own', t)] += 1; continue
            sb, sab, tier, post = x[c[f'{F}_state_bayes']], x[c[f'{F}_state_AB']], x[c[f'{F}_tier']], x[c[f'{F}_posterior']]
            lab = {'1': 'ALT', '0': 'REF', 'NA': 'missing'}
            R[F][('gap_bayes_' + lab[sb], t)] += 1; R[F][('gap_AB_' + lab[sab], t)] += 1; R[F][('gap_all', t)] += 1
            if tier in ('A', 'B', 'C', '-') and post != 'NA':
                p = float(post)
                for lo, hi in bins:
                    if lo <= p < hi: R[F][(f'gap_tier{"A" if tier == "A" else "BC-"}_post_{lo}-{hi if hi <= 1 else 1}', t)] += 1; break
with open(fo, 'w') as o:
    o.write('founder\tclass\ttruth_ALT\ttruth_REF\tn\tfrac_truth_ALT\n')
    for F in truth:
        classes = sorted({k[0] for k in R[F]})
        for cl in classes:
            a, r = R[F][(cl, 'ALT')], R[F][(cl, 'REF')]; n = a + r
            line = f"{F}\t{cl}\t{a}\t{r}\t{n}\t{a / n:.3f}" if n else f"{F}\t{cl}\t0\t0\t0\tNA"
            o.write(line + '\n'); print('[union_bench_score] ' + line)
