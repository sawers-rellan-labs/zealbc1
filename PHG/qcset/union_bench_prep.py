#!/usr/bin/env python3
"""QC-set benchmark of the marker-union gap filling (docs/PLAN_marker_union_pilot.md), prep step.
Random-subsamples each founder's discovered tier-A sites (QC-set step 4) to the size discovered in the real pilot, then writes the union
in the format of PHG/bin/union_sites.py (chrom pos ref alt n_donors donors multiallelic), so the same count-once + dhd_bayes chain runs.
Usage: union_bench_prep.py OUTDIR seed FOUNDER=step4/FOUNDER.sites.tsv.gz:N [FOUNDER=...:N]"""
import sys, gzip, os, random, collections
out, seed = sys.argv[1], int(sys.argv[2]); os.makedirs(out, exist_ok=True); random.seed(seed)
carriers = collections.defaultdict(list)
for spec in sys.argv[3:]:
    f, rest = spec.split('=', 1); path, n = rest.rsplit(':', 1); n = int(n); A = []
    with gzip.open(path, 'rt') as g:
        h = g.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in g:
            x = l.rstrip('\n').split('\t')
            if x[c['tier']] == 'A': A.append((x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']]))
    S = random.sample(A, min(n, len(A))); print(f"[union_bench_prep] {f}: tier A {len(A)} -> subsample {len(S)}")
    for k in S: carriers[k].append(f)
alts = collections.defaultdict(set)
for k in carriers: alts[k[:2]].add(k[3])
multi = {p for p, a in alts.items() if len(a) > 1}
with gzip.open(os.path.join(out, 'union_chr10.tsv.gz'), 'wt') as o:
    o.write('chrom\tpos\tref\talt\tn_donors\tdonors\tmultiallelic\n')
    for k in sorted(carriers, key=lambda z: (z[0], z[1], z[3])):
        o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{len(carriers[k])}\t{','.join(sorted(carriers[k]))}\t{int(k[:2] in multi)}\n")
nd = collections.Counter(len(v) for v in carriers.values())
print(f"[union_bench_prep] union alleles {len(carriers)} | positions {len(alts)} | multi-allelic positions {len(multi)} | by n donors {dict(nd)}")
