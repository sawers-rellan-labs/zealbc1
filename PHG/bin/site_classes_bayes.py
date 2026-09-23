#!/usr/bin/env python3
"""Site classes for the per-site test of one donor's gap filling (scored by compare_discovery_sets.py score):
  own_ALT           the donor's own discovered sites (baseline for a real ALT; random sample)
  gap_REF           gap sites called REF from the donor's reads (tier ref; random sample)
  gap_tierA         gap sites with tier A in the donor's reads
  gap_post_<bin>    gap sites with tier B / C / '-' binned by Bayes posterior: <0.5, 0.5-0.9, 0.9-0.95, 0.95-0.99, 0.99-0.999, >=0.999
Usage: site_classes_bayes.py DONOR dhd_bayes.tsv.gz OUTDIR [n_base=5000]"""
import sys, gzip, os, random
d, fb, out = sys.argv[1:4]; nb = int(sys.argv[4]) if len(sys.argv) > 4 else 5000; os.makedirs(out, exist_ok=True)
bins = [(0.0, 0.5, 'lt0.5'), (0.5, 0.9, '0.5-0.9'), (0.9, 0.95, '0.9-0.95'), (0.95, 0.99, '0.95-0.99'), (0.99, 0.999, '0.99-0.999'), (0.999, 1.01, 'ge0.999')]
cl = {'own_ALT': [], 'gap_REF': [], 'gap_tierA': []}; cl.update({f'gap_post_{b[2]}': [] for b in bins})
with gzip.open(fb, 'rt') as f:
    h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
    for l in f:
        x = l.rstrip('\n').split('\t'); k = tuple(x[:4]); src, tier, post = x[c[f'{d}_src']], x[c[f'{d}_tier']], x[c[f'{d}_posterior']]
        if src == 'own': cl['own_ALT'].append(k); continue
        if tier == 'ref': cl['gap_REF'].append(k)
        elif tier == 'A': cl['gap_tierA'].append(k)
        elif tier in ('B', 'C', '-') and post != 'NA':
            p = float(post)
            for lo, hi, name in bins:
                if lo <= p < hi: cl[f'gap_post_{name}'].append(k); break
random.seed(1)
for b in ('own_ALT', 'gap_REF'): cl[b] = random.sample(cl[b], min(nb, len(cl[b])))
with open(os.path.join(out, 'classes.tsv'), 'w') as o, open(os.path.join(out, 'targets.tsv'), 'w') as t:
    o.write('chrom\tpos\tref\talt\tclass\n'); pos = set()
    for name, ks in cl.items():
        for k in ks: o.write('\t'.join(k) + f'\t{name}\n'); pos.add((k[0], int(k[1])))
    for p in sorted(pos, key=lambda z: z[1]): t.write(f'{p[0]}\t{p[1]}\n')
print(f"[site_classes_bayes] {d}: " + ' | '.join(f"{k} {len(v)}" for k, v in cl.items()))
