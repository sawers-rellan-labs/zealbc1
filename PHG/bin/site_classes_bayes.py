#!/usr/bin/env python3
"""Site classes for the per-site test of the count-once allele calls of one donor (scored by compare_discovery_sets.py score):
  bayes_ALT_from_BC   Bayes ALT where the joint step-4 tier is B or C (admitted by the prior)
  bayes_REF_from_unk  Bayes REF where the joint tier is not 'ref' (admitted by the prior)
  ownA_now_ref        the donor's own tier-A discovery sites that the joint step 4 calls 'ref'
  base_ALT_tierA      baseline: joint tier A (random sample)
  base_REF_tierref    baseline: joint tier ref (random sample)
Usage: site_classes_bayes.py DONOR dhd_bayes.tsv.gz dhd_joint.tsv.gz OUTDIR [n_base=5000]"""
import sys, gzip, os, random
d, fb, fj, out = sys.argv[1:5]; nb = int(sys.argv[5]) if len(sys.argv) > 5 else 5000; os.makedirs(out, exist_ok=True)
def rows(p):
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f: yield c, l.rstrip('\n').split('\t')
J = {}
for c, x in rows(fj): J[(x[0], x[1], x[2], x[3])] = (x[c[f'{d}_tier']], d in x[c['discovered_in']].split(','))
cl = {k: [] for k in ('bayes_ALT_from_BC', 'bayes_REF_from_unk', 'ownA_now_ref', 'base_ALT_tierA', 'base_REF_tierref')}
for c, x in rows(fb):
    k = (x[0], x[1], x[2], x[3]); st = x[c[f'{d}_state']]; tier, own = J[k]
    if st == '1' and tier in ('B', 'C'): cl['bayes_ALT_from_BC'].append(k)
    if st == '0' and tier != 'ref': cl['bayes_REF_from_unk'].append(k)
    if own and tier == 'ref': cl['ownA_now_ref'].append(k)
    if tier == 'A': cl['base_ALT_tierA'].append(k)
    if tier == 'ref' and not own: cl['base_REF_tierref'].append(k)   # disjoint from ownA_now_ref (one class per site)
random.seed(1)
for b in ('base_ALT_tierA', 'base_REF_tierref'): cl[b] = random.sample(cl[b], min(nb, len(cl[b])))
with open(os.path.join(out, 'classes.tsv'), 'w') as o, open(os.path.join(out, 'targets.tsv'), 'w') as t:
    o.write('chrom\tpos\tref\talt\tclass\n'); pos = set()
    for name, ks in cl.items():
        for k in ks: o.write('\t'.join(k) + f'\t{name}\n'); pos.add((k[0], int(k[1])))
    for p in sorted(pos, key=lambda z: z[1]): t.write(f'{p[0]}\t{p[1]}\n')
print(f"[site_classes_bayes] {d}: " + ' | '.join(f"{k} {len(v)}" for k, v in cl.items()))
