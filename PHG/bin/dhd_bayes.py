#!/usr/bin/env python3
"""Donor allele (DHd) at every biallelic union site — own sites fixed, gaps filled from new read counts (docs/PLAN_marker_union_pilot.md).
Scope (user, 2026-09-23):
  own site (the donor is in the union row's discovered_in)  -> ALT, as discovered; never re-called
  gap site (discovered only in other donors)                -> from the donor's count-once reads (joint step-4 table):
      REF  if tier 'ref'  (the donor's own reads show B73; no prior involved)
      ALT  if posterior >= --alt, where the prior is raised by the OTHER donors' ALT at that site:
             k = prior donors with the site discovered (ALT),  m = prior donors with a genotype there (discovered ALT, or gap REF)
             pi = (w*mu_d + k) / (w + m),   mu_d = tier A / (tier A + tier ref) over the donor's gap sites
             posterior = logistic(LLR_d + logit pi);  sites flagged hidepth / af_gt_half are never promoted
      missing otherwise
  Prior donors: all other donors, or with --taxa (DONOR=taxon,...) and --same-taxon only the other donors of the same taxon.
Also writes the heuristic A+B gap state (ALT = tier A, or tier B with ALT in >1 BC1 sample) for comparison.
Output: OUTDIR/dhd_bayes_<chrom>.tsv.gz  per donor: src (own/gap), tier, LLR, k, m, prior, posterior, state_bayes, state_AB
        OUTDIR/dhd_bayes_summary.tsv
Usage: dhd_bayes.py OUTDIR union.tsv.gz DONOR=joint_step4/DONOR.sites.tsv.gz [...] [--w 2] [--alt 0.95] [--taxa D1=Zx,D2=Zd --same-taxon]"""
import gzip, os, math, collections, argparse

ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('union'); ap.add_argument('tables', nargs='+')
ap.add_argument('--w', type=float, default=2.0); ap.add_argument('--alt', type=float, default=0.95)
ap.add_argument('--taxa', default=''); ap.add_argument('--same-taxon', action='store_true')
A = ap.parse_args(); out = A.outdir; os.makedirs(out, exist_ok=True)
tabs = dict(t.split('=', 1) for t in A.tables); donors = list(tabs)
taxa = dict(x.split('=') for x in A.taxa.split(',') if x)
prior_set = {d: [e for e in donors if e != d and (not A.same_taxon or taxa.get(e) == taxa.get(d))] for d in donors}
def load(p):
    t = {}
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t')
            t[(x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']])] = (x[c['tier']], int(x[c['n']]), int(x[c['n_pools_alt']]),
                                                                              float(x[c['LLR']]), x[c['flags']])
    return t
T = {d: load(p) for d, p in tabs.items()}; NONE = ('absent', 0, 0, 0.0, '.')
union = []
with gzip.open(A.union, 'rt') as f:
    f.readline()
    for l in f:
        x = l.rstrip('\n').split('\t')
        if x[6] == '0': union.append(((x[0], int(x[1]), x[2], x[3]), set(x[5].split(','))))
def gap_ref(d, k, disc): return d not in disc and T[d].get(k, NONE)[0] == 'ref'
mu = {}
for d in donors:
    a = r = 0
    for k, disc in union:
        if d in disc: continue
        t = T[d].get(k, NONE)[0]; a += t == 'A'; r += t == 'ref'
    mu[d] = (a + 0.5) / (a + r + 1)
    print(f"[dhd_bayes] {d}: mu = {mu[d]:.3f} (gap tier A {a}, gap tier ref {r}) | prior donors {prior_set[d] or 'none'} | w = {A.w} | alt >= {A.alt}")
logit = lambda p: math.log(p / (1 - p))
def logistic(z): return 1 / (1 + math.exp(-z)) if z > -700 else 0.0
S = {d: collections.Counter() for d in donors}
chrom = union[0][0][0] if union else 'chr'
with gzip.open(os.path.join(out, f'dhd_bayes_{chrom}.tsv.gz'), 'wt') as o:
    o.write('chrom\tpos\tref\talt\tdiscovered_in\t' + '\t'.join(
        f'{d}_src\t{d}_tier\t{d}_LLR\t{d}_k\t{d}_m\t{d}_prior\t{d}_posterior\t{d}_state_bayes\t{d}_state_AB' for d in donors) + '\n')
    for k, disc in union:
        cells = []
        for d in donors:
            tier, n, npa, llr, flags = T[d].get(k, NONE)
            if d in disc:
                cells += ['own', tier, f'{llr:.2f}', '.', '.', '.', '.', '1', '1']; S[d]['own_ALT'] += 1; continue
            kk = sum(1 for e in prior_set[d] if e in disc)
            mm = kk + sum(1 for e in prior_set[d] if gap_ref(e, k, disc))
            pi = min(max((A.w * mu[d] + kk) / (A.w + mm), 1e-6), 1 - 1e-6)
            bad = any(f in flags.split(',') for f in ('hidepth', 'af_gt_half'))
            post = logistic(llr + logit(pi)) if tier != 'absent' else float('nan')
            if tier == 'ref': sb = '0'
            elif tier != 'absent' and not bad and post >= A.alt: sb = '1'
            else: sb = 'NA'
            sab = '1' if tier == 'A' or (tier == 'B' and npa >= 2) else '0' if tier == 'ref' else 'NA'
            S[d][f'gap_bayes_{sb}'] += 1; S[d][f'gap_AB_{sab}'] += 1
            if sb == '1': S[d][f'gap_bayes_ALT_tier_{tier}'] += 1
            cells += ['gap', tier, f'{llr:.2f}', str(kk), str(mm), f'{pi:.3f}', 'NA' if post != post else f'{post:.4f}', sb, sab]
        o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{','.join(sorted(disc))}\t" + '\t'.join(cells) + '\n')
n = len(union)
with open(os.path.join(out, 'dhd_bayes_summary.tsv'), 'w') as o:
    ks = sorted({k for d in donors for k in S[d]}); o.write('donor\tmu\t' + '\t'.join(ks) + '\n')
    for d in donors: o.write(f"{d}\t{mu[d]:.3f}\t" + '\t'.join(str(S[d][k]) for k in ks) + '\n')
for d in donors:
    s = S[d]; g = n - s['own_ALT']
    print(f"[dhd_bayes] {d}: own ALT {s['own_ALT']} | gaps {g}: Bayes ALT {s['gap_bayes_1']} (tier A {s['gap_bayes_ALT_tier_A']}, "
          f"B {s['gap_bayes_ALT_tier_B']}, C {s['gap_bayes_ALT_tier_C']}, - {s['gap_bayes_ALT_tier_-']}) / REF {s['gap_bayes_0']} / missing {s['gap_bayes_NA']}"
          f" || heuristic A+B: ALT {s['gap_AB_1']} / REF {s['gap_AB_0']} / missing {s['gap_AB_NA']}")
