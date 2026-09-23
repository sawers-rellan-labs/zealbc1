#!/usr/bin/env python3
"""Empirical-Bayes donor allele at every biallelic union site, from the joint step-4 tables (count-once).
Step 4 gives each donor d at site s a log Bayes factor LLR(d,s) (donor carries the ALT vs sequencing error) and uses a flat prior 0.5.
Here the prior comes from the OTHER donors' calls at the same site (leave-one-out), Beta-binomial:
    k = other donors with tier A at s,  m = other donors with tier A or tier ref at s
    pi(d,s) = (alpha_d + k) / (alpha_d + beta_d + m),   alpha_d = w*mu_d, beta_d = w*(1-mu_d)
    mu_d = donor d's sharing rate: tier A / (tier A + tier ref) at the union sites discovered only by other donors
    posterior(d,s) = logistic(LLR + logit pi)
State: ALT if posterior >= --alt (0.95); REF if posterior <= --ref (0.05) and n >= 12; unknown otherwise. Sites flagged hidepth or
af_gt_half (paralog/CNV signatures) stay unknown whatever the prior. Also cross-tabulates against the fixed rules of dhd_joint.py.
Output: OUTDIR/dhd_bayes_<chrom>.tsv.gz, OUTDIR/dhd_bayes_summary.tsv
Usage: dhd_bayes.py OUTDIR union.tsv.gz DONOR=step4/DONOR.sites.tsv.gz [DONOR=...] [--w 2] [--alt 0.95] [--ref 0.05]"""
import sys, gzip, os, math, collections

import argparse
ap = argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('union'); ap.add_argument('tables', nargs='+')
ap.add_argument('--w', type=float, default=2.0); ap.add_argument('--alt', type=float, default=0.95); ap.add_argument('--ref', type=float, default=0.05)
A = ap.parse_args(); W, PALT, PREF = A.w, A.alt, A.ref
out, uf = A.outdir, A.union; tabs = dict(t.split('=', 1) for t in A.tables)
def load(p):
    t = {}
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t')
            t[(x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']])] = (x[c['tier']], int(x[c['n']]), int(x[c['n_pools_alt']]),
                                                                              float(x[c['LLR']]), x[c['flags']])
    return t
T = {d: load(p) for d, p in tabs.items()}; donors = list(tabs)
union = []
with gzip.open(uf, 'rt') as f:
    f.readline()
    for l in f:
        x = l.rstrip('\n').split('\t')
        if x[6] == '0': union.append(((x[0], int(x[1]), x[2], x[3]), set(x[5].split(','))))
NONE = ('absent', 0, 0, 0.0, '.')
mu = {}
for d in donors:
    a = r = 0
    for k, disc in union:
        if d in disc: continue
        t = T[d].get(k, NONE)[0]; a += t == 'A'; r += t == 'ref'
    mu[d] = (a + 0.5) / (a + r + 1); print(f"[dhd_bayes] {d}: sharing rate mu = {mu[d]:.3f} (tier A {a}, tier ref {r} at other donors' sites) | w = {W}")
logit = lambda p: math.log(p / (1 - p))
def logistic(z): return 1 / (1 + math.exp(-z)) if z > -700 else 0.0
S = {d: collections.Counter() for d in donors}; X = {d: collections.Counter() for d in donors}
chrom = union[0][0][0] if union else 'chr'
with gzip.open(os.path.join(out, f'dhd_bayes_{chrom}.tsv.gz'), 'wt') as o:
    o.write('chrom\tpos\tref\talt\tdiscovered_in\t' + '\t'.join(f'{d}_LLR\t{d}_k\t{d}_m\t{d}_prior\t{d}_posterior\t{d}_state' for d in donors) + '\n')
    for k, disc in union:
        cells = []
        for d in donors:
            tier, n, npa, llr, flags = T[d].get(k, NONE)
            kk = sum(1 for e in donors if e != d and T[e].get(k, NONE)[0] == 'A')
            mm = sum(1 for e in donors if e != d and T[e].get(k, NONE)[0] in ('A', 'ref'))
            pi = (W * mu[d] + kk) / (W + mm); pi = min(max(pi, 1e-6), 1 - 1e-6)
            post = logistic(llr + logit(pi)) if tier != 'absent' else pi
            bad = any(f in flags.split(',') for f in ('hidepth', 'af_gt_half'))
            st = 'NA' if tier == 'absent' or bad else '1' if post >= PALT else '0' if (post <= PREF and n >= 12) else 'NA'
            fixA = '1' if tier == 'A' else '0' if tier == 'ref' else 'NA'
            fixAB = '1' if tier == 'A' or (tier == 'B' and npa >= 2) else '0' if tier == 'ref' else 'NA'
            S[d][st] += 1; X[d][(fixA, st)] += 1; X[d][('AB' + fixAB, st)] += 1
            if st == '1' and tier != 'A': S[d][f'ALT_from_tier_{tier}'] += 1
            cells += [f'{llr:.2f}', str(kk), str(mm), f'{pi:.3f}', f'{post:.4f}', st]
        o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{','.join(sorted(disc))}\t" + '\t'.join(cells) + '\n')
n = len(union); lab = {'1': 'ALT', '0': 'REF', 'NA': 'unknown'}
with open(os.path.join(out, 'dhd_bayes_summary.tsv'), 'w') as o:
    o.write('donor\tmu\tALT\tREF\tunknown\tALT_not_tierA\n')
    for d in donors:
        o.write(f"{d}\t{mu[d]:.3f}\t{S[d]['1']}\t{S[d]['0']}\t{S[d]['NA']}\t{S[d]['1'] - X[d][('1', '1')]}\n")
for d in donors:
    s = S[d]; pct = lambda v: f"{v} ({100 * v / n:.1f}%)"
    print(f"[dhd_bayes] {d}: ALT {pct(s['1'])} | REF {pct(s['0'])} | unknown {pct(s['NA'])} | ALT from tier " +
          ', '.join(f"{t} {s[f'ALT_from_tier_{t}']}" for t in ('B', 'C', '-', 'ref')))
    for rule, pre in (('rule A', ''), ('rule A+B', 'AB')):
        print(f"[dhd_bayes] {d} {rule} vs Bayes (rows fixed rule, cols Bayes ALT/REF/unknown): " +
              ' | '.join(f"{lab[f]}: " + '/'.join(str(X[d][(pre + f, b)]) for b in ('1', '0', 'NA')) for f in ('1', '0', 'NA')))
