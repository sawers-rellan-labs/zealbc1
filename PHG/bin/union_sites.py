#!/usr/bin/env python3
"""Union of tier-A sites over donors (docs/PLAN_marker_union_pilot.md step 1).
Input: per-donor step-4 tables (<donor>.sites.tsv.gz). Key = (chrom, pos, ref, alt).
Output (OUTDIR): union_<chrom>.tsv.gz  chrom pos ref alt n_donors donors (tier-A carriers) multiallelic(0/1)
                 per_donor.tsv   donor | tierA | shared (tier A in >=1 other donor) | private | other-only union sites |
                                 of those: in own step-4 table (has counts, any tier) | absent (needs the second counting pass)
                 need_counts_<donor>.tsv  chrom pos (union sites absent from the donor's own step-4 table)
Usage: union_sites.py OUTDIR DONOR=path/to/DONOR.sites.tsv.gz [DONOR=...]"""
import sys, gzip, os, collections

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
tabs = dict(a.split('=', 1) for a in sys.argv[2:])
tierA, present = {}, {}
for d, p in tabs.items():
    A, P = set(), set()
    with gzip.open(p, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t'); k = (x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']])
            P.add(k[:2])
            if x[c['tier']] == 'A': A.add(k)
    tierA[d], present[d] = A, P
    print(f"[union] {d}: tier A {len(A)} | positions in step-4 table {len(P)}")
carriers = collections.defaultdict(list)
for d, A in tierA.items():
    for k in A: carriers[k].append(d)
alts_at = collections.defaultdict(set)
for k in carriers: alts_at[k[:2]].add(k[3])
multi = {p for p, a in alts_at.items() if len(a) > 1}
chroms = sorted({k[0] for k in carriers})
for ch in chroms:
    with gzip.open(os.path.join(out, f'union_{ch}.tsv.gz'), 'wt') as o:
        o.write('chrom\tpos\tref\talt\tn_donors\tdonors\tmultiallelic\n')
        for k in sorted((k for k in carriers if k[0] == ch), key=lambda z: (z[1], z[3])):
            o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{len(carriers[k])}\t{','.join(sorted(carriers[k]))}\t{int(k[:2] in multi)}\n")
upos = {k[:2] for k in carriers}
print(f"[union] union alleles {len(carriers)} | positions {len(upos)} | multi-allelic positions {len(multi)}")
nd = collections.Counter(len(v) for v in carriers.values())
print('[union] alleles by number of tier-A donors: ' + ', '.join(f"{n}: {nd[n]}" for n in sorted(nd)))
with open(os.path.join(out, 'per_donor.tsv'), 'w') as o:
    o.write('donor\ttierA\tshared\tprivate\tother_only_positions\tother_only_in_own_table\tother_only_need_counts\n')
    for d in tabs:
        A = tierA[d]; sh = sum(1 for k in A if len(carriers[k]) > 1)
        own = {k[:2] for k in A}; other = upos - own
        inown = sum(1 for p in other if p in present[d]); need = sorted(p for p in other if p not in present[d])
        with open(os.path.join(out, f'need_counts_{d}.tsv'), 'w') as n:
            for p in need: n.write(f"{p[0]}\t{p[1]}\n")
        row = (d, len(A), sh, len(A) - sh, len(other), inown, len(need))
        o.write('\t'.join(map(str, row)) + '\n')
        print(f"[union] {d}: tier A {len(A)} | shared {sh} | private {len(A)-sh} | other-only positions {len(other)} "
              f"| already in own step-4 table {inown} | need second pass {len(need)}")
