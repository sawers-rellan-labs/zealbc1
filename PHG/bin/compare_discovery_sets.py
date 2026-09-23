#!/usr/bin/env python3
"""Compare two step-4 tier-A site sets of one donor (e.g. an earlier and a rerun discovery) — checks 1-3 of the 2026-09-23 review.
  prep : overlap of tier-A alleles (key chrom,pos,ref,alt) -> shared / lost (old only) / gained (new only); writes the pileup targets
         (all lost + gained, plus a random sample of shared as the baseline).
  score: per-site test with reads from the donor's BC2S3 lines, classified by an INDEPENDENT ancestry track (RTIGER SNP50K segments):
         alt fraction inside the line's donor segments (state 1/2) vs inside its B73 segments (state 0). A real donor allele shows alt
         reads inside donor segments and ~none in B73 segments; a false allele shows ~none inside, or alt everywhere (artifact).
         Also: B73 artifact sites = tier-A sites where the run's B73 'donor' row is tier A/B/C.
Usage: compare_discovery_sets.py prep  OLD.sites.tsv.gz NEW.sites.tsv.gz OUTDIR [n_shared=5000]
       compare_discovery_sets.py score OUTDIR ad.tsv rtiger_segments.csv OLD_B73.sites.tsv.gz NEW_B73.sites.tsv.gz"""
import sys, gzip, csv, random, re, os, bisect, collections, statistics

def tierA(path):
    s = {}
    with gzip.open(path, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[c['tier']] == 'A': s[(x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']])] = x
    return s

def b73_flagged(path):
    s = set()
    with gzip.open(path, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[c['tier']] in ('A', 'B', 'C'): s.add((x[c['chrom']], int(x[c['pos']]), x[c['ref']], x[c['alt']]))
    return s

mode = sys.argv[1]
if mode == 'prep':
    old, new, out = tierA(sys.argv[2]), tierA(sys.argv[3]), sys.argv[4]; nsh = int(sys.argv[5]) if len(sys.argv) > 5 else 5000
    os.makedirs(out, exist_ok=True)
    sh = sorted(set(old) & set(new)); lost = sorted(set(old) - set(new)); gained = sorted(set(new) - set(old))
    random.seed(1); shs = sorted(random.sample(sh, min(nsh, len(sh))))
    print(f"[check 1] tier A old {len(old)} | new {len(new)} | shared {len(sh)} | lost (old only) {len(lost)} | gained (new only) {len(gained)}")
    lostpos = {(k[0], k[1]) for k in lost}; gainpos = {(k[0], k[1]) for k in gained}
    print(f"[check 1] lost/gained at the same position with a different allele: {len(lostpos & gainpos)}")
    with open(os.path.join(out, 'classes.tsv'), 'w') as o:
        o.write('chrom\tpos\tref\talt\tclass\n')
        for cl, ks in (('shared', shs), ('lost', lost), ('gained', gained)):
            for k in ks: o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{cl}\n")
    with open(os.path.join(out, 'targets.tsv'), 'w') as o:
        for k in sorted({(k[0], k[1]) for k in shs + lost + gained}, key=lambda z: z[1]): o.write(f"{k[0]}\t{k[1]}\n")
    with open(os.path.join(out, 'overlap.tsv'), 'w') as o:
        o.write(f"old\tnew\tshared\tlost\tgained\n{len(old)}\t{len(new)}\t{len(sh)}\t{len(lost)}\t{len(gained)}\n")
    sys.exit(0)

out, adf, rtf, ob73, nb73 = sys.argv[2:7]
cls = {}
with open(os.path.join(out, 'classes.tsv')) as f:
    next(f)
    for l in f:
        c, p, r, a, k = l.rstrip('\n').split('\t'); cls[(c, int(p), r, a)] = k
# independent ancestry: RTIGER SNP50K segments per line
seg = collections.defaultdict(list)
with open(rtf) as f:
    for r in csv.DictReader(f): seg[r['name']].append((int(r['start_bp']), int(r['end_bp']), int(r['state'])))
for k in seg: seg[k].sort()
starts = {k: [s[0] for s in v] for k, v in seg.items()}
def state(line, pos):
    if line not in seg: return None
    i = bisect.bisect_right(starts[line], pos) - 1
    if i < 0: return None
    s, e, st = seg[line][i]
    return st if pos <= e else None
acc = {}   # key -> [alt_in, n_in, alt_out, n_out]
samples = None; nolines = set()
with open(adf) as f:
    for l in f:
        x = l.rstrip('\n').split('\t')
        if l.startswith('#'):
            samples = [os.path.basename(re.sub(r'^\[\d+\]', '', c).rsplit(':', 1)[0]).replace('.cram', '') for c in x[4:]]
            nolines = {s for s in samples if s not in seg}; continue
        c, p, r = x[0], int(x[1]), x[2]; alts = [] if x[3] == '.' else x[3].split(',')
        cand = [(c, p, r, a) for a in 'ACGT' if (c, p, r, a) in cls]
        for key in cand:
            ai = alts.index(key[3]) + 1 if key[3] in alts else None
            v = acc.setdefault(key, [0, 0, 0, 0])
            for smp, ad in zip(samples, x[4:]):
                if ad in ('.', ''): continue
                cnt = [int(z) if z != '.' else 0 for z in ad.split(',')]; rc = cnt[0]; ac = cnt[ai] if ai is not None and ai < len(cnt) else 0
                st = state(smp, p)
                if st is None or rc + ac == 0: continue
                if st > 0: v[0] += ac; v[1] += rc + ac
                else: v[2] += ac; v[3] += rc + ac
if nolines: print(f"[check 2] samples without RTIGER segments (ignored): {sorted(nolines)}")
fo, fn = b73_flagged(ob73), b73_flagged(nb73)
rows = []
for cl in ('shared', 'lost', 'gained'):
    ks = [k for k, v in cls.items() if v == cl]
    vs = [acc[k] for k in ks if k in acc]
    ins = [v[0] / v[1] for v in vs if v[1] >= 3]; outs = [v[2] / v[3] for v in vs if v[3] >= 3]
    any_in = sum(1 for v in vs if v[1] >= 3 and v[0] > 0); zero_in = sum(1 for v in vs if v[1] >= 3 and v[0] == 0)
    b_old = sum(1 for k in ks if k in fo); b_new = sum(1 for k in ks if k in fn)
    md = lambda z: f"{statistics.median(z):.3f}" if z else 'NA'
    mn = lambda z: f"{statistics.mean(z):.3f}" if z else 'NA'
    rows.append((cl, len(ks), len(vs), len(ins), md(ins), mn(ins), any_in, zero_in, len(outs), mn(outs), b_old, b_new))
hdr = ('class', 'sites', 'with_reads', 'n_in>=3', 'median_altfrac_in_donor_seg', 'mean_altfrac_in_donor_seg', 'alt_seen_in_seg',
       'no_alt_in_seg', 'n_out>=3', 'mean_altfrac_in_B73_seg', 'B73_flag_oldrun', 'B73_flag_newrun')
with open(os.path.join(out, 'site_test.tsv'), 'w') as o:
    o.write('\t'.join(hdr) + '\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
print('[check 2+3] ' + '\t'.join(hdr))
for r in rows: print('[check 2+3] ' + '\t'.join(map(str, r)))
