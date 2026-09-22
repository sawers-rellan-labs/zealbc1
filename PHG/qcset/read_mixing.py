#!/usr/bin/env python3
"""read_mixing — QC set design B, step 4. Build every simulated sample of ONE founder from pre-aligned read sources by
per-tract mixing, in ONE walk per source file, with every read assigned to at most one destination (disjoint samples).

Destinations (per founder):
  BC1 pools   : pool<p>            depth --pool-depth (15x); per tract donor share k/12 (bc1_pool<p>_dosage.bed), B73 share 1-k/12
  sweep       : <line>_lam<λ>      depth λ (one per --lambdas); per tract state 0 -> B73, 1 -> donor:B73 1:1, 2 -> donor
  (witness = samtools merge of the sweep samples, B73 control = the real CRAM: both done by the sbatch, not here)
Sources: --donor-cram (one or more files, nominal depth --donor-depth EACH) and --b73-cram (one or more, --b73-depth EACH).
Assignment: u = hash(qname) in [0,1); at the fragment's position x the destinations own consecutive sub-intervals of [0,1)
of width target_depth_i(x) / total_nominal_depth(source); the read goes to the owner of u, or nowhere. Both mates share the
qname and the fragment start min(pos, mpos), so a pair is never split. Reads are kept only if they overlap the lowcopy BED
± --flank. Output: <out>/<dest>.<source>.bam (coordinate order preserved) + <out>/samples.tsv + <out>/demand.tsv.
--plan-only computes the demand table (max Σ target / nominal per source) from the breakpoint tables and exits: run it
BEFORE aligning to size the sources (demand > 1 = oversubscribed, the walk refuses to run).
"""
import argparse, bisect, csv, hashlib, os, sys, time
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument('--founder', required=True)
ap.add_argument('--breakpoints', required=True, help='breakpoint_sim output dir')
ap.add_argument('--bed', required=True, help='lowcopy ranges BED'); ap.add_argument('--flank', type=int, default=1000)
ap.add_argument('--chrom', default='chr10'); ap.add_argument('--chrlen', type=int, default=152435371)
ap.add_argument('--lambdas', default='0.05,0.1,0.2,0.4,0.8,1.2'); ap.add_argument('--pool-depth', type=float, default=15.0)
ap.add_argument('--pool-n', type=int, default=12, help='haploid genomes per pool (6 plants)')
ap.add_argument('--donor-cram', nargs='*', default=[]); ap.add_argument('--donor-depth', type=float, default=22.0)
ap.add_argument('--b73-cram', nargs='*', default=[]); ap.add_argument('--b73-depth', type=float, default=20.0)
ap.add_argument('--reference', help='B73 FASTA (CRAM decoding)')
ap.add_argument('--out', required=True); ap.add_argument('--only', help='comma list of destinations to write (unit test)')
ap.add_argument('--plan-only', action='store_true'); ap.add_argument('--salt', default='qcsetB')
A = ap.parse_args()
lambdas = [float(x) for x in A.lambdas.split(',')]
CH = A.chrom
chrnum = CH.replace('chr', '')

def read_tsv(path):
    with open(path) as fh: return list(csv.DictReader(fh, delimiter='\t'))

# --- destinations: list of (name, kind, tracts) with tracts = sorted [(start0, end0, donor_share)] covering the chromosome
dests = []
seg = [r for r in read_tsv(os.path.join(A.breakpoints, 'bc2s3_truth_segments.tsv')) if r['chr'] in (CH, chrnum)]
lines = sorted({r['name'] for r in seg})
for ln in lines:
    tr = sorted((int(r['start_bp']) - 1, int(r['end_bp']), int(r['state']) / 2.0) for r in seg if r['name'] == ln)
    # segments carry marker coordinates: extend the first/last to the chromosome ends so every read has a tract
    tr[0] = (0, tr[0][1], tr[0][2]); tr[-1] = (tr[-1][0], A.chrlen, tr[-1][2])
    for lam in lambdas: dests.append((f"{A.founder}_{ln}_lam{lam:g}", 'sweep', lam, tr))
p = 1
while os.path.exists(os.path.join(A.breakpoints, f'bc1_pool{p}_dosage.bed')):
    tr = sorted((int(x[1]), int(x[2]), int(x[3]) / A.pool_n) for x in (l.split() for l in open(os.path.join(A.breakpoints, f'bc1_pool{p}_dosage.bed'))) if x[0] == CH)
    dests.append((f"{A.founder}_pool{p}", 'pool', A.pool_depth, tr)); p += 1
if not lines or p == 1: sys.exit(f"no BC2S3 segments / pool BEDs for {CH} in {A.breakpoints}")
if A.only:
    keep = set(A.only.split(',')); dests = [d for d in dests if d[0] in keep]
    if not dests: sys.exit(f"--only matched nothing")

def share_at(tr, x):
    i = bisect.bisect_right([t[0] for t in tr], x) - 1
    return tr[max(i, 0)][2]

# --- piecewise-constant demand: cut points = union of all tract boundaries
cuts = sorted({0, A.chrlen} | {t[0] for d in dests for t in d[3]} | {t[1] for d in dests for t in d[3]})
cuts = [c for c in cuts if 0 <= c <= A.chrlen]
nsrc = {'donor': max(1, len(A.donor_cram)) * A.donor_depth, 'b73': max(1, len(A.b73_cram)) * A.b73_depth}
# cum[src][interval i] = list of cumulative fractions (len = ndest) in dest order
cum = {'donor': [], 'b73': []}; demand = []
for i in range(len(cuts) - 1):
    x = cuts[i]; row = {'donor': [], 'b73': []}; acc = {'donor': 0.0, 'b73': 0.0}
    for name, kind, depth, tr in dests:
        s = share_at(tr, x)
        for src, want in (('donor', depth * s), ('b73', depth * (1 - s))):
            acc[src] += want / nsrc[src]; row[src].append(acc[src])
    for src in row: cum[src].append(row[src])
    demand.append((cuts[i], cuts[i + 1], acc['donor'] * nsrc['donor'], acc['b73'] * nsrc['b73']))
os.makedirs(A.out, exist_ok=True)
with open(os.path.join(A.out, 'demand.tsv'), 'w') as fh:
    fh.write("chrom\tstart\tend\tdonor_depth_needed\tb73_depth_needed\n")
    for s, e, dd, bd in demand: fh.write(f"{CH}\t{s}\t{e}\t{dd:.2f}\t{bd:.2f}\n")
with open(os.path.join(A.out, 'samples.tsv'), 'w') as fh:
    fh.write("sample\tkind\tfounder\tline_or_pool\tlambda_or_depth\n")
    for name, kind, depth, tr in dests:
        fh.write(f"{name}\t{kind}\t{A.founder}\t{name[len(A.founder) + 1:].split('_lam')[0]}\t{depth:g}\n")
maxd = max(d[2] for d in demand); maxb = max(d[3] for d in demand)
wmax_d = sum((e - s) for s, e, dd, bd in demand if dd > nsrc['donor']) / A.chrlen
wmax_b = sum((e - s) for s, e, dd, bd in demand if bd > nsrc['b73']) / A.chrlen
print(f"[read_mixing] {A.founder}: {len(dests)} destinations ({len(lines)} lines x {len(lambdas)} lambdas + {p - 1} pools), {len(cuts) - 1} intervals")
print(f"[read_mixing] donor demand max {maxd:.1f}x vs {nsrc['donor']:.0f}x available ({wmax_d * 100:.1f}% of chr oversubscribed); "
      f"B73 demand max {maxb:.1f}x vs {nsrc['b73']:.0f}x available ({wmax_b * 100:.1f}% oversubscribed)")
if A.plan_only: sys.exit(0)
if maxd > nsrc['donor'] or maxb > nsrc['b73']: sys.exit("[read_mixing] oversubscribed: add source slices (see demand.tsv)")

# --- the walk
import pysam
bed = sorted((max(0, int(x[1]) - A.flank), int(x[2]) + A.flank) for x in (l.split() for l in open(A.bed)) if x and x[0] == CH)
merged = []
for s, e in bed:
    if merged and s <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    else: merged.append((s, e))
bstart = [s for s, e in merged]
def in_lowcopy(rs, re_):
    i = bisect.bisect_right(bstart, re_ - 1) - 1
    return i >= 0 and merged[i][1] > rs and merged[i][0] < re_
cstart = cuts[:-1]
def uniform(qname, src):
    h = hashlib.blake2b(f"{A.salt}|{src}|{qname}".encode(), digest_size=8).digest()
    return int.from_bytes(h, 'little') / 2 ** 64
names = [d[0] for d in dests]
def walk(path, src):
    tag = os.path.basename(path).rsplit('.', 1)[0]
    t0 = time.time(); n = kept = 0; counts = defaultdict(int)
    with pysam.AlignmentFile(path, reference_filename=A.reference) as fin:
        outs = {nm: pysam.AlignmentFile(os.path.join(A.out, f"{nm}.{tag}.bam"), 'wb', template=fin) for nm in names}
        for r in fin.fetch(CH):
            n += 1
            if n % 5_000_000 == 0:
                el = (time.time() - t0) / 60; print(f"[read_mixing] {tag}: {n / 1e6:.0f} M reads | kept {kept / 1e6:.2f} M | {el:.1f} min", flush=True)
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            rs, re_ = r.reference_start, r.reference_end or r.reference_start + 1
            if not in_lowcopy(rs, re_): continue
            x = rs if (r.mate_is_unmapped or r.next_reference_id != r.reference_id) else min(rs, r.next_reference_start)
            i = bisect.bisect_right(cstart, x) - 1
            u = uniform(r.query_name, src); c = cum[src][i]
            j = bisect.bisect_right(c, u)
            if j < len(c):
                outs[names[j]].write(r); kept += 1; counts[names[j]] += 1
        for o in outs.values(): o.close()
    print(f"[read_mixing] {tag} ({src}): {n} reads walked, {kept} assigned in {(time.time() - t0) / 60:.1f} min; "
          + " ".join(f"{k}={v}" for k, v in sorted(counts.items())[:8]) + (" ..." if len(counts) > 8 else ""), flush=True)
    with open(os.path.join(A.out, 'assigned_counts.tsv'), 'a') as fh:
        for k, v in sorted(counts.items()): fh.write(f"{tag}\t{src}\t{k}\t{v}\n")
for f in A.donor_cram: walk(f, 'donor')
for f in A.b73_cram: walk(f, 'b73')
print(f"[read_mixing] done -> {A.out}")
