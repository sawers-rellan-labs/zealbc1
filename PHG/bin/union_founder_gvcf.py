#!/usr/bin/env python3
"""Union founder gVCF per donor (docs/PLAN_marker_union_pilot.md step 4) from the gap-filling table (dhd_bayes.py).
At every union site the donor's state_bayes decides the record:
  1  (own ALT, or gap ALT at posterior >= cut-off) -> alt record        REF ALT,<NON_REF>  GT:AD:DP:PL 1:0,d,0:d:90,90,0
  0  (gap REF: tier ref from the donor's own reads) -> reference record (inside a ref block, or a 1-bp block if outside the covered BED)
  NA (missing)                                       -> NO record: the ref block is split around the site
Reference blocks run over the donor's covered stretches (covered_blocks.sh) and are split around alt records and missing sites; the format
is pilot_step5_donor_gvcf.py's, so PHG loads it like the earlier founders. Positions with two union alleles (multi-allelic) are missing.
Also writes <out>.alt.tsv (chrom pos ref alt) for the pseudo-assembly (bcftools consensus).
Usage: union_founder_gvcf.py --dhd dhd_bayes_chr10.tsv.gz --donor Zx.0570_P2 --covered-bed covered.bed --fai B73.fa.fai --out D.g.vcf
                             [--chrom chr10] [--depth 8] [--name FOUNDER_NAME]"""
import argparse, gzip, bisect, collections
ap = argparse.ArgumentParser()
ap.add_argument('--dhd', required=True); ap.add_argument('--donor', required=True); ap.add_argument('--covered-bed', required=True)
ap.add_argument('--fai', required=True); ap.add_argument('--out', required=True); ap.add_argument('--chrom', default='chr10')
ap.add_argument('--depth', type=int, default=8, help='nominal DP written in the records'); ap.add_argument('--name', help='sample name (default donor)')
A = ap.parse_args(); NAME = A.name or A.donor; ASM = f"{NAME}_{A.chrom}"; d = A.depth

state = collections.defaultdict(list)            # pos -> [(ref, alt, state)]
with gzip.open(A.dhd, 'rt') as f:
    h = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(h)}; sc = c[f'{A.donor}_state_bayes']
    for l in f:
        x = l.rstrip('\n').split('\t')
        if x[c['chrom']] == A.chrom: state[int(x[c['pos']])].append((x[c['ref']], x[c['alt']], x[sc]))
alt, ref, miss = {}, set(), set()
for p, v in state.items():
    if len(v) > 1: miss.add(p)
    elif v[0][2] == '1': alt[p] = v[0][:2]
    elif v[0][2] == '0': ref.add(p)
    else: miss.add(p)

fields = {l.split()[0]: l.split() for l in open(A.fai)}
_, ln, off, lb, lw = fields[A.chrom][0], int(fields[A.chrom][1]), int(fields[A.chrom][2]), int(fields[A.chrom][3]), int(fields[A.chrom][4])
fh = open(A.fai[:-4], 'rb')
def base(p1):
    p = p1 - 1; fh.seek(off + (p // lb) * lw + (p % lb)); return fh.read(1).decode().upper()
iv = sorted((int(x[1]), int(x[2])) for x in (l.split() for l in open(A.covered_bed)) if x and x[0] == A.chrom)   # 0-based half-open
starts = [s for s, e in iv]
def inside(p1):
    i = bisect.bisect_right(starts, p1 - 1) - 1
    return i >= 0 and p1 - 1 < iv[i][1]
# REF sites outside the covered stretches still get a record (1-bp interval); ALT sites outside them get their own interval too
extra = sorted(p for p in list(ref) + list(alt) if not inside(p))
iv = sorted(iv + [(p - 1, p) for p in extra])

n = collections.Counter()
def block(o, s1, e1):
    o.write(f"{A.chrom}\t{s1}\t.\t{base(s1)}\t<NON_REF>\t.\t.\tASM_Chr={ASM};ASM_Start={s1};ASM_End={e1};ASM_Strand=+;END={e1}\tGT:AD:DP:PL\t0:{d},0:{d}:0,90,90\n"); n['ref_blocks'] += 1
events = sorted(set(alt) | miss)
with open(A.out, 'w') as o, open(A.out + '.alt.tsv', 'w') as t:
    o.write("##fileformat=VCFv4.2\n##source=zeal union_founder_gvcf (dhd_bayes state: 1 alt record, 0 reference, NA no record)\n")
    o.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n')
    for k_, t_ in (('ASM_Chr', 'String'), ('ASM_Start', 'Integer'), ('ASM_End', 'Integer'), ('ASM_Strand', 'String')):
        o.write(f'##INFO=<ID={k_},Number=1,Type={t_},Description="pseudo-assembly coordinate">\n')
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=AD,Number=3,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    o.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n')
    o.write(f"##contig=<ID={A.chrom},length={ln}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{NAME}\n")
    for s, e in iv:
        cur = s + 1
        lo, hi = bisect.bisect_right(events, s), bisect.bisect_right(events, e)
        for p in events[lo:hi]:
            if p > cur: block(o, cur, p - 1)
            if p in alt:
                r, a = alt[p]
                o.write(f"{A.chrom}\t{p}\t.\t{r}\t{a},<NON_REF>\t60\t.\tASM_Chr={ASM};ASM_Start={p};ASM_End={p};ASM_Strand=+\tGT:AD:DP:PL\t1:0,{d},0:{d}:90,90,0\n")
                t.write(f"{A.chrom}\t{p}\t{r}\t{a}\n"); n['alt'] += 1
            else: n['missing_holes'] += 1
            cur = p + 1
        if cur <= e: block(o, cur, e)
print(f"[union_founder] {NAME}: union positions {len(state)} | ALT {len(alt)} (records {n['alt']}) | REF {len(ref)} | missing {len(miss)} "
      f"(holes in blocks {n['missing_holes']}) | covered intervals {len(starts)} + {len(extra)} 1-bp | ref blocks {n['ref_blocks']}")
