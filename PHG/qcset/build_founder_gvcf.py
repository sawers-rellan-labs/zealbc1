#!/usr/bin/env python3
"""build_founder_gvcf — QC set design B, step 3: the PERFECT founder = the founder assembly's SNPs vs B73 (anchorwave MAF -> SNP
table `chrom pos ref alt`) restricted to the lowcopy BED, written as a haploid gVCF in the pilot's founder style
(pilot_step5_donor_gvcf.py) so PHG loads it exactly like the discovered founders A / A+B.
  alt record : REF ALT,<NON_REF>  GT:AD:DP:PL 1:0,d,0:d:90,90,0
  ref block  : REF <NON_REF> END=..  GT:AD:DP:PL 0:d,0:d:0,90,90     over every BED interval, split around alt records
(PL keeps the pilot founders 3-value layout on purpose: byte-parity with the A / A+B gVCFs PHG has already loaded.)
Also writes <out>.alt.tsv (chrom pos ref alt) = the truth allele set for benchmarking discovery (false / missed founder alleles).
Usage: build_founder_gvcf.py --snps Gigi_vs_B73_chr10_snps.tsv --donor Gigi_PERFECT --bed union_chr10.bed --fai B73.fa.fai
                             --out Gigi_PERFECT.g.vcf [--chrom chr10] [--depth 30]
"""
import argparse, bisect, sys
ap = argparse.ArgumentParser()
ap.add_argument('--snps', required=True, help='TSV chrom pos ref alt (1-based), no header')
ap.add_argument('--donor', required=True); ap.add_argument('--bed', required=True); ap.add_argument('--fai', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--chrom', default='chr10')
ap.add_argument('--depth', type=int, default=30, help='nominal DP written in every record (the founder is error-free)')
ap.add_argument('--asm-chr', help='ASM_Chr value (default <donor>_<chrom>)')
A = ap.parse_args()
ASM = A.asm_chr or f"{A.donor}_{A.chrom}"; d = A.depth

fields = {l.split()[0]: l.split() for l in open(A.fai)}
if A.chrom not in fields: sys.exit(f"{A.chrom} not in {A.fai}")
_, ln, off, lb, lw = fields[A.chrom][0], int(fields[A.chrom][1]), int(fields[A.chrom][2]), int(fields[A.chrom][3]), int(fields[A.chrom][4])
fh = open(A.fai[:-4], 'rb')
def base(pos1):
    p = pos1 - 1; fh.seek(off + (p // lb) * lw + (p % lb)); return fh.read(1).decode().upper()

iv = sorted((int(x[1]), int(x[2])) for x in (l.split() for l in open(A.bed)) if x and x[0] == A.chrom)   # 0-based half-open
starts = [s for s, e in iv]
def which(pos1):
    i = bisect.bisect_right(starts, pos1 - 1) - 1
    return i if i >= 0 and pos1 - 1 < iv[i][1] else -1

alts = {}; n_all = n_out = n_dup = n_refmis = n_nonsnp = 0
for l in open(A.snps):
    x = l.split()
    if not x or x[0] != A.chrom: continue
    n_all += 1
    pos, ref, alt = int(x[1]), x[2].upper(), x[3].upper()
    if len(ref) != 1 or len(alt) != 1 or ref == alt or alt not in 'ACGT': n_nonsnp += 1; continue
    if which(pos) < 0: n_out += 1; continue
    if pos in alts: n_dup += 1; continue
    if base(pos) != ref: n_refmis += 1; continue
    alts[pos] = (ref, alt)
altpos = sorted(alts)

nblk = nalt = 0
with open(A.out, 'w') as o, open(A.out + '.alt.tsv', 'w') as t:
    o.write("##fileformat=VCFv4.2\n##source=zeal_qcset build_founder_gvcf PERFECT founder (assembly SNPs vs B73 in lowcopy ranges)\n")
    o.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n')
    for k_, t_ in (('ASM_Chr', 'String'), ('ASM_Start', 'Integer'), ('ASM_End', 'Integer'), ('ASM_Strand', 'String')):
        o.write(f'##INFO=<ID={k_},Number=1,Type={t_},Description="pseudo-assembly coordinate">\n')
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    o.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    o.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n')
    o.write(f"##contig=<ID={A.chrom},length={ln}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{A.donor}\n")
    def block(a, b):   # 1-based inclusive
        global nblk
        o.write(f"{A.chrom}\t{a}\t.\t{base(a)}\t<NON_REF>\t.\t.\tASM_Chr={ASM};ASM_Start={a};ASM_End={b};ASM_Strand=+;END={b}\tGT:AD:DP:PL\t0:{d},0:{d}:0,90,90\n"); nblk += 1
    for s, e in iv:
        cur = s + 1
        lo = bisect.bisect_right(altpos, s); hi = bisect.bisect_right(altpos, e)
        for p in altpos[lo:hi]:
            if p > cur: block(cur, p - 1)
            r, a = alts[p]
            o.write(f"{A.chrom}\t{p}\t.\t{r}\t{a},<NON_REF>\t90\t.\tASM_Chr={ASM};ASM_Start={p};ASM_End={p};ASM_Strand=+\tGT:AD:DP:PL\t1:0,{d},0:{d}:90,90,0\n"); nalt += 1
            t.write(f"{A.chrom}\t{p}\t{r}\t{a}\n"); cur = p + 1
        if cur <= e: block(cur, e)
print(f"{A.donor}: snps_in={n_all} non_snp={n_nonsnp} outside_bed={n_out} dup_pos={n_dup} ref_mismatch={n_refmis} "
      f"-> alt_records={nalt} ref_blocks={nblk} over {len(iv)} intervals ({sum(e - s for s, e in iv) / 1e6:.1f} Mb)")
