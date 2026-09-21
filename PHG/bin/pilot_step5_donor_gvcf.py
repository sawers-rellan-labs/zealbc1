#!/usr/bin/env python3
"""Pilot step 5 — synthetic double-haploid founder gVCF per donor from step-4 site tables (PHG create-maf-vcf gVCF style).

Usage: pilot_step5_donor_gvcf.py --sites <donor>.sites.tsv.gz --donor NAME --bed genes.bed --fai B73.fa.fai --out <donor>.g.vcf
         [--tiers A,B] [--min-depth 8] [--chrom chr10]
Haploid GT as in PHG's assembly gVCFs: alt record  REF ALT,<NON_REF>  GT:AD:DP:PL 1:0,d,0:d:90,90,0
                                         ref block   REF <NON_REF> END=.. GT:AD:DP:PL 0:d,0:d:0,90,90
Ref blocks: emitted over each gene interval where the donor has >= 1 record with depth >= --min-depth (coverage proxy),
split around alt records. Intervals with no such evidence get no record (no haplotype for this founder there).
QUAL of alt records = 10*LLR/ln(10) (phred-like). Sites of tiers not selected are neither alt nor excluded from blocks.
Also writes <out>.alt.tsv (chrom pos ref alt LLR tier) for the consensus step.
"""
import gzip, argparse, math, bisect
ap = argparse.ArgumentParser()
ap.add_argument('--sites', required=True); ap.add_argument('--donor', required=True); ap.add_argument('--bed', required=True)
ap.add_argument('--fai', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--tiers', default='A,B'); ap.add_argument('--min-depth', type=int, default=8); ap.add_argument('--chrom', default='chr10')
ap.add_argument('--covered-bed', help='depth-based covered stretches (from pilot_step5_covered_blocks.sh); if given, reference blocks are emitted over these instead of the record-coverage proxy')
ap.add_argument('--asm-chr', help='ASM_Chr value = contig name of the pseudo-assembly (default <donor>_<chrom>)')
ap.add_argument('--min-pools-alt', type=int, default=2,
                help='founder allele needs alt reads in at least this many of the donor\'s pools (user rule 2026-09-19: >1 pool; per-site test: 1-pool alleles ~50%% contradicted, multi-pool ~15%%)')
ap.add_argument('--accession-sites', nargs='*', default=[],
                help='step-4 site tables of OTHER donors of the same accession; an allele below --min-pools-alt is still kept if any of them shows alt reads at the site (fallback for 1-2 sample donors)')
A = ap.parse_args(); tiers = set(A.tiers.split(','))
acc_alt = set()   # (pos, alt) with alt reads in a same-accession donor
for f_ in A.accession_sites:
    with gzip.open(f_, 'rt') as fh:
        h_ = fh.readline().rstrip('\n').split('\t'); ci = {k: i for i, k in enumerate(h_)}
        for l in fh:
            x = l.rstrip('\n').split('\t')
            if x[ci['chrom']] == A.chrom and int(x[ci['a']]) > 0: acc_alt.add((int(x[ci['pos']]), x[ci['alt']]))
n_drop = n_acc = 0
ASM = A.asm_chr or f"{A.donor}_{A.chrom}"
length = {l.split()[0]: int(l.split()[1]) for l in open(A.fai)}
iv = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(A.covered_bed or A.bed) if l.split()[0] == A.chrom)  # 0-based half-open
alts = {}; cov = [False]*len(iv); dep = [0]*len(iv); starts = [s for s, e in iv]
def which(pos1):
    i = bisect.bisect_right(starts, pos1-1) - 1
    return i if i >= 0 and pos1-1 < iv[i][1] else -1
with gzip.open(A.sites, 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t'); c = {k: i for i, k in enumerate(hdr)}
    for l in f:
        x = l.rstrip('\n').split('\t')
        if x[c['chrom']] != A.chrom: continue
        pos = int(x[c['pos']]); n = int(x[c['n']]); i = which(pos)
        if i >= 0 and (A.covered_bed or n >= A.min_depth): cov[i] = True; dep[i] = max(dep[i], n)
        if x[c['tier']] in tiers:
            npa = int(x[c['n_pools_alt']])
            if npa < A.min_pools_alt:
                if (pos, x[c['alt']]) in acc_alt: n_acc += 1
                else: n_drop += 1; continue
            alts[pos] = (x[c['ref']], x[c['alt']], float(x[c['LLR']]), n, int(x[c['a']]), x[c['tier']])
with open(A.out, 'w') as o, open(A.out + '.alt.tsv', 'w') as t:
    o.write("##fileformat=VCFv4.2\n##source=zeal_pilot_step5 synthetic double-haploid founder from pooled BC1 counts\n")
    o.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n')
    o.write('##INFO=<ID=TIER,Number=1,Type=String,Description="step-4 confidence tier">\n')
    for k_, t_ in (('ASM_Chr','String'),('ASM_Start','Integer'),('ASM_End','Integer'),('ASM_Strand','String')):
        o.write(f'##INFO=<ID={k_},Number=1,Type={t_},Description="pseudo-assembly coordinate">\n')
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=AD,Number=3,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    o.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n')
    o.write(f"##contig=<ID={A.chrom},length={length[A.chrom]}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{A.donor}\n")
    # reference base for block starts is needed: read from fasta? Not available here -> use 'N' placeholder is invalid for PHG;
    # instead take REF from the FASTA via the .fai-indexed file next to it.
    fa = A.fai[:-4]; import os
    fh = open(fa, 'rb'); fields = {l.split()[0]: l.split() for l in open(A.fai)}
    name, ln, off, lb, lw = fields[A.chrom][0], int(fields[A.chrom][1]), int(fields[A.chrom][2]), int(fields[A.chrom][3]), int(fields[A.chrom][4])
    def base(pos1):
        p = pos1 - 1; fh.seek(off + (p // lb) * lw + (p % lb)); return fh.read(1).decode().upper()
    nblk = nalt = 0; altpos = sorted(alts)
    for k, (s, e) in enumerate(iv):
        if A.covered_bed: cov[k] = True; dep[k] = max(dep[k], A.min_depth)
        if not cov[k]: continue
        cur = s + 1  # 1-based
        for p in [q for q in altpos if s < q <= e]:
            if p > cur:
                o.write(f"{A.chrom}\t{cur}\t.\t{base(cur)}\t<NON_REF>\t.\t.\tASM_Chr={ASM};ASM_Start={cur};ASM_End={p-1};ASM_Strand=+;END={p-1}\tGT:AD:DP:PL\t0:{dep[k]},0:{dep[k]}:0,90,90\n"); nblk += 1
            r, a, llr, n, ac, tier = alts[p]; q = max(1, int(round(10*llr/math.log(10))))
            o.write(f"{A.chrom}\t{p}\t.\t{r}\t{a},<NON_REF>\t{q}\t.\tASM_Chr={ASM};ASM_Start={p};ASM_End={p};ASM_Strand=+;TIER={tier}\tGT:AD:DP:PL\t1:0,{ac},0:{n}:90,90,0\n"); nalt += 1
            t.write(f"{A.chrom}\t{p}\t{r}\t{a}\t{llr:.2f}\t{tier}\n"); cur = p + 1
        if cur <= e:
            o.write(f"{A.chrom}\t{cur}\t.\t{base(cur)}\t<NON_REF>\t.\t.\tASM_Chr={ASM};ASM_Start={cur};ASM_End={e};ASM_Strand=+;END={e}\tGT:AD:DP:PL\t0:{dep[k]},0:{dep[k]}:0,90,90\n"); nblk += 1
print(f"{A.donor}: intervals={len(iv)} covered={sum(cov)} ref_blocks={nblk} alt_records={nalt} (tiers {A.tiers}; min_pools_alt={A.min_pools_alt}: dropped {n_drop} single-pool alleles, kept {n_acc} via same-accession support)")
