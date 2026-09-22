#!/usr/bin/env python3
"""bcftools AD table -> RTIGER poolseq counts (pilot inline script, 2026-09-21).
Usage: ad_to_counts.py alleles.tsv.gz ad.tsv counts.tsv   (alleles: chrom pos ref,alt; ad.tsv from `bcftools query -H -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'`)
Sample names: bcftools' header tokens `[n]SAMPLE:AD` -> SAMPLE (path-named RG-less BAMs are reduced to their basename)."""
import re, gzip, sys
alleles_f, ad_f, out_f = sys.argv[1:4]
fa = {}
for l in gzip.open(alleles_f, 'rt'):
    c, p, a = l.rstrip('\n').split('\t'); r, alt = a.split(','); fa[int(p)] = (r, alt)
samples = None; n = 0
with open(out_f, 'w') as out:
    out.write("SAMPLE\tCONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n")
    for l in open(ad_f):
        x = l.rstrip('\n').split('\t')
        if l[0] == '#':
            samples = [re.sub(r'^\[\d+\]', '', c.rsplit(':', 1)[0]).split('/')[-1].replace('_sorted_alignment.bam', '').replace('.cram', '') for c in x[4:]]; continue
        pos = int(x[1]); chrom = x[0]
        if pos not in fa: continue
        r, alt = fa[pos]; alts = [] if x[3] == '.' else x[3].split(','); ai = alts.index(alt) + 1 if alt in alts else None
        for smp, ad in zip(samples, x[4:]):
            if ad in ('.', ''): continue
            c = [int(v) if v != '.' else 0 for v in ad.split(',')]; rc = c[0]; ac = c[ai] if (ai is not None and ai < len(c)) else 0
            if rc + ac > 0: out.write(f"{smp}\t{chrom}\t{pos}\t{rc}\t{ac}\t{r}\t{alt}\n"); n += 1
print(f"counts: {n} sample x site observations with reads, {len(samples or [])} samples")
