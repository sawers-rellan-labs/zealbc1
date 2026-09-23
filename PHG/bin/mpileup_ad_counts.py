#!/usr/bin/env python3
"""Per-pool (n, a) at CRISP's alleles from `bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'` over an mpileup of extra pools
(the B73 controls in donor_discovery_chr10.sbatch), for pilot_step4_postfilter_llr.py --extra-counts.
Usage: mpileup_ad_counts.py query.tsv sites.tsv bams.txt out.tsv
  sites.tsv = chrom pos ref alt (the CRISP records); sample name = bam basename without extension, in bams.txt order (= query columns).
  a = AD of CRISP's alt allele (0 if mpileup did not see it); n = ref + a (other alleles ignored, as step 4 does)."""
import sys, os
q, sites, bams, out = sys.argv[1:5]
names = [os.path.basename(l.strip()).rsplit('.', 1)[0] for l in open(bams) if l.strip()]
want = {}
for l in open(sites):
    c, p, r, a = l.split()[:4]; want.setdefault((c, int(p)), []).append((r, a))
nrec = 0
with open(q) as f, open(out, 'w') as o:
    o.write('chrom\tpos\tref\talt\tsample\tn\ta\n')
    for l in f:
        x = l.rstrip('\n').split('\t'); c, p, r = x[0], int(x[1]), x[2]
        alts = [] if x[3] == '.' else x[3].split(','); ads = x[4:]
        if len(ads) != len(names): sys.exit(f"{len(ads)} AD columns vs {len(names)} bams")
        for cr, ca in want.get((c, p), []):
            if cr != r: continue
            ai = alts.index(ca) + 1 if ca in alts else None
            for smp, ad in zip(names, ads):
                v = [int(t) if t not in ('.', '') else 0 for t in ad.split(',')] if ad not in ('.', '') else [0]
                k = v[ai] if ai is not None and ai < len(v) else 0
                o.write(f"{c}\t{p}\t{cr}\t{ca}\t{smp}\t{v[0] + k}\t{k}\n")
            nrec += 1
print(f"[mpileup_ad_counts] {nrec} records x {len(names)} pools -> {out}")
