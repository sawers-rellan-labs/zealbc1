#!/usr/bin/env python3
"""Per-pool (n, a) counts (mpileup_ad_counts.py output: chrom pos ref alt sample n a) -> a minimal VCF in the CRISP layout that
pilot_step4_postfilter_llr.py reads (FORMAT ADf:ADr:ADb; the counts go in ADf, ADr/ADb = 0,0). Lets the second counting pass of the
marker union pilot be classified by the unchanged step-4 code. Samples missing at a record get 0,0.
Usage: counts_to_crisp_vcf.py counts.tsv out.vcf"""
import sys, collections
cin, vout = sys.argv[1:3]
rec = collections.OrderedDict(); samples = []
with open(cin) as f:
    f.readline()
    for l in f:
        c, p, r, a, s, n, k = l.rstrip('\n').split('\t')
        if s not in samples: samples.append(s)
        rec.setdefault((c, int(p), r, a), {})[s] = (int(n) - int(k), int(k))
with open(vout, 'w') as o:
    o.write('##fileformat=VCFv4.2\n##source=counts_to_crisp_vcf (second counting pass; counts in ADf)\n')
    o.write('##FORMAT=<ID=ADf,Number=2,Type=Integer,Description="ref,alt reads (all strands)">\n')
    o.write('##FORMAT=<ID=ADr,Number=2,Type=Integer,Description="unused (0,0)">\n##FORMAT=<ID=ADb,Number=2,Type=Integer,Description="unused (0,0)">\n')
    o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
    for (c, p, r, a), v in sorted(rec.items(), key=lambda z: (z[0][0], z[0][1])):
        cells = [f"{v.get(s, (0, 0))[0]},{v.get(s, (0, 0))[1]}:0,0:0,0" for s in samples]
        o.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tADf:ADr:ADb\t" + '\t'.join(cells) + '\n')
print(f"[counts_to_crisp_vcf] {len(rec)} records x {len(samples)} samples -> {vout}")
