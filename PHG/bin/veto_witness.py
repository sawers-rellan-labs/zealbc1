#!/usr/bin/env python3
"""Witness veto (pilot filter_bc2s3pool.py, 2026-09-21): keep CRISP records whose witness pool (the donor's merged BC2S3 lines)
shows >= 1 ALT read (ADf+ADr+ADb alt counts); drop the rest. Usage: veto_witness.py in.vcf.gz out.vcf <witness_sample_name>"""
import sys, gzip
vin, vout, pool = sys.argv[1:4]; kept = dropped = 0
with gzip.open(vin, 'rt') as f, open(vout, 'w') as o:
    for l in f:
        if l.startswith('##'): o.write(l); continue
        x = l.rstrip('\n').split('\t')
        if l.startswith('#'):
            o.write(l)
            if pool not in x: sys.exit(f"witness sample {pool} not in VCF columns: {x[9:]}")
            idx = x.index(pool); continue
        fmt = x[8].split(':'); vals = x[idx].split(':'); alt = 0
        for k in ('ADf', 'ADr', 'ADb'):
            if k in fmt:
                v = vals[fmt.index(k)].split(',')
                if len(v) >= 2 and v[1] not in ('.', ''): alt += int(v[1])
        if alt >= 1: o.write(l); kept += 1
        else: dropped += 1
print(f"records kept (witness ALT>=1): {kept} | dropped (no ALT read in the mapping population): {dropped}")
