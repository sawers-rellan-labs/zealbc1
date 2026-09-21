#!/usr/bin/env python3
"""Pilot step 4 — post-filter + per-donor pooled likelihood on CRISP record-level output (BC1 6-plant pools, 12 haplotypes).

Usage: pilot_step4_postfilter_llr.py --vcf crisp.vcf[.gz] --map demux_qc.tsv --out OUTDIR
         [--til18 TIL18.tsv --gigi Gigi.tsv --schnable schn.tsv --mgdb mgdb.tsv] [--eps0 0.005] [--prior 0.5]
Per record (biallelic SNPs only): per pool n_i (ref+alt reads, ADf+ADr+ADb) and a_i (alt reads).
Per donor d (pools i in d):  H1: j_i~Binom(6,1/2) carriers, f_j=j/12, p_j=f_j(1-eps)+(1-f_j)eps, L1_i=sum_j w_j Binom(a_i|n_i,p_j)
                              H0: L0_i=Binom(a_i|n_i,eps);  LLR_d = sum_i log L1_i - log L0_i
eps_site: from pools of confidently non-carrier donors (LLR < -3 at eps0) when they hold >= 20 reads, else eps0 (flag no_zero_class).
Flags: hidepth (summed depth > 2x median), af_gt_half (donor alt fraction > 0.5, binomial p<0.01), inconsistent (multi-sample donor:
       one pool a>=3 while another pool has n>=10 and a=0), single_sample.
Tiers (posterior with prior pi): A: LLR>=6.9 & >=2 pools with alt & no flags | B: LLR>=4.6 & no hidepth/af flags | C: LLR>=2.2 | ref: LLR<=-4 & n>=12.
Outputs: OUTDIR/<donor>.sites.tsv.gz (all records with donor depth>0), OUTDIR/summary.tsv, OUTDIR/pool_qc.tsv, OUTDIR/run_info.txt
"""
import sys, gzip, math, argparse, collections, statistics
ap = argparse.ArgumentParser()
ap.add_argument('--vcf', required=True); ap.add_argument('--map', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--til18'); ap.add_argument('--gigi'); ap.add_argument('--schnable'); ap.add_argument('--mgdb')
ap.add_argument('--eps0', type=float, default=0.005); ap.add_argument('--prior', type=float, default=0.5)
ap.add_argument('--plants', type=int, default=6)
A = ap.parse_args()
import os; os.makedirs(A.out, exist_ok=True)
PL = A.plants; H = 2 * PL
LOGW = [math.lgamma(PL+1) - math.lgamma(j+1) - math.lgamma(PL-j+1) + PL*math.log(0.5) for j in range(PL+1)]
def logbinom(n, k, p):
    p = min(max(p, 1e-9), 1-1e-9)
    return math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1) + k*math.log(p) + (n-k)*math.log(1-p)
def logsumexp(v):
    m = max(v); return m + math.log(sum(math.exp(x-m) for x in v))
def llr_pool(n, a, eps):
    if n == 0: return 0.0
    l1 = logsumexp([LOGW[j] + logbinom(n, a, (j/H)*(1-eps) + (1-j/H)*eps) for j in range(PL+1)])
    return l1 - logbinom(n, a, eps)
def binom_sf_half(n, a):  # P(X >= a | n, 0.5)
    return sum(math.exp(logbinom(n, k, 0.5)) for k in range(a, n+1))
# ---- sample -> donor
smap = {}; taxon = {}
with open(A.map) as f:
    hdr = f.readline().rstrip('\n').split('\t'); iS, iD, iT = hdr.index('Sample_Id'), hdr.index('donor'), hdr.index('taxon')
    for l in f:
        x = l.rstrip('\n').split('\t'); smap[x[iS]] = x[iD]; taxon[x[iD]] = x[iT]
def loadset(p):
    S = set()
    if p:
        for l in open(p):
            c, pos, r, a = l.split()[:4]; S.add((c, int(pos), r, a))
    return S
TIL, GIG, SCH, MG = loadset(A.til18), loadset(A.gigi), loadset(A.schnable), loadset(A.mgdb)
# ---- parse VCF
op = gzip.open if A.vcf.endswith('.gz') else open
recs = []  # (chrom,pos,ref,alt, {pool:(n,a)})
pools = []
with op(A.vcf, 'rt') as f:
    for l in f:
        if l.startswith('##'): continue
        if l.startswith('#CHROM'): pools = l.rstrip('\n').split('\t')[9:]; continue
        x = l.rstrip('\n').split('\t')
        if len(x[3]) != 1 or len(x[4]) != 1 or ',' in x[4]: continue
        fmt = x[8].split(':'); 
        try: iF, iR, iB = fmt.index('ADf'), fmt.index('ADr'), fmt.index('ADb')
        except ValueError: continue
        cnt = {}
        for p, s in zip(pools, x[9:]):
            v = s.split(':'); r = a = 0
            for k in (iF, iR, iB):
                if k < len(v) and ',' in v[k]:
                    q = v[k].split(','); 
                    try: r += int(q[0]); a += int(q[1])
                    except ValueError: pass
            cnt[p] = (r + a, a)
        recs.append((x[0], int(x[1]), x[3], x[4], cnt))
donors = sorted(set(smap[p] for p in pools if p in smap)); dpools = {d: [p for p in pools if smap.get(p) == d] for d in donors}
tot = [sum(n for n, a in r[4].values()) for r in recs]; med = statistics.median(tot) if tot else 0
info = [f"records(biallelic SNPs)={len(recs)} pools={len(pools)} donors={len(donors)} median_total_depth={med} eps0={A.eps0} prior={A.prior}"]
# ---- per record, per donor
out = {d: gzip.open(f"{A.out}/{d}.sites.tsv.gz", 'wt') for d in donors}
cols = "chrom\tpos\tref\talt\tn\ta\tn_pools\tn_pools_alt\teps\tLLR\tposterior\ttier\tflags\tin_TIL18\tin_Gigi\tin_schnable\tin_mgdb26\tpool_counts\n"
for d in donors: out[d].write(cols)
summ = collections.defaultdict(collections.Counter); poolqc = collections.defaultdict(lambda: [0, 0, 0, 0])  # pool: n_sites_depth, alt_reads, ref_reads, af>0.5 sites
logit_prior = math.log(A.prior/(1-A.prior))
for chrom, pos, ref, alt, cnt in recs:
    T = sum(n for n, a in cnt.values()); hidepth = T > 2*med
    llr0 = {d: sum(llr_pool(*cnt[p], A.eps0) for p in dpools[d]) for d in donors}
    zn = sum(cnt[p][0] for d in donors if llr0[d] < -3 for p in dpools[d]); za = sum(cnt[p][1] for d in donors if llr0[d] < -3 for p in dpools[d])
    if zn >= 20: eps = max((za + 0.5)/(zn + 1), 0.001); nz = ''
    else: eps = A.eps0; nz = 'no_zero_class'
    key = (chrom, pos, ref, alt); prov = (key in TIL, key in GIG, key in SCH, key in MG)
    for d in donors:
        ps = dpools[d]; n = sum(cnt[p][0] for p in ps); a = sum(cnt[p][1] for p in ps)
        if n == 0: continue
        npa = sum(1 for p in ps if cnt[p][1] > 0); LLR = sum(llr_pool(*cnt[p], eps) for p in ps)
        post = 1/(1+math.exp(-(LLR + logit_prior))) if LLR > -700 else 0.0
        flags = []
        if hidepth: flags.append('hidepth')
        if n >= 4 and a/n > 0.5 and binom_sf_half(n, a) < 0.01: flags.append('af_gt_half')
        if len(ps) == 1: flags.append('single_sample')
        elif any(cnt[p][1] >= 3 for p in ps) and any(cnt[p][0] >= 10 and cnt[p][1] == 0 for p in ps): flags.append('inconsistent')
        if nz: flags.append(nz)
        hard = {'hidepth', 'af_gt_half'}
        if LLR >= 6.9 and npa >= 2 and not (set(flags) & (hard | {'inconsistent'})): tier = 'A'
        elif LLR >= 4.6 and not (set(flags) & hard): tier = 'B'
        elif LLR >= 2.2: tier = 'C'
        elif LLR <= -4 and n >= 12: tier = 'ref'
        else: tier = '-'
        summ[d][tier] += 1; summ[d]['tested'] += 1
        if tier in 'ABC':
            tx = taxon[d][:2] if d else ''
            summ[d][f'{tier}_inTIL18'] += prov[0]; summ[d][f'{tier}_inGigi'] += prov[1]; summ[d][f'{tier}_inSchnable'] += prov[2]; summ[d][f'{tier}_inMgdb'] += prov[3]
        for p in ps:
            q = poolqc[p]; q[0] += cnt[p][0] > 0; q[1] += cnt[p][1]; q[2] += cnt[p][0] - cnt[p][1]
            if cnt[p][0] >= 4 and cnt[p][1]/cnt[p][0] > 0.5 and binom_sf_half(cnt[p][0], cnt[p][1]) < 0.01: q[3] += 1
        out[d].write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{n}\t{a}\t{len(ps)}\t{npa}\t{eps:.4f}\t{LLR:.2f}\t{post:.4f}\t{tier}\t{','.join(flags) or '.'}\t{int(prov[0])}\t{int(prov[1])}\t{int(prov[2])}\t{int(prov[3])}\t{';'.join(f'{p}:{cnt[p][1]}/{cnt[p][0]}' for p in ps)}\n")
for d in donors: out[d].close()
with open(f"{A.out}/summary.tsv", 'w') as f:
    f.write("donor\ttaxon\tn_pools\ttested\tA\tB\tC\tref\tA_inAsm\tB_inAsm\tC_inAsm\tA_inSchnable\tA_inMgdb\tasm_used\n")
    for d in donors:
        s = summ[d]; asm = 'TIL18' if taxon[d] == 'mexicana' else ('Gigi' if taxon[d] == 'diploperennis' else 'none')
        ia = lambda t: (s[f'{t}_inTIL18'] if asm == 'TIL18' else s[f'{t}_inGigi'] if asm == 'Gigi' else 0)
        f.write(f"{d}\t{taxon[d]}\t{len(dpools[d])}\t{s['tested']}\t{s['A']}\t{s['B']}\t{s['C']}\t{s['ref']}\t{ia('A')}\t{ia('B')}\t{ia('C')}\t{s['A_inSchnable']}\t{s['A_inMgdb']}\t{asm}\n")
with open(f"{A.out}/pool_qc.tsv", 'w') as f:
    f.write("pool\tdonor\tsites_with_depth\talt_reads\tref_reads\talt_frac\tsites_af_gt_half\n")
    for p in pools:
        q = poolqc[p]; f.write(f"{p}\t{smap.get(p,'?')}\t{q[0]}\t{q[1]}\t{q[2]}\t{q[1]/max(q[1]+q[2],1):.4f}\t{q[3]}\n")
open(f"{A.out}/run_info.txt", 'w').write('\n'.join(info) + '\n')
print('\n'.join(info)); print(open(f"{A.out}/summary.tsv").read()); print(open(f"{A.out}/pool_qc.tsv").read())
