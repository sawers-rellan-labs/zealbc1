import csv,sys,collections,statistics as st
U=sys.argv[1]; OUT=sys.argv[2]; CHRLEN=152435371
counts=collections.defaultdict(lambda: collections.defaultdict(int))   # counts[line][pos] = ref+alt
sites=set()
for r in csv.DictReader(open(f"{U}/rtiger_v5/counts.tsv"),delimiter='\t'):
    p=int(r['POSITION']); sites.add(p); counts[r['SAMPLE']][p]+=int(r['REF_COUNT'])+int(r['ALT_COUNT'])
sites=sorted(sites); lines=sorted(counts)
tot={l:sum(counts[l].values()) for l in lines}
print(f"lines {len(lines)} | tier-A sites with counts {len(sites)} | reads at sites per line: median {st.median(tot.values()):.0f} min {min(tot.values())} max {max(tot.values())}")
union=[(int(x[1]),int(x[2])) for x in (l.split('\t') for l in open(f"{U}/union_chr10.bed"))]
def bins(w): return [(s,min(s+w,CHRLEN)) for s in range(0,CHRLEN,w)]
schemes={"union_ranges":union,"bin_250kb":bins(250_000),"bin_1Mb":bins(1_000_000)}
import bisect
w=csv.writer(open(f"{OUT}/reads_per_bin_summary.tsv","w"),delimiter='\t')
w.writerow(["scheme","n_bins","bins_with_site","sites_per_bin_median","reads_per_bin_mean","d10","d20","d30","d40","d50","d60","d70","d80","d90","frac_bins_zero_reads"])
for name,B in schemes.items():
    nsite=[]; vals=[]
    for a,b in B:
        i=bisect.bisect_left(sites,a+1); j=bisect.bisect_right(sites,b); ss=sites[i:j]; nsite.append(len(ss))
        for l in lines: vals.append(sum(counts[l].get(p,0) for p in ss))
    vals.sort(); n=len(vals); dec=[vals[int(n*q/10)] for q in range(1,10)]
    row=[name,len(B),sum(1 for k in nsite if k>0),st.median(nsite),round(st.mean(vals),2)]+dec+[round(sum(1 for v in vals if v==0)/n,3)]
    w.writerow(row); print("\t".join(map(str,row)))
    # per-line mean for the three coverage tertiles
    per_line={l:st.mean(sum(counts[l].get(p,0) for p in sites[bisect.bisect_left(sites,a+1):bisect.bisect_right(sites,b)]) for a,b in B) for l in lines}
    o=sorted(per_line.items(),key=lambda x:x[1]); k=len(o)
    print(f"   per-line mean reads/bin: lowest {o[0][0]} {o[0][1]:.2f} | median {o[k//2][0]} {o[k//2][1]:.2f} | highest {o[-1][0]} {o[-1][1]:.2f}")
