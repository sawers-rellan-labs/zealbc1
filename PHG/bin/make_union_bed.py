import sys
ranges_f,nonrep_f,out=sys.argv[1:4]
iv=[(int(x[1]),int(x[2])) for x in (l.split() for l in open(ranges_f)) if x[0]=="chr10" and not x[3].startswith("intergenic")]
iv+=[(int(x[1]),int(x[2])) for x in (l.split() for l in open(nonrep_f)) if x[0]=="chr10"]
iv.sort(); m=[]
for s,e in iv:
    if m and s-m[-1][1]<=200: m[-1][1]=max(m[-1][1],e)
    else: m.append([s,e])
m=[(s,e) for s,e in m if e-s>=500]
with open(out,"w") as o:
    for s,e in m: o.write(f"chr10\t{s}\t{e}\tu_chr10:{s}-{e}\t0\t.\n")
L=sorted(e-s for s,e in m); n=len(L)
print(f"union BED: {n} ranges, {sum(L)/1e6:.1f} Mb, median {L[n//2]} bp, deciles {[L[int(n*q/10)] for q in range(1,10,2)]}")
