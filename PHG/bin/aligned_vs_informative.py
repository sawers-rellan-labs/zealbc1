import sys,csv,re,glob,os,bisect,subprocess,statistics as st,collections
U,G,BAMDIR,OUT=sys.argv[1:5]; CHRLEN=152435371
ranges=[(int(x[1]),int(x[2])) for x in (l.rstrip('\n').split('\t') for l in open(f"{U}/union_chr10.bed"))]
starts=[a for a,b in ranges]
# hapid -> (range index, sample) from the two founder hVCFs
hap={}
for fn in glob.glob(f"{G}/hvcf/*.h.vcf"):
    for l in open(fn):
        if l.startswith('##ALT'):
            m=re.search(r'ID=([0-9a-f]+).*SampleName=([^,>]+).*Regions=chr10:(\d+)-(\d+)',l)
            if m:
                i=bisect.bisect_right(starts,int(m.group(3))-1)-1   # hvcf regions are 1-based; bed 0-based
                hap[m.group(1)]=(i,m.group(2))
samples=sorted(os.path.basename(p).split('_chr10')[0] for p in glob.glob(f"{G}/readmap/*_readMapping.txt"))
print(f"ranges {len(ranges)} | hapids {len(hap)} | samples {len(samples)}")
R=len(ranges); inf={s:[0]*R for s in samples}; amb={s:[0]*R for s in samples}; multi={s:0 for s in samples}
for s in samples:
    for l in open(f"{G}/readmap/{s}_chr10_R1_readMapping.txt"):
        if l[0]=='#' or l.startswith('HapIds'): continue
        ids,k=l.rstrip('\n').split('\t'); k=int(k); h=[hap[i] for i in ids.split(',') if i in hap]
        if not h: continue
        rs={x[0] for x in h}; ss={x[1] for x in h}
        if len(rs)>1: multi[s]+=k
        elif len(ss)==1: inf[s][rs.pop()]+=k
        else: amb[s][rs.pop()]+=k
# aligned reads per range from the BAMs (MAPQ>=20, chr10 union ranges)
ali={}
for s in samples:
    b=glob.glob(f"{BAMDIR}/filtered_S*/{s}_sorted_alignment.bam")
    if not b: print("no BAM for",s); continue
    out=subprocess.run(["samtools","bedcov","-c","-Q","20",f"{U}/union_chr10.bed",b[0]],capture_output=True,text=True,check=True).stdout
    ali[s]=[int(l.split('\t')[-1]) for l in out.rstrip('\n').split('\n')]
samples=[s for s in samples if s in ali]
w=csv.writer(open(f"{OUT}/aligned_vs_informative_per_range.tsv","w"),delimiter='\t'); w.writerow(["sample","range_idx","start","end","aligned","informative","ambiguous"])
for s in samples:
    for i,(a,b) in enumerate(ranges): w.writerow([s,i,a,b,ali[s][i],inf[s][i],amb[s][i]])
def dec(v): v=sorted(v); n=len(v); return [v[int(n*q/10)] for q in range(1,10)]
def summarize(name,binner,B):
    A=[];I=[];M=[]
    for s in samples:
        a=[0]*len(B); ii=[0]*len(B); m=[0]*len(B)
        for r,(x,y) in enumerate(ranges):
            j=binner(x); a[j]+=ali[s][r]; ii[j]+=inf[s][r]; m[j]+=amb[s][r]
        A+=a; I+=ii; M+=m
    ratio=sum(I)/max(1,sum(A))
    print(f"{name}\tbins {len(B)}\taligned mean {st.mean(A):.1f} dec {dec(A)}\tinformative mean {st.mean(I):.2f} dec {dec(I)}\tambiguous mean {st.mean(M):.2f}\tinformative/aligned {ratio:.3f}\tcells informative=0 {sum(1 for v in I if v==0)/len(I):.3f}")
summarize("union_ranges",lambda x: bisect.bisect_right(starts,x)-1,ranges)
for w_,name in [(250_000,"bin_250kb"),(1_000_000,"bin_1Mb")]:
    B=list(range(0,CHRLEN,w_)); summarize(name,lambda x,w_=w_: x//w_,B)
print("multi-range read pairs per sample (median):",st.median(multi.values()))
