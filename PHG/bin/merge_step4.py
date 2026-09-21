import sys,gzip
a,b,out=sys.argv[1:4]   # a = gene-space table (preferred), b = non-repeat table
seen=set(); n_a=n_b=0
with gzip.open(out,"wt") as o:
    with gzip.open(a,"rt") as f:
        hdr=f.readline(); o.write(hdr)
        for l in f: seen.add(int(l.split("\t",2)[1])); o.write(l); n_a+=1
    with gzip.open(b,"rt") as f:
        f.readline()
        for l in f:
            if int(l.split("\t",2)[1]) not in seen: o.write(l); n_b+=1
print(f"merged: {n_a} from gene-space + {n_b} new from non-repeat = {n_a+n_b}")
