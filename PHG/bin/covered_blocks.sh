#!/bin/bash
# Pilot step 5 helper (2026-09-17): depth-based "covered" intervals per donor over a BED (reference blocks only where the donor's
# pooled reads are present). Usage: covered_blocks.sh <donor> <bed> <out.bed> <cram...>   [MINDP=8] [SAM=samtools]
# Summed depth over the donor's pools (MAPQ>=20) >= MINDP -> merged BED of covered stretches.
set -eo pipefail
donor=$1; bed=$2; out=$3; shift 3; MINDP=${MINDP:-8}; SAM=${SAM:-samtools}
$SAM depth -a -b "$bed" -Q 20 "$@" | awk -v m=$MINDP 'BEGIN{OFS="\t"} {s=0; for(i=3;i<=NF;i++) s+=$i; if(s>=m){ if(c==$1 && $2==e+1){e=$2} else { if(c!="") print c, st-1, e; c=$1; st=$2; e=$2 } } } END{if(c!="") print c, st-1, e}' > "$out"
echo "$donor: covered stretches $(wc -l < $out), $(awk '{s+=$3-$2} END{printf "%.2f", s/1e6}' $out) Mb at summed depth >= $MINDP"
