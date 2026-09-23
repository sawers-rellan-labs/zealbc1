#!/bin/bash
#SBATCH --job-name=discovery_perdonor
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/discovery_perdonor_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/discovery_perdonor_%j.err
# Per-donor chr10 discovery = steps 1-4 of the Zd.0040_P1 v4 chain (pilot_1B_chr10/crisp_perdonor_Zd0040_bc2s3pool/j1.sbatch), with the donor
# as a parameter and the MAPQ20-realigned BC2S3 CRAMs (results/bc2s3_realign/cram, job 910286) as the witness pool instead of the old BAMs.
#   [1] witness pool: merge the donor's BC2S3 lines (checks excluded), chr10, one RG
#   [2] CRISP: donor BC1 pools + 2 B73 pools (zero class for eps_site in step 4) + witness pool, lowcopy union BED, -p 12
#   [3] witness veto (veto_witness.py)   [4] step 4 (pilot_step4_postfilter_llr.py; witness mapped as its own 'donor')
# Env: DONOR (e.g. Zx.0100_P4)  [OUT=$Z/results/pilot_union_chr10/discovery/$DONOR] [LINES=file of line ids; default = the donor's rows
#      of meta/bc2s3_batch1_realign_lines.tsv] [CHECKS=PN10_SID893 (comma-sep, excluded)] [DEMUXQC=$Z/results/gate2/demux_qc/demux_qc.tsv]
# Run:  DONOR=Zx.0100_P4 sbatch --parsable --export=ALL $Z/code/PHG/bin/discovery_perdonor.sbatch
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONOR:?set DONOR}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; R=$Z/results; P=$R/pilot_1B_chr10; U=$P/union; B=$R/crisp_bench; BIN=$Z/code/PHG/bin
SAM=$Z/envs/phgv2-conda/bin/samtools; FA=$Z/reference/B73.fa
OUT=${OUT:-$R/pilot_union_chr10/discovery/$DONOR}; DEMUXQC=${DEMUXQC:-$R/gate2/demux_qc/demux_qc.tsv}; CHECKS=${CHECKS:-PN10_SID893}
TAG=${DONOR//./}; W=${TAG}_BC2S3_chr10
export TMPDIR=$Z/work/tmp_$SLURM_JOB_ID; mkdir -p $TMPDIR $OUT; cd $OUT
echo "[discovery] donor $DONOR | out $OUT | code $(git -C $Z/code log --oneline -1 2>/dev/null || echo unknown)"

echo "== [1] witness pool  $(date)"
if [ -n "$LINES" ]; then cp $LINES lines.txt; else awk -F'\t' -v d=$DONOR '$2==d{print $1}' $Z/code/meta/bc2s3_batch1_realign_lines.tsv > lines.txt; fi
: > line_crams.txt
for s in $(cat lines.txt); do
  case ",$CHECKS," in *",$s,"*) echo "  skip check $s"; continue;; esac
  f=$R/bc2s3_realign/cram/$s.cram; [ -s $f ] || { echo "missing $f"; exit 1; }; echo $f >> line_crams.txt
done
echo "lines in witness pool: $(wc -l < line_crams.txt)"
$SAM merge -@8 --reference $FA -R chr10 -o $TMPDIR/merged.bam -b line_crams.txt
$SAM addreplacerg -@8 -r "ID:$W" -r "SM:$W" -r "PL:ILLUMINA" -m overwrite_all -o $OUT/$W.bam $TMPDIR/merged.bam && rm -f $TMPDIR/merged.bam
$SAM index -@8 $OUT/$W.bam; $SAM flagstat -@8 $OUT/$W.bam | head -1

echo "== [2] CRISP  $(date)"
awk -F'\t' -v d=$DONOR 'NR>1 && $4==d{print $3}' $DEMUXQC > bc1_samples.txt
[ -s bc1_samples.txt ] || { echo "no BC1 samples for $DONOR in $DEMUXQC"; exit 1; }
: > bams.txt
for s in $(cat bc1_samples.txt); do
  f=""; for d in $R/cram $R/align_membench; do [ -s $d/$s.cram ] && { f=$d/$s.cram; break; }; done
  [ -n "$f" ] || { echo "no CRAM for $s"; exit 1; }; echo $f >> bams.txt
done
printf "%s\n" $R/b73_control/ERR3288215/B73_ERR3288215.cram $R/b73_control/skim10/B73_skim10.bam $OUT/$W.bam >> bams.txt; cat bams.txt
set +e; /usr/bin/time -v $Z/envs/crisp/bin/CRISP.binary --bams bams.txt --ref $FA --bed $U/union_chr10.bed --regions chr10 -p 12 --sm 0 --mmq 20 --filterreads 0 --minc 2 --VCF crisp_all.vcf > crisp.log 2> crisp.time; set -e   # CRISP exits 1 on success
grep -E "Elapsed|Maximum resident" crisp.time; echo "records: $(grep -vc '^#' crisp_all.vcf)"; bgzip -f crisp_all.vcf; tabix -f -p vcf crisp_all.vcf.gz

echo "== [3] witness veto  $(date)"
python3 $BIN/veto_witness.py crisp_all.vcf.gz crisp_filtered.vcf $W; bgzip -f crisp_filtered.vcf; tabix -f -p vcf crisp_filtered.vcf.gz

echo "== [4] step 4  $(date)"
{ head -1 $DEMUXQC; awk -F'\t' -v d=$DONOR 'NR>1 && $4==d' $DEMUXQC
  printf "B73\t0\tB73_ERR3288215\tB73\tmaize\t0\t15.5\t0\t0\t0\t0\t0\tFALSE\nB73\t0\tB73_skim10\tB73\tmaize\t0\t5.7\t0\t0\t0\t0\t0\tFALSE\n"
  printf "BC2S3\t0\t$W\t${TAG}_BC2S3\tmapping_pop\t0\t15\t0\t0\t0\t0\t0\tFALSE\n"; } > map.tsv
python3 $BIN/pilot_step4_postfilter_llr.py --vcf crisp_filtered.vcf.gz --map map.tsv --til18 $B/TIL18_vs_B73_chr10_snps.tsv --gigi $B/Gigi_vs_B73_chr10_snps.tsv \
  --schnable $P/panel_positions/wideseq_chr10_snps.tsv --mgdb $P/panel_positions/mgdb26hq_chr10_snps.tsv --out step4 2>&1 | tail -2
cat step4/summary.tsv
zcat step4/$DONOR.sites.tsv.gz | awk -F'\t' 'NR>1{t[$12]++} END{for(k in t) print "tier", k, t[k]}' | sort
rm -rf $TMPDIR; echo "DONE $(date)"
