# Plan (DRAFT) — pipeline v2: one Nextflow pipeline from raw libraries to marker layers

Status: **proposal, 2026-09-24**, for discussion; nothing implemented. Storage section (§5) pending the disk audit
(`nilhmm/bin/audit_du.sbatch`, job 946049, results in `ZEAL/results/audit_du_20260924/`).

## 1. Why
The chr10 pilots ran as standalone sbatch scripts next to the `nilhmm/` pipeline. The same step was implemented more than once and the
copies drifted: B73 pools inside vs outside CRISP, `mpileup -I` present in one pileup and missing in another, a demux QC table overwritten
by every pool run, no duplicate removal anywhere. Each fix had to be applied in several places. One module per step removes that class
of error. The known issues this plan must settle are in §4.

## 2. Principles
1. **One module per step**, used by every entry; no standalone copies of a module's command.
2. **Separate entries per stage** (the existing `--entry` dispatcher): editing a downstream module never puts upstream tasks at risk.
3. **Costly, reusable outputs live in a permanent store (`storeDir`), not in `work/`**: CRAMs, per-pool demux QC. A task whose stored
   output exists is skipped whatever changed in the module; rerunning it is a deliberate act (delete the stored output).
4. **Hash hygiene** (why cosmetic edits reran costly steps):
   - comments and notes outside the `script:` block (Groovy `//`), never as bash `#` inside it;
   - threads and memory read from Slurm at run time (`-@ \$SLURM_CPUS_PER_TASK`), not `${task.cpus}` / `${task.memory}` in the script,
     so reallocating resources does not change the task hash;
   - `cache 'lenient'` (path + size, not mtime) on processes with large inputs.
5. **Temporaries die inside the task** (merged witness pool, sort temp, demux FASTQs if merged with ALIGN; see §5).
6. Tools: minibwa, samtools, bcftools, CRISP, nilHMM, PHG — the ones in use; environments prebuilt and referenced by prefix.

## 3. Stages (entries) and modules
| # | entry | modules | per | main output (store) |
|---|---|---|---|---|
| 1 | `demux` | DEMUX (cutadapt exact inline, `-e 0 --no-indels`), DEMUX_QC | pool | per-sample FASTQ (transient), `demux_qc/<pool>.tsv` (store, one file per pool) |
| 2 | `align` | ALIGN (minibwa -x sr) → **MARKDUP** → CRAM (MAPQ 20, `-F 0x904`, duplicates flagged or removed) → MOSDEPTH | sample (BC1 sample, BC2S3 line, B73 pool) | `cram/<sample>.cram` (store) |
| 3 | `discovery` | WITNESS_POOL → CRISP (BC1 samples + witness only) → VETO → B73_COUNTS (`mpileup -I`) → STEP4 | donor × chr | `step4/<donor>.sites.tsv.gz` |
| 4 | `union` | UNION (tier-A sites of the donor set; multi-allelic dropped) | donor set × chr | `union/<set>_<chr>.tsv.gz` |
| 5 | `count_once` | COUNT_SAMPLE (`mpileup -I -T union`, one task per sample) → JOINT_STEP4 → GAP_FILL (`dhd_bayes`) | sample / donor set × chr | donor allele table |
| 6 | `layer1` | LINE_COUNTS → RTIGER (design BC2S3, rigidity 500) | donor × chr | ancestry segments per line |
| 7 | `layer2` | FOUNDER (gVCF → pseudo-assembly) → PHG_DB → PHG_IMPUTE (pairwise: B73 + donor, that donor's lines; F = 0, stay 0.99999) → RASTERIZE | donor × chr | genotypes at the union sites |
| 8 | `report` | PAINT, summary tables, KS / single-locus checks | donor × chr | paintings, tables |

Donor sets for stages 4–5 are named in a run card (`docs/runs/<run>.md`: purpose, donors with BC1 count / lines / coverage, exclusions),
and each entry checks the run card before starting.

## 4. Known issues and where each is settled
| # | issue (found 2026-09-20 → 24) | settled in | decision needed |
|---|---|---|---|
| 1 | No duplicate removal (BC1, lines, B73 pools). Pooled-caller benchmark and GATK best practices remove/mark PCR duplicates [1, 2]; CRISP paper silent [3] | stage 2 MARKDUP | tool (samtools markdup [4] needs collate/fixmate; Picard MarkDuplicates works on coordinate-sorted); mark vs remove; rerun existing CRAMs? |
| 2 | CRISP run with `--filterreads 0` (its mismatch filter off; the CRISP paper used ≤ 3 mismatches, MAPQ ≥ 20, base quality ≥ 17 [3]) | stage 3 CRISP | turn it back on? |
| 3 | Insertion records at a SNP position overwrite its counts (ALT → 0) unless `mpileup -I` | one COUNTS helper used by stages 3, 5, 6 | make the helper skip indel records itself |
| 4 | Demux QC table overwritten by every pool run | stage 1 DEMUX_QC | one file per pool (store) |
| 5 | Witness veto depends on witness depth (10 lines at 0.4x kept 20% of records) | stage 3 VETO | keep "≥ 1 ALT read" or make it depth-aware |
| 6 | Tiers depend on the count source (CRISP vs mpileup disagreed at ~15% of own tier-A sites) | stages 3 vs 5 | which counts define tiers |
| 7 | RTIGER rigidity fixed at 500 vs a density-scaled rule | stage 6 | confirm 500 |
| 8 | Union / count-once / gap filling / layer 1 exist only as standalone scripts (`PHG/bin/`) | stages 4–6 | port as modules |

## 5. Storage, caching and cleanup — PENDING the disk audit
To be written from the audit: size and file count per write location (`ZEAL/{code,work,reference,envs}`, every `results/<subdir>`, every
Nextflow `work/`, the PHG databases, `/share/maize/frodrig4/{conda,tmp}`). Questions it must answer:
- `workDir` on `/share/maize/frodrig4/...` (2 TB, not persistent) instead of `${params.outdir}/work` on `/rsstu` (today's setting);
- demux FASTQs (~150 GB per pool): capped concurrent DEMUX + deletion after ALIGN, or one DEMUX+ALIGN task per pool writing FASTQs to
  its own scratch;
- what moves to `storeDir` (CRAMs, demux QC, step-4 tables?) and where;
- which existing `results/` subdirectories are kept, archived or deleted (pilot version soup, superseded runs, stub outputs);
- routine: `nextflow clean -f -but <last run>` after each successful run, plus a size report.

## 6. Supervision
Per run: own launch dir and `workDir`; a post-run check (work size, failed tasks, published outputs); monitors, not sleep loops; long runs
watched with the session kept open (`/loop`), acting only as the run card allows.

## 7. Open decisions (summary)
1–7 of §4; storage rules of §5; which existing CRAMs and tables are reused vs regenerated after MARKDUP.

## References
1. Huang HW, Mullikin JC, Hansen NF. Evaluation of variant detection software for pooled next-generation sequence data.
   *BMC Bioinformatics* 2015;16:235. doi:10.1186/s12859-015-0624-y
2. Van der Auwera GA, Carneiro MO, Hartl C, Poplin R, et al. From FastQ data to high-confidence variant calls: the Genome Analysis
   Toolkit best practices pipeline. *Curr Protoc Bioinformatics* 2013;43:11.10.1–11.10.33. doi:10.1002/0471250953.bi1110s43
3. Bansal V. A statistical method for the detection of variants from next-generation resequencing of DNA pools.
   *Bioinformatics* 2010;26(12):i318–i324. PMC2881398
4. Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience* 2021;10(2):giab008. (not re-checked this session)
