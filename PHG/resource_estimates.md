# PHG — measured resources (chr10 pilot)

Measured on hazel (`seff` / job resource summary), phg-debug. Requests are tightened to peak + headroom
as each step runs. Full-genome extrapolation anchors live in `../agent/PHG_PILOT_chr10.md`.

## Inputs confirmed on hazel (2026-09-12)
- B73 v5 FASTA + `.fai` + gene GFF (`Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3`), chr10 = `chr10` (152.4 Mb).
- 5 taxon references in `/rsstu/.../BZea/ref/`, all `chr1..chr10` + `alt-scaf_*`:
  Zx-TIL18 (chr10 144.6 Mb), Zv-TIL11, Zd-Gigi (chr10 147.2 Mb), Zh-RIMHU001 (all plain `.fa`, no `.fai`);
  **Zl-RIL003** downloaded from MaizeGDB, md5-verified, as bgzipped `.fa.gz`+`.fai`+`.gzi` (chr10 213.1 Mb).

## Stage: PHG_setup — pilot (chr10) measured

| step | job | cpu | mem peak | wall | notes |
|---|---|---|---|---|---|
| `02_subset_chr10` (samtools faidx 2 taxa + B73 + GFF) | 811363 | 2 | **6.25 GB** | **30 s** | faidx of two 2.5 GB genomes + chr10 extract; I/O not a bottleneck at chr10 scale. Request was 8 GB (78% used). |
| `03_align_taxa_to_b73` — Gigi, `proali` **`-w default`** | 811436_1 | 8 | **OOM >128 GB** | — | divergent (diploperennis) chr10; unbounded WFA window OOM'd even at 128 GB |
| `03_align_taxa_to_b73` — Gigi, `proali` **`-w 50000`** | 811459_1 | 8 | **32.7 GB** | **14:20** | COMPLETED, MAF 472 MB. Bounding the WFA window fixed the OOM (~4× less mem). Request was 200 G (16% used) → tighten to ~48 G. |
| `03_align_taxa_to_b73` — TIL18, `proali` **`-w default`** | 811436_0 | 8 | **94.0 GB** | 36:51 | COMPLETED, MAF 442 MB; ref_bp_aligned 148,721,054; 10 blocks |
| `03_align_taxa_to_b73` — TIL18, `proali` **`-w 50000`** | 811571_0 | 8 | **34.8 GB** | 12:16 | COMPLETED, MAF 446 MB; **ref_bp_aligned 148,721,054; 10 blocks — IDENTICAL to unbounded** |
| `create-maf-vcf` (chr10) | — | — | — | — | pending |
| `load-vcf` (chr10) | — | — | — | — | pending |

**AnchorWave memory driver:** proali peak is the base-level WFA of inter-anchor (intergenic) blocks, bounded by
`-w` (window width, default 100000). At default `-w`, unbounded chr10 proali ran 94 GB (TIL18) to >128 GB OOM
(Gigi) — taxon-dependent and unpredictable. Inversions are handled at the anchoring level and are NOT the
memory driver. Do not label taxa collinear/divergent without grounding — size from measured peaks.

**DECIDED: bound all aligns with `-w 50000`.** Same-genome test (TIL18): `-w 50000` gave **identical**
alignment (ref_bp_aligned 148,721,054; 10 blocks) as unbounded, at **~3× less memory (94→35 GB) and ~3× less
time (37→12 min)**. So bounding is quality-neutral on chr10 and strictly cheaper. Universal recipe: `-w 50000`,
`--mem ~48 GB`, uniform for all aligns (no per-taxon OOM risk). Definitive confirmation still pending: compare
gene-range haplotype hashes out of `create-maf-vcf` (bounded vs unbounded MAF).

**Full-genome faidx note:** a full-genome `samtools faidx` (all chromosomes) will read the whole 2.5–3.3 GB
FASTA; peak was 6.25 GB at chr10 extract, so budget ~16 GB for the full-genome index step.
