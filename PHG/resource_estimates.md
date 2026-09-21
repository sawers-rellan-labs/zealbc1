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
| DB stage 1 `initdb→prepare→agc→create-ranges→create-ref-vcf` | 812418 | 8 | 2.8 GB | 58 s | COMPLETED. agc 96 MB; **4,685 chr10 ranges** (gene+intergenic, pad 500); B73 ref haplotypes loaded. prepare-assemblies appends ` sampleName=X`, keeps orig contig (B73 stays `chr10`). |
| DB stage 2 `create-maf-vcf (2 taxa) → load-vcf` | 812482 | 8 | **12.8 GB** | **5:13** | COMPLETED. TIL18+Gigi h.vcf/g.vcf created + loaded. **Graph now = 3 samples (B73+TIL18+Gigi).** MAF files named `<sample>.maf`; query contigs match AGC → sequence resolves. `create-maf-vcf` doesn't mkdir its `-o` (must pre-create). |

**PILOT VALIDATED END-TO-END (chr10):** subset → AnchorWave `proali` → prepare-assemblies → agc-compress →
create-ranges → create-ref-vcf → create-maf-vcf → load-vcf → 3-founder graph. The full phg build recipe works.

**AnchorWave memory driver:** proali peak is the base-level WFA of inter-anchor (intergenic) blocks, bounded by
`-w` (window width, default 100000). At default `-w`, unbounded chr10 proali ran 94 GB (TIL18) to >128 GB OOM
(Gigi) — taxon-dependent and unpredictable. Inversions are handled at the anchoring level and are NOT the
memory driver. Do not label taxa collinear/divergent without grounding — size from measured peaks.

**DECIDED: bound all aligns with `-w 50000 -t 8` (`--mem 48 G`).** Confirmed quality-neutral on BOTH a
less-divergent (TIL18) and the most-divergent available (Gigi) taxon — same MAF coverage to the base pair:

| taxon | recipe | ref_bp aligned | mem | wall |
|---|---|---|---|---|
| TIL18 | `-w 50000 -t8`  | 148,721,054 (10 blocks) | 35 GB | 12 min |
| TIL18 | `-w 100000 -t8` | 148,721,054 (10 blocks) — identical | 94 GB | 37 min |
| Gigi  | `-w 50000 -t8`  | 139,125,172 (9 blocks) | 33 GB | 14 min |
| Gigi  | `-w 100000 -t4` | 139,125,172 (9 blocks) — identical | 67 GB | 54 min |
| Gigi  | `-w 100000 -t8 -M100` | 139,125,172 (9 blocks) — identical | 108 GB | 23 min |

Gigi's 91.3% ref coverage is **real divergence, not window truncation** (identical at full `-w`). So `-w 50000`
loses zero alignment and is cheapest on both axes — it wins outright.

**`-M` / anchorwave 1.3.1 proved unnecessary for the align.** `-M` works (feasible `-M100 -t8` held ~108 GB and
completed; `-M/-t` must be ≥ ~10 GB/worker, and `-M` is predictive so give `--mem` headroom), but since `-w 50000`
is quality-neutral and far cheaper, we don't need it. `anchorwave13` (1.3.1) stays as the align env (newer, `-M`
available if ever needed); the recipe is plain `proali -w 50000 -t 8`.

**Memory model = `-t × per-thread-WFA(-w)`, not just `-w`.** proali parallelizes inter-anchor region alignment
across threads, so peak ≈ (worst-case per-thread WFA memory, set by `-w`) × `-t`. Evidence: Gigi `-t 8` unbounded
OOM'd in 6:37 (8 threads spiking on big divergent gaps simultaneously); bounded `-t 8` ≈ 33–35 GB (~4–5 GB/thread).
Implications for the FULL run:
- chr10's ~35 GB / 48 GB does **NOT** extrapolate to full genome — the whole genome has more/larger inter-anchor
  gaps, so per-thread memory is higher and `×t` amplifies it. **Measure full-genome bounded proali before sizing
  the 100-align allocation.**
- `-t` is a memory/speed dial: if full-genome bounded `-t 8` runs hot, `-t 4` roughly halves the peak (slower).
  Keep `-w 50000` (quality-neutral) and tune `-t`/`--mem` from the full-genome measurement.

**Full-genome faidx note:** a full-genome `samtools faidx` (all chromosomes) will read the whole 2.5–3.3 GB
FASTA; peak was 6.25 GB at chr10 extract, so budget ~16 GB for the full-genome index step.
