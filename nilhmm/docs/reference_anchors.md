# BZea cohort-QC reference anchors (one per taxon)

`nilhmm/bin/cohort_qc.R` places one **reference "pole" sample per taxon** in each PCA (`fixed_anchor`). The taxon is **forced** from this map — the SNP50K panel metadata mislabels the TIL lines (mexicana TILs shown as `Zv`), and the diploperennis Gigi/Momo accession is not in the panel.

## Chosen anchors

| taxon | anchor sample | PanAnd assembly match | het¹ | basis |
|---|---|---|---|---|
| parviglumis | `TIL11` | Zv-TIL11 | — | PanAnd reference, present in panel |
| mexicana | `TIL25` | (Zx-TIL25) | — | Zx-**TIL18 absent** from panel; TIL25 is the mexicana TIL present (metadata mislabels it Zv) |
| luxurians | `RIL003` | Zl-RIL003 | — | PanAnd reference, present |
| huehuetenangensis | `RIMH001` | Zh-RIMHU001 | — | PanAnd reference, present |
| diploperennis | `Ame2317` | (none) | **0.052** | Gigi/Momo (PI 462368) **absent** from panel; Ame2317 = least-heterozygous of the 5 diplo refs |

¹ `Proportion_Heterozygous` from `schnable2023_metadata.tab`. All reference samples are **0-missing at SNP50K**, so missingness cannot pick a representative; heterozygosity breaks the diploperennis tie.

## diploperennis candidates (why Ame2317, and the alternative)

| sample | accession | locality | lat, long | het |
|---|---|---|---|---|
| **Ame2317** ✔ | Ames 2317 | — | — | 0.052 |
| Ame21884 | Ames 21884 (2265A) | — | 19.59, −104.26 | 0.088 |
| 5D3 | CIMMYTMA 9476 | Las Joyas | 19.37, −104.12 | 0.097 |
| **5D7** | CIMMYTMA 10003 | **San Miguel** | 19.35, −104.15 | 0.091 |
| RRA-14 | — (no accession) | — | — | 0.077 |

- **Ame2317** chosen for the **cleanest pole** (least within-individual heterozygosity; Zd is obligate-outcrossing, so every individual is heterozygous).
- **5D7 (CIMMYTMA 10003, "San Miguel")** is the geographic alternative: it is from **Cerro de San Miguel, Sierra de Manantlán, Jalisco** — the *type locality of PI 462368*, i.e. the Gigi/Momo source accession. Prefer it if geographic fidelity to the PanAnd Zd lineage matters more than low het. To switch: set `diploperennis = "5D7"` in `cohort_qc.R`'s `fixed_anchor`.
- **RRA-14** has no genebank accession/locality (all `-` in the passport) — likely a **Rellán-Álvarez / Hufford field collection** (`RRA` ≈ Rubén Rellán-Álvarez; a Bajío collection with Matt Hufford). **Taxon is genetically confirmed Zd** — its PCA coordinates in `schnable2023_metadata.tab` (PC1 ≈ −0.073, PC2 ≈ 0.033) sit in the diploperennis cluster with the other diplo refs. **Locality unconfirmed:** if the "Bajío" recollection is right, note that diploperennis is normally Sierra de Manantlán / western-Jalisco montane, not the Bajío (annual mexicana/parviglumis country) — pin the collection site down with Ruben. Not used as an anchor.

## Data sources

- **SNP50K genotypes + sample classification:** `/rsstu/users/r/rrellan/BZea/bzeaseq/50K/results/joint/bzea_50K_cohort_ref.vcf.gz` and `bzea_50K_cohort_ref_metadata.csv` (`is_reference`, `is_B73`, `maizegdb_prefix`, `taxa_label`).
- **Per-sample heterozygosity:** `…/bzeaseq/schnable2023/schnable2023_metadata.tab` (`Proportion_Heterozygous`); reference samples are Chen et al. 2022 wild relatives.
- **Passport (accession IDs, localities):** `…/bzeaseq/chen2022_passport.tab`. PI 462368 appears there as `PI46-2368` (Zea diploperennis) — the Gigi/Momo source accession — but is **not** among the panel's diploperennis reference samples.

## References — Zd Gigi/Momo lineage (PI 462368)

Both `Zd-Gigi` and `Zd-Momo` (PanAnd) are separate clonal individuals of USDA NPGS **PI 462368** (Cerro de San Miguel, Sierra de Manantlán, Jalisco, Mexico):
- PanAnd genome project — *Genetics* 227(1):iyae036 <https://academic.oup.com/genetics/article/227/1/iyae036/7641224>
- PI 462368 germplasm context — PMC8421791 <https://pmc.ncbi.nlm.nih.gov/articles/PMC8421791/>
- *Zea diploperennis* description — Iltis, Doebley, Guzmán & Pazy, *Science* 203(4376):186 (1979) <https://www.science.org/doi/10.1126/science.203.4376.186>
- Zd-Momo assembly download — <https://download.maizegdb.org/Zd-Momo-REFERENCE-PanAnd-1.0/>

## Caveat

The panel's `maizegdb_prefix` is **unreliable for the `TIL*` lines** (mexicana TILs labeled `Zv`). Anything keying on that prefix for TILs must override it, as `fixed_anchor` does by forcing taxon from the map.
