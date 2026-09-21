#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// nilhmm Phase 1, two halves:
//   (1) BUILD THE MASK — demux the 32 pooled BC1 libraries into per-plant reads, align + genotype,
//       and per F1 donor collect the sites where a BC1 plant is het: <donor>.hd.tsv.gz.
//   (2) CALL DOSAGE — apply the masks to the ALREADY-EXISTING merged BC2S3 counts and run binhmm.
//       BC2S3 is already aligned + counted; half 2 does NOT map or count.

include { DEMUX            } from './modules/demux'
include { DEMUX_QC         } from './modules/demux_qc'
include { INDEX_REF        } from './modules/index_ref'
include { ALIGN            } from './modules/align'
include { GENOTYPE         } from './modules/genotype'
include { MOSDEPTH         } from './modules/mosdepth'          // per-plant mapped coverage
include { MULTIQC          } from './modules/multiqc'           // aggregate coverage + demux QC
include { QC_INTROGRESSION } from './modules/qc_introgression'
include { BUILD_HD         } from './modules/build_hd'          // per-F1 mask
include { BINHMM_DOSAGE    } from './modules/binhmm_dosage'     // mask -> exclude -> binHMM

// (pool, [R1s], [R2s]) from DEMUX -> (Sample_Id, R1, R2), pairing by Sample_Id. Shared by both entries.
def flatten_reads(ch) {
    ch.flatMap { pool, r1s, r2s ->
        def l1 = (r1s instanceof List) ? r1s : [r1s]
        def l2 = (r2s instanceof List) ? r2s : [r2s]
        def m2 = l2.collectEntries { f -> [(f.name.replaceFirst(/_R2\.fq\.gz$/, '')): f] }
        l1.collect { f1 -> def s = f1.name.replaceFirst(/_R1\.fq\.gz$/, ''); tuple(s, f1, m2[s]) }
    }
}

// comma/space list param -> Set, or null (= keep all)
def keep_set(p) { def k = (p ?: '').toString().trim(); k ? (k.split(/[,\s]+/) as List) : null }

workflow {

    // ---- (1) BUILD THE MASK from the BC1 plants -------------------------
    // sample -> donor / taxon from the well map (Sample_Id is the key through the whole pipeline).
    def well = file(params.bc1_well_map).splitCsv(header: true)
    donor_of = well.collectEntries { r -> [(r.Sample_Id): r.donor] }
    taxon_of = well.collectEntries { r -> [(r.Sample_Id): r.taxon] }

    bc_fasta = file(params.bc_fasta)
    well_map = file(params.bc1_well_map)

    // index the reference once (minibwa + faidx); ALIGN gates on it. .first() = broadcast value.
    ref_ready = INDEX_REF(Channel.value(params.reference)).ready.first()

    // per-pool inputs: (pool, [R1 lanes], [R2 lanes]). --pools '1B,2C' restricts pools (Gate 1/testing).
    def pool_set = keep_set(params.pools)
    libs = Channel.fromPath(params.bc1_libraries)
        .splitCsv(header: true)
        .filter { r -> pool_set == null || pool_set.contains(r.pool) }
        .map { r -> tuple(r.pool,
                          files("${params.bc1_rawdata}/${r.raw_dir}/*_1.fq.gz"),
                          files("${params.bc1_rawdata}/${r.raw_dir}/*_2.fq.gz")) }

    DEMUX(libs, bc_fasta, well_map)

    // flatten each pool's per-plant FASTQs into (Sample_Id, R1, R2), pairing by Sample_Id
    reads = flatten_reads(DEMUX.out.reads)

    ALIGN(reads, ref_ready)
    GENOTYPE(ALIGN.out)
    MOSDEPTH(ALIGN.out)                                          // per-plant mapped coverage

    qc_in = GENOTYPE.out.map { sample, vcf, csi ->
        tuple(sample, donor_of[sample], taxon_of[sample], vcf, csi)
    }
    QC_INTROGRESSION(qc_in)

    masks = BUILD_HD(                                           // emits tuple(donor, mask_file)
        QC_INTROGRESSION.out
            .map { sample, donor, vcf, csi, qc -> tuple(donor, vcf, csi, qc) }
            .groupTuple(by: 0)
    )

    // demux balance QC (flag, don't block) — aggregates every pool's cutadapt json
    DEMUX_QC(DEMUX.out.json.collect(), well_map)

    // one per-sample coverage/QC report across all plants: mosdepth (mapped coverage) + cutadapt (demux)
    MULTIQC(MOSDEPTH.out.mix(DEMUX.out.json).collect())

    // ---- (2) CALL DOSAGE on the existing BC2S3 counts -------------------
    // One binhmm run over the whole cohort: merged counts + sample->donor map + every donor mask.
    // Independent of the mask half; runs only once its inputs exist (the sample->donor map is built
    // from the pedigree separately), so Gate 0 can validate the mask half on its own.
    if (file(params.bc2s3_counts).exists() && file(params.bc2s3_samples).exists()) {
        BINHMM_DOSAGE(
            file(params.bc2s3_counts),                             // merged allelic_counts50K.tsv
            file(params.bc2s3_samples),                            // sample,donor map
            masks.map { donor, mask_file -> mask_file }.collect()  // all <donor>.hd.tsv.gz
        )
    }
}

// ---- BC2S3 batch 2 (1.2x, row-pooled with the same 12 inline barcodes) ---------------------
// `nextflow run main.nf -entry demux_bc2s3_batch2`: DEMUX every row library (32 = V21A..V24H), then
// ALIGN + MOSDEPTH only the samples in --samples (comma list of Sample_Id = P<Plot_id>; '' = all 384).
// Reuses DEMUX/ALIGN/MOSDEPTH unchanged: the batch-2 well map has the same first four columns
// (pool,column,barcode,Sample_Id) + taxon, so DEMUX's rename and demux_qc.R work as for BC1.
// No GENOTYPE/QC/BUILD_HD here: those are BC1 (H_d) steps; BC2S3 lines go to the poolseq/PHG route.
workflow demux_bc2s3_batch2 {
    bc_fasta = file(params.bc_fasta)
    well_map = file(params.bc2s3_well_map)

    ref_ready = INDEX_REF(Channel.value(params.reference)).ready.first()

    def pool_set   = keep_set(params.pools)
    def sample_set = keep_set(params.samples)
    libs = Channel.fromPath(params.bc2s3_libraries)
        .splitCsv(header: true)
        .filter { r -> pool_set == null || pool_set.contains(r.pool) }
        .map { r -> tuple(r.pool,
                          files("${params.bc2s3_rawdata}/${r.raw_dir}/*_1.fq.gz"),
                          files("${params.bc2s3_rawdata}/${r.raw_dir}/*_2.fq.gz")) }

    DEMUX(libs, bc_fasta, well_map)

    reads = flatten_reads(DEMUX.out.reads)
        .filter { s, r1, r2 -> sample_set == null || sample_set.contains(s) }

    ALIGN(reads, ref_ready)
    MOSDEPTH(ALIGN.out)                                          // per-line depth = the coverage label

    DEMUX_QC(DEMUX.out.json.collect(), well_map)
    MULTIQC(MOSDEPTH.out.mix(DEMUX.out.json).collect())
}
