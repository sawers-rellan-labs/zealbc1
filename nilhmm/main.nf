#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// nilhmm Phase 1, two halves:
//   (1) BUILD THE MASK — align+genotype the 384 BC1 plants, and per F1 donor collect the sites where
//       a BC1 plant is het (the informative sites for that F1): <donor>.hd.tsv.gz.
//   (2) CALL DOSAGE — take the ALREADY-EXISTING merged BC2S3 allelic counts, drop the sites not in
//       each line's F1 mask, and run nilhmm binhmm. BC2S3 is already aligned+counted; half 2 does
//       NOT map or count.

include { INDEX_REF        } from './modules/index_ref'
include { ALIGN            } from './modules/align'
include { GENOTYPE         } from './modules/genotype'
include { QC_INTROGRESSION } from './modules/qc_introgression'
include { BUILD_HD         } from './modules/build_hd'          // per-F1 mask
include { BINHMM_DOSAGE    } from './modules/binhmm_dosage'     // mask -> exclude -> binHMM

workflow {

    // ---- (1) BUILD THE MASK from the BC1 plants -------------------------
    def bc1_sheet = file(params.bc1_samplesheet).splitCsv(header: true)
    donor_of = bc1_sheet.collectEntries { r -> [(r.sample): r.donor] }
    taxon_of = bc1_sheet.collectEntries { r -> [(r.sample): r.taxon] }

    // index the reference once (minibwa + faidx); ALIGN gates on it. .first() = broadcast value.
    ref_ready = INDEX_REF(Channel.value(params.reference)).ready.first()

    bc1_reads = Channel.fromPath(params.bc1_samplesheet)
        .splitCsv(header: true)
        .map { r -> tuple(r.sample, file(r.fastq_1), file(r.fastq_2)) }

    ALIGN(bc1_reads, ref_ready)
    GENOTYPE(ALIGN.out)

    qc_in = GENOTYPE.out.map { sample, vcf, csi ->
        tuple(sample, donor_of[sample], taxon_of[sample], vcf, csi)
    }
    QC_INTROGRESSION(qc_in)

    masks = BUILD_HD(                                           // emits tuple(donor, mask_file)
        QC_INTROGRESSION.out
            .map { sample, donor, vcf, csi, qc -> tuple(donor, vcf, csi, qc) }
            .groupTuple(by: 0)
    )

    // ---- (2) CALL DOSAGE on the existing BC2S3 counts -------------------
    // One binhmm run over the whole cohort: the merged counts file + the sample->donor map + every
    // donor mask. The R script reads the counts once and keeps each line to its F1's mask sites.
    BINHMM_DOSAGE(
        file(params.bc2s3_counts),                             // merged allelic_counts50K.tsv
        file(params.bc2s3_samples),                            // sample,donor map
        masks.map { donor, mask_file -> mask_file }.collect()  // all <donor>.hd.tsv.gz
    )
}
