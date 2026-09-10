#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// The point of Branch A is the per-F1 MASK of informative sites. Branch B applies it to the
// ALREADY-EXISTING BC2S3 counts to exclude the non-informative sites per F1, then binHMM.
// BC2S3 is already aligned and counted (rsstu .../BZea/bzeaseq) — Branch B does NOT map or count.

include { ALIGN            } from './modules/align'
include { GENOTYPE         } from './modules/genotype'
include { QC_INTROGRESSION } from './modules/qc_introgression'
include { BUILD_HD         } from './modules/build_hd'          // per-F1 mask
include { BINHMM_DOSAGE    } from './modules/binhmm_dosage'     // mask -> exclude -> binHMM

workflow {

    // ---- Branch A: BC1 individuals -> per-F1 mask -----------------------
    def bc1_sheet = file(params.bc1_samplesheet).splitCsv(header: true)
    donor_of = bc1_sheet.collectEntries { r -> [(r.sample): r.donor] }
    taxon_of = bc1_sheet.collectEntries { r -> [(r.sample): r.taxon] }

    bc1_reads = Channel.fromPath(params.bc1_samplesheet)
        .splitCsv(header: true)
        .map { r -> tuple(r.sample, file(r.fastq_1), file(r.fastq_2)) }

    ALIGN(bc1_reads)
    GENOTYPE(ALIGN.out)

    qc_in = GENOTYPE.out.map { sample, vcf, csi ->
        tuple(sample, donor_of[sample], taxon_of[sample], vcf, csi)
    }
    QC_INTROGRESSION(qc_in)

    mask = BUILD_HD(                                            // emits tuple(donor, mask_file)
        QC_INTROGRESSION.out
            .map { sample, donor, vcf, csi, qc -> tuple(donor, vcf, csi, qc) }
            .groupTuple(by: 0)
    )

    // ---- Branch B: existing BC2S3 counts -> exclude non-informative -> binHMM ----
    // bc2s3_counts.csv : sample,donor,counts   (counts = existing per-line allelic counts)
    bc2_counts = Channel.fromPath(params.bc2s3_counts)
        .splitCsv(header: true)
        .map { r -> tuple(r.donor, r.sample, file(r.counts)) }

    dosage_in = bc2_counts
        .combine(mask, by: 0)                                  // join line to its F1's mask
        .map { donor, sample, counts, mask_file -> tuple(sample, donor, counts, mask_file) }

    BINHMM_DOSAGE(dosage_in)
}
