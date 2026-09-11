process BINHMM_DOSAGE {
    tag    "${sample}"
    label  'rstats'
    publishDir "${params.outdir}/bc2s3_dosage", mode: 'copy'

    cpus   2
    memory '16 GB'
    time   '4h'

    input:
    tuple val(sample), val(donor), path(counts), path(mask)

    output:
    path "${sample}.dosage.tsv.gz"

    script:
    // Take the EXISTING BC2S3 counts, keep only this F1's informative (mask) sites, then run
    // nilhmm's binned Gaussian HMM (caller="binhmm"). No alignment, no re-counting.
    // Emission: Gaussian now; beta-binomial over BIN counts is a later swap (not bbnil).
    """
    binhmm_dosage.R --counts ${counts} --mask ${mask} --sample ${sample} --donor ${donor} \
      --design ${params.design} --bin-size ${params.bin_size} --out ${sample}.dosage.tsv.gz
    """
}
