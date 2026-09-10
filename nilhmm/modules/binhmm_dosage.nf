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
    // Take the EXISTING BC2S3 counts, keep only this F1's informative (mask) sites — excluding
    // the non-informative ones — then bin and run the binHMM. No alignment, no re-counting.
    // Emission: Gaussian now; beta-binomial on (n_donor, n_total) as the planned upgrade.
    """
    binhmm_dosage.R --counts ${counts} --mask ${mask} --conc ${params.binhmm_conc} \
      --out ${sample}.dosage.tsv.gz
    """
}
