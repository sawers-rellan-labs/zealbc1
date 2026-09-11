process BINHMM_DOSAGE {
    tag    "cohort"
    label  'rstats'
    publishDir "${params.outdir}/bc2s3_dosage", mode: 'copy'

    cpus   4
    memory '32 GB'
    time   '8h'

    input:
    path counts               // merged allelic_counts50K.tsv (SAMPLE column)
    path samples              // sample,donor map (selects which lines to call)
    path masks                // all <donor>.hd.tsv.gz, staged flat into the work dir

    output:
    path "bc2s3_dosage.tsv.gz"

    script:
    // ONE nilhmm call over the whole cohort: read the merged counts once, restrict each line to its
    // F1 donor's mask sites, then caller="binhmm" (binned Gaussian HMM; dispatches per sample on
    // `name`). No split, no fan-out — binhmm is serial but cheap; memory (hold the merged table) is
    // the real reservation. --masks-dir . picks up the staged *.hd.tsv.gz.
    """
    binhmm_dosage.R \
      --counts ${counts} \
      --samples ${samples} \
      --masks-dir . \
      --design ${params.design} \
      --bin-size ${params.bin_size} \
      --threads ${task.cpus} \
      --out bc2s3_dosage.tsv.gz
    """
}
