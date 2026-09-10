process BUILD_HD {
    tag    "${donor}"
    label  'rstats'
    publishDir "${params.outdir}/hd", mode: 'copy'

    cpus   2
    memory '16 GB'
    time   '4h'

    input:
    tuple val(donor), path(vcfs), path(csis), path(qcs)

    output:
    tuple val(donor), path("${donor}.hd.tsv.gz")

    script:
    // per donor: drop QC-failed plants, union ALT within called teosinte blocks (cross-ear concordance)
    """
    build_hd.R --donor ${donor} --vcfs "${vcfs}" --qc "${qcs}" --out ${donor}.hd.tsv.gz
    """
}
