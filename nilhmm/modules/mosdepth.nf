process MOSDEPTH {
    tag    "${sample}"
    label  'qc'
    publishDir "${params.outdir}/mosdepth", mode: 'copy'

    cpus   4
    memory '4 GB'
    time   '1h'

    input:
    tuple val(sample), path(cram), path(crai)

    output:
    path "${sample}.mosdepth.*"

    script:
    // Mapped genome coverage per plant (fast, light). --fasta decodes the CRAM. If a SNP50K sites BED
    // is available, --by adds a second track = depth at the informative sites (what matters for het calling).
    """
    BY=""
    if [ -n "${params.sites_bed}" ] && [ -s "${params.sites_bed}" ]; then BY="--by ${params.sites_bed}"; fi
    mosdepth -t ${task.cpus} --no-per-base --fasta ${params.reference} \$BY ${sample} ${cram}
    """

    stub:
    """
    : > ${sample}.mosdepth.summary.txt
    : > ${sample}.mosdepth.global.dist.txt
    """
}
