process MULTIQC {
    tag    "cohort"
    label  'qc'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    cpus   1
    memory '4 GB'
    time   '30m'

    input:
    path qc_files          // mosdepth summaries + cutadapt demux jsons (staged flat)

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    // One per-sample report across all 384: MultiQC's mosdepth module -> coverage; cutadapt module -> demux.
    """
    multiqc -f -o . .
    """

    stub:
    """
    : > multiqc_report.html
    mkdir -p multiqc_data
    """
}
