process DEMUX_QC {
    tag    "cohort"
    label  'rstats'
    publishDir "${params.outdir}/demux_qc", mode: 'copy'

    cpus   1
    memory '4 GB'
    time   '30m'

    input:
    path jsons          // all <pool>.cutadapt.json
    path well_map

    output:
    path "demux_qc.tsv"
    path "demux_balance.png", optional: true

    script:
    // Flag, don't block: per-pool/per-column read counts + untrimmed% from the cutadapt jsons;
    // flags wells far below their pool median (pipetting) and pools with high untrimmed%.
    """
    Rscript "${projectDir}/bin/demux_qc.R" --jsons "${jsons}" --well-map ${well_map} --out demux_qc.tsv
    """

    stub:
    """
    : > demux_qc.tsv
    """
}
