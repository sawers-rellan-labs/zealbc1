process ALIGN {
    tag    "${sample}"
    label  'align'
    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram*'

    cpus   8
    memory '24 GB'
    time   '6h'

    input:
    tuple val(sample), path(r1), path(r2)
    val ready                          // gate: reference is minibwa-indexed + faidx'd

    output:
    tuple val(sample), path("${sample}.cram"), path("${sample}.cram.crai")

    script:
    // CRAM is the durable per-plant product (~2-3 GB @ ~9x); the demuxed FASTQs are disposable (work/).
    // reference read by absolute path (must be minibwa-indexed once by INDEX_REF).
    """
    set -euo pipefail
    minibwa map -t ${task.cpus} ${params.reference} ${r1} ${r2} \
      | samtools sort -@ 2 -o aln.bam
    samtools view -T ${params.reference} -C -F 0x904 -q ${params.mapq} -o ${sample}.cram aln.bam
    samtools index ${sample}.cram
    rm -f aln.bam
    """

    stub:
    """
    : > ${sample}.cram
    : > ${sample}.cram.crai
    """
}
