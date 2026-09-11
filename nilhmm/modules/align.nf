process ALIGN {
    tag    "${sample}"
    label  'align'
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.bam*'

    cpus   8
    memory '24 GB'
    time   '6h'

    input:
    tuple val(sample), path(r1), path(r2)
    val ready                          // gate: reference is minibwa-indexed + faidx'd

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    // reference read by absolute path (rather than staged); must be minibwa-indexed once:
    //   minibwa index ${params.reference}   -> <ref>.l2b, <ref>.mbw   (18N RAM; add -l for low-mem)
    """
    minibwa map -t ${task.cpus} ${params.reference} ${r1} ${r2} \
      | samtools sort -@ 2 -o aln.bam
    samtools view -b -F 0x904 -q ${params.mapq} aln.bam > ${sample}.bam
    samtools index ${sample}.bam
    """
}
