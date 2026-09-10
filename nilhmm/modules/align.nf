process ALIGN {
    tag    "${sample}"
    label  'align'
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.bam*'

    cpus   8
    memory '24 GB'
    time   '6h'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    // reference is read by absolute path (pre-indexed on /rsstu) rather than staged.
    """
    bwa-mem2 mem -t ${task.cpus} ${params.reference} ${r1} ${r2} \
      | samtools sort -@ 2 -o aln.bam
    samtools view -b -F 0x904 -q ${params.mapq} aln.bam > ${sample}.bam
    samtools index ${sample}.bam
    """
}
