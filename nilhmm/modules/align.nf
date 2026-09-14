process ALIGN {
    tag    "${sample}"
    label  'align'
    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram*'

    cpus   8
    // OOM auto-escalation (the LSF "increase RAM on retry" idiom): start at 24 GB (covers the ~18.6 GB
    // peak observed for S_1B_10), and on an OOM-kill Nextflow resubmits the SAME task at 48, then 72 GB
    // (maxRetries=2 in nextflow.config). Static 16 GB was the bug — retries re-OOM'd at the same request.
    memory { 24.GB * task.attempt }
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
    # -x sr = short-read preset (correct for Illumina 150bp; -x adap invokes long-read DP → heavier).
    # No -K batch cap: tested -K 100m,250m and it slowed mapping ~2x WITHOUT bounding memory. Peak comes
    # from trace.txt on this run (only clean prior datum: 18.6 GB for S_1B_10 under -x adap).
    minibwa map -x sr -t ${task.cpus} ${params.reference} ${r1} ${r2} \
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
