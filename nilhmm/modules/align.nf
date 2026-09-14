process ALIGN {
    tag    "${sample}"
    label  'align'
    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram*'

    cpus   8
    // Memory = fixed floor (reference index ~5 GB, loaded once) + a VARIABLE part from alignment/chaining
    // buffers that scales with how hard reads are to align (divergence, multi-mapping), NOT with file size.
    // So it's neither linear-in-size (reads are streamed) nor flat (S_1B_10 peaked 18.6 GB, S_1B_4 >24 GB
    // in the same -x adap run) — a floor plus a spiky, content-driven term. BC1 is ~half divergent teosinte,
    // which inflates that term unevenly across samples. Start 32 GB (above the observed range), escalate to
    // 64/96 on an OOM-kill for spiky outliers (maxRetries=2). trace.txt gives the real per-sample peaks.
    memory { 32.GB * task.attempt }
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
