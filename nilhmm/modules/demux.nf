process DEMUX {
    tag    "${pool}"
    label  'align'          // cutadapt lives in envs/assembly.yml
    publishDir { "${params.outdir}/demux/${pool}" }, mode: 'copy', pattern: '*.{json,unknown_R*.fq.gz}'

    cpus   8
    memory '16 GB'     // measured MaxRSS ~8 GB (cutadapt -j8 -Z); headroom for the full run
    time   '6h'

    input:
    tuple val(pool), path(r1s), path(r2s)
    path bc_fasta
    path well_map

    output:
    tuple val(pool), path("S_*_R1.fq.gz"), path("S_*_R2.fq.gz"), emit: reads
    path "${pool}.cutadapt.json",                                emit: json
    path "${pool}.unknown_R{1,2}.fq.gz",                         emit: unknown

    script:
    // Inline-inline exact-match demux. -e 0 --no-indels is NON-NEGOTIABLE: barcode min Hamming
    // distance is 2, so any mismatch misassigns between adjacent wells; the symmetric barcode
    // (matched independently on R1 and R2) is the redundancy. -Z = fast gzip (the real cost).
    // Untrimmed reads are kept (not discarded) so DEMUX_QC can count them. Lanes are concatenated
    // (sorted identically for R1/R2) so pair order is preserved. Per-column outputs are renamed to
    // Sample_Id via the well map (this pool only).
    """
    set -euo pipefail
    R1LANES=\$(ls *_L*_1.fq.gz | sort); R2LANES=\$(ls *_L*_2.fq.gz | sort)

    # Shared cutadapt call; inputs appended below. -e 0 --no-indels is non-negotiable (barcode min
    # Hamming distance 2); -Z fast gzip; symmetric ^file: barcodes matched on R1 and R2 independently.
    CA="cutadapt -j ${task.cpus} -Z -e 0 --no-indels --pair-adapters --action=trim \
      -g ^file:${bc_fasta} -G ^file:${bc_fasta} \
      --untrimmed-output ${pool}.unknown_R1.fq.gz --untrimmed-paired-output ${pool}.unknown_R2.fq.gz \
      --json ${pool}.cutadapt.json \
      -o ${pool}_{name}_R1.fq.gz -p ${pool}_{name}_R2.fq.gz"

    if [ "${params.subsample}" -gt 0 ]; then
        # Gate 1: first N read pairs to small files (head stops the stream early).
        set +o pipefail
        zcat \$R1LANES | head -n \$(( ${params.subsample} * 4 )) | gzip > R1.sub.fq.gz
        zcat \$R2LANES | head -n \$(( ${params.subsample} * 4 )) | gzip > R2.sub.fq.gz
        set -o pipefail
        \$CA R1.sub.fq.gz R2.sub.fq.gz
    else
        # STREAM the lanes straight into cutadapt via process substitution -- no ~230 GB concatenated
        # intermediate written to (and re-read from) the contended /rsstu partition. Measured on pool
        # 1B: the old cat-first path did ~490 GB read + ~494 GB write; this removes ~half of both.
        \$CA <(zcat \$R1LANES) <(zcat \$R2LANES)
    fi

    awk -F, -v p="${pool}" 'NR>1 && \$1==p {print \$2"\\t"\$4}' ${well_map} | while IFS=\$'\\t' read -r col sid; do
      mv "${pool}_\${col}_R1.fq.gz" "\${sid}_R1.fq.gz"
      mv "${pool}_\${col}_R2.fq.gz" "\${sid}_R2.fq.gz"
    done
    rm -f R1.sub.fq.gz R2.sub.fq.gz    # subsample temp (streaming path writes no concat)
    """

    stub:
    """
    awk -F, -v p="${pool}" 'NR>1 && \$1==p {print \$4}' ${well_map} | while read -r sid; do
      : > "\${sid}_R1.fq.gz"; : > "\${sid}_R2.fq.gz"
    done
    : > ${pool}.cutadapt.json
    : > ${pool}.unknown_R1.fq.gz; : > ${pool}.unknown_R2.fq.gz
    """
}
