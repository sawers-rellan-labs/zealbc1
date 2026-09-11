process DEMUX {
    tag    "${pool}"
    label  'align'          // cutadapt lives in envs/assembly.yml
    publishDir "${params.outdir}/demux/${pool}", mode: 'copy', pattern: '*.{json,unknown_R*.fq.gz}'

    cpus   8
    memory '8 GB'
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
    cat \$(ls *_L*_1.fq.gz | sort) > R1.fq.gz
    cat \$(ls *_L*_2.fq.gz | sort) > R2.fq.gz

    cutadapt -j ${task.cpus} -Z -e 0 --no-indels --pair-adapters --action=trim \
      -g ^file:${bc_fasta} -G ^file:${bc_fasta} \
      --untrimmed-output ${pool}.unknown_R1.fq.gz --untrimmed-paired-output ${pool}.unknown_R2.fq.gz \
      --json ${pool}.cutadapt.json \
      -o ${pool}_{name}_R1.fq.gz -p ${pool}_{name}_R2.fq.gz \
      R1.fq.gz R2.fq.gz

    awk -F, -v p="${pool}" 'NR>1 && \$1==p {print \$2"\\t"\$4}' ${well_map} | while IFS=\$'\\t' read -r col sid; do
      mv "${pool}_\${col}_R1.fq.gz" "\${sid}_R1.fq.gz"
      mv "${pool}_\${col}_R2.fq.gz" "\${sid}_R2.fq.gz"
    done
    rm -f R1.fq.gz R2.fq.gz
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
