process INDEX_REF {
    tag    "index_ref"
    label  'align'

    cpus   8
    memory '64 GB'     // measured: peak RSS 39 GB for B73 v5 (2.3 Gb); 64 GB gives headroom
    time   '30m'       // measured: minibwa index B73 v5 = 3m56s (8 cpu). 30m is ample headroom

    input:
    val reference

    output:
    val reference, emit: ready     // emitted only after indexing succeeds -> gates ALIGN

    script:
    // Writes index files next to the reference (absolute path on /rsstu). Skips if already present.
    """
    set -euo pipefail
    if [ ! -s "${reference}.mbw" ]; then
        minibwa index -t ${task.cpus} "${reference}"        # -> <ref>.l2b, <ref>.mbw
    else
        echo "minibwa index present; skipping"
    fi
    if [ ! -s "${reference}.fai" ]; then
        samtools faidx "${reference}"                        # -> <ref>.fai (for bcftools mpileup)
    else
        echo "faidx present; skipping"
    fi
    """

    stub:
    """
    echo "stub INDEX_REF ${reference}"
    """
}
