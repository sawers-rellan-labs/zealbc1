process QC_INTROGRESSION {
    tag    "${sample}"
    label  'rstats'
    publishDir "${params.outdir}/qc", mode: 'copy'

    cpus   2
    memory '16 GB'
    time   '2h'

    input:
    tuple val(sample), val(donor), val(taxon), path(vcf), path(csi)

    output:
    tuple val(sample), val(donor), path(vcf), path(csi), path("${sample}.qc.tsv")

    script:
    // contamination QC: taxon match, hom-teosinte VAF tail, segment structure -> pass/fail flag
    """
    qc_introgression.R --vcf ${vcf} --sample ${sample} --donor ${donor} --taxon ${taxon} \
      --out ${sample}.qc.tsv
    """
}
