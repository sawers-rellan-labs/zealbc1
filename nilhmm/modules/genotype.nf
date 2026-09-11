process GENOTYPE {
    tag    "${sample}"
    label  'call'
    publishDir "${params.outdir}/genotype", mode: 'copy'

    cpus   4
    memory '16 GB'
    time   '4h'

    input:
    tuple val(sample), path(cram), path(crai)

    output:
    tuple val(sample), path("${sample}.vcf.gz"), path("${sample}.vcf.gz.csi")

    script:
    // ML genotype only at the bzeaseq biallelic sites (mpileup reads CRAM with -f the reference).
    """
    bcftools mpileup -f ${params.reference} -T ${params.sites} -a AD,DP -Ou ${cram} \
      | bcftools call -m -Oz -o ${sample}.vcf.gz
    bcftools index ${sample}.vcf.gz
    """

    stub:
    """
    : > ${sample}.vcf.gz
    : > ${sample}.vcf.gz.csi
    """
}
