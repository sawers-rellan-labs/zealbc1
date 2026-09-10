process GENOTYPE {
    tag    "${sample}"
    label  'call'
    publishDir "${params.outdir}/genotype", mode: 'copy'

    cpus   4
    memory '16 GB'
    time   '4h'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.vcf.gz"), path("${sample}.vcf.gz.csi")

    script:
    // genotype only at the bzeaseq biallelic sites (targets/streaming, AD+DP for the count model)
    """
    bcftools mpileup -f ${params.reference} -T ${params.sites} -a AD,DP -Ou ${bam} \
      | bcftools call -m -Oz -o ${sample}.vcf.gz
    bcftools index ${sample}.vcf.gz
    """
}
