process SAMTOOLS_MPILEUP{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}"
    label 'local'

    input:
    tuple val(meta), path(bam), path(reference)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    def args = task.ext.args ?: ""
    """
    samtools mpileup ${bam} -f ${reference} $args  > all.tsv
    """
}

process SAMTOOLS_MPILEUP_DEAM3{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}"
    label 'local'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    def args = task.ext.args ?: ""
    """
    samtools mpileup ${bam} $args  > all.tsv
    """
}