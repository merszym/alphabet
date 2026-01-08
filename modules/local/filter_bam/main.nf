process FILTER_BAM{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.id}"
    label 'local'

    input:
    tuple val(meta), path(bam), path(positions)

    output:
    tuple val(meta), path("output.deaminated3.bam"), emit: bam

    script:
    def args = task.ext.args ?: ''
    """
    filter_bam.py ${bam} ${positions}
    """
}