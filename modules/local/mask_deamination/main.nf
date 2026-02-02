process MASK_DEAMINATION{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.id}"
    label 'local'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("masked_${bam}"), path('count.txt'), emit: bam

    script:
    """
    mask_qual_scores.py ${bam} > count.txt
    """
}