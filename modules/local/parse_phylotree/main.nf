process SUMMARIZE_PHYLOTREE{
    container (workflow.containerEngine ? "merszym/anytree:nextflow" : null)
    tag "${meta.id}"
    label 'local'

    input:
    tuple val(meta), path(pileup), path(xml)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    def args = task.ext.args
    """
    main.py ${xml} ${pileup} ${meta.id} $args
    """
}