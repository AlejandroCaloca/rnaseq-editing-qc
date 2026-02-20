/*
========================================================================================
    LOCAL MODULE: HTSEQ_COUNT
========================================================================================
*/

process HTSEQ_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::htseq=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py310ha14a713_0' :
        'biocontainers/htseq:2.0.2--py310ha14a713_0' }"

    input:
    tuple val(meta), path(bam)
    path gtf

    output:
    tuple val(meta), path("*.counts.txt"), emit: counts
    path "versions.yml",                   emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stranded = params.htseq_stranded ?: 'reverse'
    def mode = params.htseq_count_mode ?: 'union'
    """
    # Sort BAM by name for HTSeq
    samtools sort -n -@ ${task.cpus} ${bam} -o ${prefix}.namesorted.bam

    # Run HTSeq-count
    htseq-count \\
        --format=bam \\
        --stranded=${stranded} \\
        --mode=${mode} \\
        --additional-attr=gene_name \\
        --counts_output=${prefix}.counts.txt \\
        ${prefix}.namesorted.bam \\
        ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(python -c "import HTSeq; print(HTSeq.__version__)")
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """
}
