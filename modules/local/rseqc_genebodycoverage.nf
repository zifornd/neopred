process RSEQC_GENEBODYCOVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.gene_body_coverage.r"), emit: rscript
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //gene_body_cvg_qc:- gives RNA-seq read coverage over the gene body
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    geneBody_coverage.py \\
        -i $bam \\
        -r $bed \\
        -f png \\
        > ${prefix}.gene_body_coverage.r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(geneBody_coverage.py --version | sed -e "s/geneBody_coverage.py //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_body_coverage.r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_distribution.py --version | sed -e "s/read_distribution.py //g")
    END_VERSIONS
    """
}
