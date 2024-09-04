process RSEQC_GENEBODYCOVERAGEPLOT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "perl=5.32.1=7_hd590300_perl5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.32' :
        'quay.io/biocontainers/perl:5.32' }"

    input:
    tuple val(meta), path(rscript)

    output:
    path("geneBodyCoverage.r"), emit: g_rscript
    path("geneBodyCoverage.curves.pdf"), emit: png_curves
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //plot_gene_body_cvg:- plots gene body coverage
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_gene_body_cvg.pl --rfile geneBodyCoverage.r --curves_png geneBodyCoverage.curves.pdf ${rscript}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch geneBodyCoverage.r
    touch geneBodyCoverage.curves.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_distribution.py --version | sed -e "s/read_distribution.py //g")
    END_VERSIONS
    """
}
