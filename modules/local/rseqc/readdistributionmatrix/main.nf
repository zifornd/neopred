process RSEQC_READDISTRIBUTIONMATRIX {
    label 'process_medium'

    conda (params.enable_conda ? "perl=5.32.1=7_hd590300_perl5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.32' :
        'quay.io/biocontainers/perl:5.32' }"

    input:
    path(txt)

    output:
    path("read_distrib.matrix.tab"), emit: matrix_tab
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_distrib_matrix.pl -f ${txt} 1>read_distrib.matrix.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch read_distrib.matrix.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
