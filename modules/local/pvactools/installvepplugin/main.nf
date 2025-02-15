process PVACTOOLS_INSTALLVEPPLUGIN {
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/pvactools:4.1.1' :
        'docker.io/griffithlab/pvactools:4.1.1' }"

    output:
    path "*.pm"            , emit: results
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pvacseq install_vep_plugin ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
