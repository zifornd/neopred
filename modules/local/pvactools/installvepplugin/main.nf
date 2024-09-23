process PVACTOOLS_INSTALLVEPPLUGIN {
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/pvactools:4.1.1' :
        'griffithlab/pvactools:4.1.1' }"

    output:
    path("Wildtype.pm")  , emit: results_1
    path("Frameshift.pm")  , emit: results_2

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pvacseq install_vep_plugin ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvacseq: \$(pvacseq --version | sed 's/Python //g')
    END_VERSIONS
    """
}
